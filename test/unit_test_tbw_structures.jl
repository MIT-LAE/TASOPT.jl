# Logging is not a test dependency, so use the copy in Base
const CoreLogging = Base.CoreLogging

# ============================================================================
# Shared helpers: reading aircraft from TOML text, and sizing aircraft once and sharing them
# ============================================================================

"Default input text with the wing configuration lines removed, so `edit` controls them."
default_toml_without_wing_config() = replace(
    read(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"), String),
    r"^[ \t]*(wing_config|has_strut)[ \t]*=.*\n"m => "")

"Insert `line` right below the `[table]` header in `txt`."
insert_below_header(txt, table, line) = replace(txt, "[$table]" => "[$table]\n    $line"; count = 1)

"Default input text with `wing_config = \"config\"` set in [Options]."
toml_with_wing_config(config) =
    insert_below_header(default_toml_without_wing_config(), "Options", "wing_config = \"$config\"")

"Read an aircraft from TOML text; return the aircraft and the log records emitted while reading."
function read_aircraft_from_text(txt)
    logger = Test.TestLogger(min_level = CoreLogging.Info)
    ac = mktempdir() do dir
        path = joinpath(dir, "wing_config_test.toml")
        write(path, txt)
        CoreLogging.with_logger(logger) do
            read_aircraft_model(path)
        end
    end
    return ac, logger.logs
end

read_aircraft_with_wing_config(config) = first(read_aircraft_from_text(toml_with_wing_config(config)))

"Parsed default input, with `edit!` applied to the Dict; return the aircraft and the log records."
function read_aircraft_from_default_dict(edit!)
    d = TASOPT.TOML.parsefile(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"))
    edit!(d)
    return read_aircraft_from_text(sprint(TASOPT.TOML.print, d))
end

const sized_aircraft_cache = Dict{Any,Any}()

"""
Aircraft built by `make()` and sized on the first request for `key`, then shared by all later requests.
Tests that read it must not modify it: use `deepcopy` first.
"""
function sized_once(make, key; kwargs...)
    return get!(sized_aircraft_cache, key) do
        CoreLogging.with_logger(CoreLogging.NullLogger()) do
            ac = make()
            size_aircraft!(ac; printiter = false, kwargs...)
            ac
        end
    end
end

sized_default_aircraft() = sized_once(load_default_model, :default; iter = 50)

"""
Box of the strut recomputed with `size_member_box` from the loads, geometry and inputs recorded on the sized
aircraft `ac`: the strut chord is the chord ratio times the wing chord at the strut attachment, the webs carry the
recorded shear load (jury load and own weight) and the caps the recorded own-weight bending moment.
"""
function recomputed_strut_box(ac)
    wing, strut = ac.wing, ac.wing.strut
    has_jury = ac.options.opt_wing_config == WingConfig.TrussBraced
    return TASOPT.structures.size_member_box(strut.caps, strut.webs;
        chord = strut.chord_to_wing_chord * wing.layout.root_chord * wing.inboard.λ,
        height_to_chord = strut.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
        web_to_box_height = strut.web_to_box_height,
        tension = strut.axial_force, compression = strut.axial_force_compression,
        buckling_length = TASOPT.structures.strut_buckling_length(strut.length, has_jury ? wing.jury.eta_strut : nothing),
        K = strut.end_fixity_K, shear_load = strut.shear_load, bending_moment = strut.bending_moment,
        sigfac = ac.parg[igsigfac],
        min_cap_thickness = strut.min_cap_thickness, min_web_thickness = strut.min_web_thickness)
end

"""
Box of the jury strut recomputed with `size_member_box`: a pinned member of the chord ratio times the strut chord,
carrying its recorded compression over its length, without tension or shear.
"""
function recomputed_jury_box(ac)
    strut, jury = ac.wing.strut, ac.wing.jury
    return TASOPT.structures.size_member_box(jury.caps, jury.webs;
        chord = jury.chord_to_strut_chord * strut.chord,
        height_to_chord = strut.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
        web_to_box_height = strut.web_to_box_height,
        tension = 0.0, compression = jury.axial_force, buckling_length = jury.length, K = 1.0,
        shear_load = 0.0, sigfac = ac.parg[igsigfac],
        min_cap_thickness = jury.min_cap_thickness, min_web_thickness = jury.min_web_thickness)
end

"""
Shear `V` [N] and midspan bending moment `M` [N m] of the own weight of the strut of the sized aircraft `ac`, recomputed from
its recorded weight: the weight of one strut times the largest load factor of the active sizing cases is spread over the strut
length, its component perpendicular to the strut loads a pinned beam over the buckling length of the strut.
"""
function strut_own_weight_loads(ac)
    wing, strut = ac.wing, ac.wing.strut
    sl = wing.sizing_loads
    n_max = maximum(abs(getfield(sl, name).n_computed)
                    for name in fieldnames(typeof(sl)) if getfield(sl, name) isa TASOPT.structures.WingLoadCase && getfield(sl, name).active)
    L, z = strut.length, strut.z
    q_perp = n_max * (strut.weight / 2) / L * sqrt(L^2 - z^2) / L
    has_jury = ac.options.opt_wing_config == WingConfig.TrussBraced
    l_b = TASOPT.structures.strut_buckling_length(L, has_jury ? wing.jury.eta_strut : nothing)
    return (V = q_perp * l_b / 2, M = q_perp * l_b^2 / 8)
end

"Weight [N] of the strut or jury strut `member` for both wings from the areas of its `box`, with the added weight fraction `frac`."
box_member_weight(member, box; frac = member.weight_frac_added) =
    2 * TASOPT.gee * member.length * (member.caps.material.ρ * box.area_caps + member.webs.material.ρ * box.area_webs) * (1 + frac)

# Sizing a strut-braced wing with the strut force from the compatibility of the wing and strut deflections.
# The aircraft are sized once and shared; `edit!` changes the strut first. The break is moved outboard to ηs = 0.9,
# where the Prescribed strut force sizes a nearly weightless wing.
function sized_braced_aircraft(key = :compatibility, edit! = ac -> nothing)
    return sized_once(key) do
        ac = read_aircraft_with_wing_config("strut_braced")
        ac.wing.layout.ηs = 0.9
        ac.wing.η_engine = 0.9 # the engine stays at the break
        edit!(ac)
        ac
    end
end

# Truss-braced wing like `sized_jury_aircraft("truss_braced", 0.3)` but without added weight fractions of strut and jury
function sized_without_added_weight_fraction()
    return sized_once(:no_added_weight_fraction) do
        ac = read_aircraft_with_wing_config("truss_braced")
        ac.wing.layout.ηs = 0.7
        ac.wing.jury.eta_strut = 0.3
        ac.wing.strut.weight_frac_added = 0.0
        ac.wing.jury.weight_frac_added = 0.0
        ac
    end
end

"Pitching-moment arm of the strut about the wing-box x."
strut_dx_arm(ac) = 0.25 * ac.wing.layout.span * (ac.wing.layout.ηs - ac.wing.layout.ηo) * tand(ac.wing.layout.sweep)

"Pitching-moment arm of the jury about the wing-box x: the jury hangs at eta_wing, on the swept wing."
jury_dx_arm(ac) = 0.5 * ac.wing.layout.span * (ac.wing.jury.eta_wing - ac.wing.layout.ηo) * tand(ac.wing.layout.sweep)

# Strut-braced wing with the prescribed strut force and only the pull-up case active, sized once and shared
function sized_prescribed_pullup_aircraft()
    return sized_once(:prescribed_pullup) do
        ac = read_aircraft_with_wing_config("strut_braced")
        ac.wing.strut.opt_strut_force = TASOPT.StrutForceModel.Prescribed
        for name in (:pushover, :pullup_light, :pushover_light, :taxi, :gust_up_light, :gust_down_light)
            getfield(ac.wing.sizing_loads, name).active = false
        end
        ac
    end
end

# Wing with the strut attached at ηs = 0.7 (the default ηs does not converge with a strut), sized once and shared;
# `eta_strut` sets the jury of a truss-braced wing
function sized_jury_aircraft(config, eta_strut = 0.5)
    return sized_once((:jury, config, eta_strut)) do
        ac = read_aircraft_with_wing_config(config)
        ac.wing.layout.ηs = 0.7
        config == "truss_braced" && (ac.wing.jury.eta_strut = eta_strut)
        ac
    end
end

# ============================================================================
# I/O: TOML input, deprecated keys, defaults, save and read round trips, summary and parameter dump output
# ============================================================================

# Deprecated strut keys that sit directly in [Wing], and the [Wing.Strut] keys that replace them
const deprecated_flat_strut_keys = (
    z_strut = :z_fuselage_attach,
    strut_toc = :thickness_to_chord,
    strut_local_velocity_ratio = :local_velocity_ratio,
)

mentions_deprecated_strut_key(log, key) = log.level == CoreLogging.Info && occursin(string(key), log.message)

@testset "TBW wing configuration option" verbose=true begin
    ac_default = CoreLogging.with_logger(CoreLogging.NullLogger()) do
        load_default_model()
    end

    @testset "wing_config key selects the wing configuration" begin
        @test !ac_default.options.has_strut
        @test ac_default.options.opt_wing_config == WingConfig.Cantilever
        @test_throws Exception (ac_default.options.has_strut = true)

        for (name, config, strut) in (("cantilever",   WingConfig.Cantilever,  false),
                                      ("strut_braced", WingConfig.StrutBraced, true),
                                      ("truss_braced", WingConfig.TrussBraced, true))
            ac = read_aircraft_with_wing_config(name)
            @test ac.options.opt_wing_config == config
            @test ac.options.has_strut == strut
        end

        @test_throws ErrorException read_aircraft_with_wing_config("biplane")
    end

    @testset "has_strut key is a deprecated alias" begin
        has_strut_txt = insert_below_header(default_toml_without_wing_config(), "Wing", "has_strut = true")

        ac_has_strut, logs = read_aircraft_from_text(has_strut_txt)
        @test ac_has_strut.options.opt_wing_config == WingConfig.StrutBraced
        @test any(r -> r.level == CoreLogging.Info && occursin("has_strut", r.message) &&
                       occursin("wing_config", r.message), logs)

        # wing_config takes precedence when both are given
        both_txt = insert_below_header(has_strut_txt, "Options", "wing_config = \"cantilever\"")
        ac_both, _ = read_aircraft_from_text(both_txt)
        @test ac_both.options.opt_wing_config == WingConfig.Cantilever

        # wing_config alone gives no deprecation message
        _, logs_new = read_aircraft_from_text(toml_with_wing_config("strut_braced"))
        @test !any(r -> occursin("has_strut", r.message), logs_new)
    end

    @testset "save and read keep the wing configuration" begin
        for (name, config) in (("strut_braced", WingConfig.StrutBraced), ("truss_braced", WingConfig.TrussBraced))
            mktempdir() do dir
                path = joinpath(dir, "wing_config_roundtrip.toml")
                save_aircraft_model(read_aircraft_with_wing_config(name), path)
                saved = TASOPT.TOML.parsefile(path)
                @test saved["Options"]["wing_config"] == name
                @test !haskey(saved["Wing"], "has_strut")
                @test first(read_aircraft_from_text(read(path, String))).options.opt_wing_config == config
            end
        end
    end

    @testset "Options summary and show report the wing configuration" begin
        options = ac_default.options
        summary_text = mktemp() do path, io
            redirect_stdout(() -> summary(options), io)
            close(io)
            read(path, String)
        end
        @test occursin("cantilever", summary_text)
        @test occursin("cantilever", sprint(show, options))

        options_strut = deepcopy(options)
        options_strut.opt_wing_config = WingConfig.StrutBraced
        @test occursin("strut_braced", sprint(show, options_strut))
    end

    @testset "default_tbw.toml reads as a truss-braced high wing with fuel in the wing" begin
        ac, logs = read_aircraft_from_text(read(joinpath(TASOPT.__TASOPTroot__,
                                                         "../example/defaults/default_tbw.toml"), String))
        wing = ac.wing
        @test ac.options.opt_wing_config == WingConfig.TrussBraced
        @test !any(r -> occursin("has_strut", r.message) ||
                        any(key -> mentions_deprecated_strut_key(r, key), keys(deprecated_flat_strut_keys)), logs)
        @test wing.layout.AR ≈ 16
        @test wing.layout.sweep ≈ 20
        @test 0.09 ≤ wing.inboard.cross_section.thickness_to_chord ≤ 0.11
        @test 0.55 ≤ wing.layout.ηs ≤ 0.65
        @test wing.layout.z > 0
        @test ac.options.opt_fuel == TASOPT.FuelType.JetA
        @test ac.options.has_wing_fuel
    end
end

# Inputs of the strut and of the jury strut of a truss-braced wing. The strut and the jury strut are boxes with caps
# and webs: material, chord ratio and minimum gauges are inputs.

# Strut input keys in [Wing.Strut] and the `ac.wing.strut` fields they populate
const strut_block_fields = (
    z_fuselage_attach = :z,
    y_fuselage_attach = :y_root,
    thickness_to_chord = :thickness_to_chord,
    local_velocity_ratio = :local_velocity_ratio,
    box_width_to_chord = :box_width_to_chord,
    web_to_box_height = :web_to_box_height,
    end_fixity_K = :end_fixity_K,
    weight_frac_added = :weight_frac_added,
    weight_frac_attachments = :weight_frac_attachments,
    opt_strut_force = :opt_strut_force,
)

# [Wing.Strut] body with values that differ from the defaults; distances carry units
const custom_strut_block = """
    z_fuselage_attach = "170 in"
    y_fuselage_attach = "1.5 ft"
    thickness_to_chord = 0.12
    local_velocity_ratio = 1.1
    box_width_to_chord = 0.45
    web_to_box_height = 0.7
    end_fixity_K = 0.8
    weight_frac_added = 0.5
    weight_frac_attachments = 0.08
    opt_strut_force = "prescribed"
"""

"Default input text without the strut inputs, whether given as flat [Wing] keys or as a [Wing.Strut] table."
function default_toml_without_strut_inputs(config)
    return replace(toml_with_wing_config(config),
        r"^[ \t]*(z_strut|strut_toc|strut_local_velocity_ratio)[ \t]*=.*\n"m => "",
        r"^\[Wing\.Strut\]\n(?:(?!\[).*\n)*"m => "")
end

"Default input text for `config` with a [Wing.Strut] table of `block` and flat [Wing] keys of `flat`."
function toml_with_strut_inputs(config; block = "", flat = "")
    txt = default_toml_without_strut_inputs(config)
    isempty(flat) || (txt = insert_below_header(txt, "Wing", flat))
    isempty(block) || (txt = replace(txt, "[Wing.Aero]" => "[Wing.Strut]\n$block\n[Wing.Aero]"; count = 1))
    return txt
end

"Default input text for `config` with a [Wing.Jury] table of `block` (no table if `block` is empty)."
function toml_with_jury_inputs(config, block = "")
    txt = replace(toml_with_wing_config(config), r"^\[Wing\.Jury\]\n(?:(?!\[).*\n)*"m => "")
    isempty(block) && return txt
    return replace(txt, "[Wing.Aero]" => "[Wing.Jury]\n$block\n[Wing.Aero]"; count = 1)
end

"Aircraft read from the default input for `config`, after `edit!` is applied to the parsed Dict."
function read_aircraft_with_box_inputs(edit!, config)
    return first(read_aircraft_from_default_dict() do d
        d["Options"]["wing_config"] = config
        edit!(d)
    end)
end

"The [Wing.`name`] table of the parsed input `d`, created if missing."
wing_subtable(d, name) = get!(d["Wing"], name, Dict{String,Any}())

"Name, properties and factors of an alloy, to compare materials by value."
alloy_signature(m) = (m.name, m.σmax, m.E, m.ρ, m.τmax, m.max_avg_stress, m.safety_factor)

"Inputs of the strut box: caps and webs materials, chord ratio and minimum gauges."
strut_box_inputs(strut) = (alloy_signature(strut.caps.material), alloy_signature(strut.webs.material),
                           strut.chord_to_wing_chord, strut.min_cap_thickness, strut.min_web_thickness)

"Inputs of the jury strut box: caps and webs materials, chord ratio and minimum gauges."
jury_box_inputs(jury) = (alloy_signature(jury.caps.material), alloy_signature(jury.webs.material),
                         jury.chord_to_strut_chord, jury.min_cap_thickness, jury.min_web_thickness)

@testset "strut and jury strut TOML inputs" verbose=true begin
    custom_values = (
        z_fuselage_attach = TASOPT.Distance("170 in"),
        y_fuselage_attach = TASOPT.Distance("1.5 ft"),
        thickness_to_chord = 0.12,
        local_velocity_ratio = 1.1,
        box_width_to_chord = 0.45,
        web_to_box_height = 0.7,
        end_fixity_K = 0.8,
        weight_frac_added = 0.5,
        weight_frac_attachments = 0.08,
        opt_strut_force = TASOPT.StrutForceModel.Prescribed,
    )

    @testset "[Wing.Strut] populates the strut fields" begin
        for config in ("strut_braced", "truss_braced")
            ac, logs = read_aircraft_from_text(toml_with_strut_inputs(config; block = custom_strut_block))
            for (key, field) in pairs(strut_block_fields)
                @test getfield(ac.wing.strut, field) == custom_values[key]
            end
            @test ac.wing.strut.z ≈ 170 * 0.0254
            @test ac.wing.strut.y_root ≈ 1.5 * 0.3048
            @test !any(r -> any(key -> mentions_deprecated_strut_key(r, key), keys(deprecated_flat_strut_keys)), logs)
        end
    end

    @testset "flat [Wing] strut keys are deprecated but honoured" begin
        flat = """
            z_strut = "160 in"
            strut_toc = 0.13
            strut_local_velocity_ratio = 1.2
        """
        ac, logs = read_aircraft_from_text(toml_with_strut_inputs("strut_braced"; flat))
        @test ac.wing.strut.z ≈ 160 * 0.0254
        @test ac.wing.strut.thickness_to_chord == 0.13
        @test ac.wing.strut.local_velocity_ratio == 1.2
        for key in keys(deprecated_flat_strut_keys)
            @test any(r -> mentions_deprecated_strut_key(r, key), logs)
        end

        # [Wing.Strut] wins over the flat keys; flat keys not in the block are still honoured
        partial_block = """
            z_fuselage_attach = "170 in"
            thickness_to_chord = 0.12
        """
        ac_both, _ = read_aircraft_from_text(toml_with_strut_inputs("strut_braced"; flat, block = partial_block))
        @test ac_both.wing.strut.z ≈ 170 * 0.0254
        @test ac_both.wing.strut.thickness_to_chord == 0.12
        @test ac_both.wing.strut.local_velocity_ratio == 1.2
    end

    @testset "default input gives the default strut values" begin
        ac, logs = read_aircraft_from_text(toml_with_wing_config("strut_braced"))
        strut = ac.wing.strut
        @test strut.z ≈ 154 * 0.0254
        @test strut.thickness_to_chord == 0.15
        @test strut.local_velocity_ratio == 1.0
        @test strut.box_width_to_chord == 0.5
        @test strut.web_to_box_height == 0.75
        @test strut.end_fixity_K == 1.0
        @test strut.weight_frac_added == 0.60
        @test TASOPT.structures.Strut().weight_frac_added == 0.60
        @test strut.weight_frac_attachments == 0.06
        @test TASOPT.structures.Strut().weight_frac_attachments == 0.06
        @test TASOPT.structures.JuryStrut().weight_frac_added == 0.60
        @test !any(r -> any(key -> mentions_deprecated_strut_key(r, key), keys(deprecated_flat_strut_keys)), logs)

        # the jury of a truss-braced wing has its own added weight fraction, 0.60 as well
        ac_truss = first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced")))
        @test ac_truss.wing.jury.weight_frac_added == 0.60
        @test ac_truss.wing.strut.weight_frac_added == 0.60

        # A cantilever wing does not read the strut inputs
        _, logs_cantilever = read_aircraft_from_text(toml_with_wing_config("cantilever"))
        @test !any(r -> any(key -> mentions_deprecated_strut_key(r, key), keys(deprecated_flat_strut_keys)), logs_cantilever)
    end

    @testset "save and read reproduce the strut inputs exactly" begin
        for config in ("strut_braced", "truss_braced")
            ac = first(read_aircraft_from_text(toml_with_strut_inputs(config; block = custom_strut_block)))
            mktempdir() do dir
                path = joinpath(dir, "strut_roundtrip.toml")
                save_aircraft_model(ac, path)
                saved = TASOPT.TOML.parsefile(path)
                @test saved["Wing"]["Strut"] isa AbstractDict
                @test saved["Wing"]["Strut"]["opt_strut_force"] == "prescribed"
                for key in keys(strut_block_fields)
                    @test haskey(saved["Wing"]["Strut"], string(key))
                end
                for key in keys(deprecated_flat_strut_keys)
                    @test !haskey(saved["Wing"], string(key))
                end

                ac_reread, logs = read_aircraft_from_text(read(path, String))
                for field in values(strut_block_fields)
                    @test getfield(ac_reread.wing.strut, field) == getfield(ac.wing.strut, field)
                end
                @test !any(r -> any(key -> mentions_deprecated_strut_key(r, key), keys(deprecated_flat_strut_keys)), logs)
            end
        end
    end

    @testset "opt_strut_force selects the strut force model" begin
        # compatibility is the default, with or without a [Wing.Strut] table
        for txt in (toml_with_wing_config("strut_braced"),
                    toml_with_strut_inputs("strut_braced"),
                    toml_with_strut_inputs("strut_braced"; block = "thickness_to_chord = 0.12"))
            @test first(read_aircraft_from_text(txt)).wing.strut.opt_strut_force ==
                  TASOPT.StrutForceModel.Compatibility
        end

        for (name, model) in (("prescribed",    TASOPT.StrutForceModel.Prescribed),
                              ("compatibility", TASOPT.StrutForceModel.Compatibility))
            for config in ("strut_braced", "truss_braced")
                block = "opt_strut_force = \"$name\""
                ac = first(read_aircraft_from_text(toml_with_strut_inputs(config; block)))
                @test ac.wing.strut.opt_strut_force == model
            end
        end

        @test_throws Exception read_aircraft_from_text(
            toml_with_strut_inputs("strut_braced"; block = "opt_strut_force = \"rigid\""))
    end

    al_default = StructuralAlloy("TASOPT-Al", max_avg_stress = 1.1, safety_factor = 1.5)
    al_7075 = StructuralAlloy("Al-7075", max_avg_stress = 1.2, safety_factor = 2.0)
    default_strut_box = (alloy_signature(al_default), alloy_signature(al_default), 0.5, 1.5e-3, 1.0e-3)
    default_jury_box = (alloy_signature(al_default), alloy_signature(al_default), 0.5, 1.5e-3, 1.0e-3)

    # distinct values from the defaults; the distances carry units or are in metres
    function set_box_inputs!(d)
        d["Structures"]["strut_material"] = "Al-7075"
        d["Structures"]["strut_max_avg_stress"] = 1.2
        d["Structures"]["strut_safety_factor"] = 2.0
        d["Structures"]["strut_web_material"] = "TASOPT-Al"
        d["Structures"]["strut_web_max_avg_stress"] = 1.3
        d["Structures"]["strut_web_safety_factor"] = 1.7
        strut_table = wing_subtable(d, "Strut")
        strut_table["chord_to_wing_chord"] = 0.8
        strut_table["min_cap_thickness"] = "0.1 in"
        strut_table["min_web_thickness"] = 0.0007
        jury_table = wing_subtable(d, "Jury")
        jury_table["chord_to_strut_chord"] = 0.4
        jury_table["min_cap_thickness"] = 0.003
        jury_table["min_web_thickness"] = "0.05 in"
    end
    al_web = StructuralAlloy("TASOPT-Al", max_avg_stress = 1.3, safety_factor = 1.7)
    explicit_strut_box = (alloy_signature(al_7075), alloy_signature(al_web), 0.8, 0.1 * 0.0254, 0.0007)
    explicit_jury_box = (alloy_signature(al_7075), alloy_signature(al_web), 0.4, 0.003, 0.05 * 0.0254)

    @testset "explicit keys set the strut and the jury; the webs differ from the caps and the jury follows the strut" begin
        ac_strut = read_aircraft_with_box_inputs(set_box_inputs!, "strut_braced")
        @test strut_box_inputs(ac_strut.wing.strut) == explicit_strut_box
        @test alloy_signature(ac_strut.wing.strut.caps.material) != alloy_signature(ac_strut.wing.strut.webs.material)

        # a strut-braced wing does not read [Wing.Jury]
        @test jury_box_inputs(ac_strut.wing.jury) == default_jury_box

        ac_truss = read_aircraft_with_box_inputs(set_box_inputs!, "truss_braced")
        @test strut_box_inputs(ac_truss.wing.strut) == explicit_strut_box
        @test jury_box_inputs(ac_truss.wing.jury) == explicit_jury_box
        @test alloy_signature(ac_truss.wing.jury.caps.material) == alloy_signature(ac_truss.wing.strut.caps.material)
        @test alloy_signature(ac_truss.wing.jury.webs.material) == alloy_signature(ac_truss.wing.strut.webs.material)

        # the wing caps and webs are not touched by the strut inputs
        @test alloy_signature(ac_truss.wing.inboard.caps.material) == alloy_signature(al_default)
        @test alloy_signature(ac_truss.wing.inboard.webs.material) == alloy_signature(al_default)

        # a single gauge key leaves the other at its default
        ac_one = read_aircraft_with_box_inputs("truss_braced") do d
            wing_subtable(d, "Strut")["min_cap_thickness"] = 0.002
            wing_subtable(d, "Jury")["min_web_thickness"] = 0.0004
        end
        @test (ac_one.wing.strut.min_cap_thickness, ac_one.wing.strut.min_web_thickness) == (0.002, 1.0e-3)
        @test (ac_one.wing.jury.min_cap_thickness, ac_one.wing.jury.min_web_thickness) == (1.5e-3, 0.0004)
    end

    @testset "the web keys fall back to the strut keys one by one" begin
        for config in ("strut_braced", "truss_braced")
            # no web keys: webs equal the caps
            ac = read_aircraft_with_box_inputs(config) do d
                d["Structures"]["strut_material"] = "Al-7075"
                d["Structures"]["strut_max_avg_stress"] = 1.2
                d["Structures"]["strut_safety_factor"] = 2.0
                for key in ("strut_web_material", "strut_web_max_avg_stress", "strut_web_safety_factor")
                    delete!(d["Structures"], key)
                end
            end
            @test alloy_signature(ac.wing.strut.caps.material) == alloy_signature(al_7075)
            @test alloy_signature(ac.wing.strut.webs.material) == alloy_signature(al_7075)
            config == "truss_braced" &&
                @test alloy_signature(ac.wing.jury.webs.material) == alloy_signature(al_7075)

            # only the web alloy given: its stress factors come from the strut keys
            ac_alloy = read_aircraft_with_box_inputs(config) do d
                d["Structures"]["strut_max_avg_stress"] = 1.2
                d["Structures"]["strut_safety_factor"] = 2.0
                d["Structures"]["strut_web_material"] = "Al-7075"
                delete!(d["Structures"], "strut_web_max_avg_stress")
                delete!(d["Structures"], "strut_web_safety_factor")
            end
            @test alloy_signature(ac_alloy.wing.strut.webs.material) == alloy_signature(al_7075)
            @test alloy_signature(ac_alloy.wing.strut.caps.material) ==
                  alloy_signature(StructuralAlloy("TASOPT-Al", max_avg_stress = 1.2, safety_factor = 2.0))
        end
    end

    @testset "save and read reproduce the strut and jury inputs exactly" begin
        for config in ("strut_braced", "truss_braced")
            ac = read_aircraft_with_box_inputs(set_box_inputs!, config)
            mktempdir() do dir
                path = joinpath(dir, "strut_box_roundtrip.toml")
                save_aircraft_model(ac, path)
                saved = TASOPT.TOML.parsefile(path)
                @test saved["Structures"]["strut_web_material"] == "TASOPT-Al"
                @test saved["Structures"]["strut_web_max_avg_stress"] == 1.3
                @test saved["Structures"]["strut_web_safety_factor"] == 1.7
                @test saved["Structures"]["strut_material"] == "Al-7075"
                @test saved["Wing"]["Strut"]["chord_to_wing_chord"] == 0.8
                @test saved["Wing"]["Strut"]["min_cap_thickness"] == ac.wing.strut.min_cap_thickness
                @test saved["Wing"]["Strut"]["min_web_thickness"] == 0.0007
                if config == "truss_braced"
                    @test saved["Wing"]["Jury"]["chord_to_strut_chord"] == 0.4
                    @test saved["Wing"]["Jury"]["min_cap_thickness"] == 0.003
                    @test saved["Wing"]["Jury"]["min_web_thickness"] == ac.wing.jury.min_web_thickness
                end

                ac_reread, _ = read_aircraft_from_text(read(path, String))
                @test strut_box_inputs(ac_reread.wing.strut) == strut_box_inputs(ac.wing.strut)
                @test jury_box_inputs(ac_reread.wing.jury) == jury_box_inputs(ac.wing.jury)
                config == "truss_braced" && @test jury_box_inputs(ac_reread.wing.jury) == explicit_jury_box
            end
        end
    end

    @testset "[Wing.Jury] sets eta_strut of a truss-braced wing, and save and read reproduce it exactly" begin
        # the jury attachment is 0.5 without a table
        @test first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced"))).wing.jury.eta_strut == 0.5
        @test first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced", "eta_strut = 0.3"))).wing.jury.eta_strut == 0.3

        # a strut-braced wing does not read it
        ac = first(read_aircraft_from_text(toml_with_jury_inputs("strut_braced", "eta_strut = 0.3")))
        @test ac.wing.jury.eta_strut == 0.5

        ac = first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced", "eta_strut = 0.3")))
        mktempdir() do dir
            path = joinpath(dir, "jury_roundtrip.toml")
            save_aircraft_model(ac, path)
            saved = TASOPT.TOML.parsefile(path)
            @test saved["Wing"]["Jury"]["eta_strut"] == 0.3
            @test first(read_aircraft_from_text(read(path, String))).wing.jury.eta_strut == 0.3
        end

        # only a truss-braced wing has a jury table
        mktempdir() do dir
            path = joinpath(dir, "no_jury.toml")
            save_aircraft_model(read_aircraft_with_wing_config("strut_braced"), path)
            @test !haskey(TASOPT.TOML.parsefile(path)["Wing"], "Jury")
        end
    end
end

# Wing sizing load cases in [Wing.SizingLoads.<case>] tables

const sizing_load_case_names = (:pullup, :pushover, :pullup_light, :pushover_light,
                                :taxi, :gust_up_light, :gust_down_light)

"The [Wing.SizingLoads.`case`] table of the parsed input `d`, created if missing."
sizing_load_table(d, case) =
    get!(get!(d["Wing"], "SizingLoads", Dict{String,Any}()), string(case), Dict{String,Any}())

"The inputs of load case `name` that the TOML block carries."
function case_inputs(loads, name)
    case = getfield(loads, name)
    return (case.load_factor, case.mass_state, case.active, case.tail_load_scale)
end

mentions_nlift_deprecation(log) =
    log.level == CoreLogging.Info && occursin("Nlift", log.message) && occursin("deprecat", log.message)

@testset "wing sizing load case inputs" verbose=true begin
    @testset "defaults" begin
        # local bindings, so the names cannot clash with anything TASOPT exports
        WingMassState = TASOPT.structures.WingMassState
        WingLoadKind = TASOPT.structures.WingLoadKind
        WingLoadCase = TASOPT.structures.WingLoadCase
        WingCaseLoads = TASOPT.structures.WingCaseLoads
        WingSizingLoads = TASOPT.structures.WingSizingLoads

        loads = WingSizingLoads()

        @test Tuple(n for n in fieldnames(WingSizingLoads)
                    if fieldtype(WingSizingLoads, n) === WingLoadCase) == sizing_load_case_names

        # only the pull-up case is active
        for name in sizing_load_case_names
            @test getfield(loads, name).active == (name == :pullup)
        end

        @test loads.pullup.kind == WingLoadKind.Maneuver
        @test loads.pullup.load_factor == 3.0
        @test loads.pullup.mass_state == WingMassState.DesignCruise

        expected = (pushover        = (WingLoadKind.Maneuver, -1.0, WingMassState.MTOW_MaxFuel),
                    pullup_light    = (WingLoadKind.Maneuver,  3.0, WingMassState.Light),
                    pushover_light  = (WingLoadKind.Maneuver, -1.0, WingMassState.Light),
                    taxi            = (WingLoadKind.Ground,   -2.0, WingMassState.MTOW_MaxFuel),
                    gust_up_light   = (WingLoadKind.Gust,      1.0, WingMassState.Light),
                    gust_down_light = (WingLoadKind.Gust,     -1.0, WingMassState.Light))
        for (name, (kind, load_factor, mass_state)) in pairs(expected)
            case = getfield(loads, name)
            @test case.kind == kind
            @test case.load_factor == load_factor
            @test case.mass_state == mass_state
        end

        @test isempty(loads.results)
        @test loads.governing_root == :pullup
        @test loads.governing_break == :pullup
        @test loads.governing_strut == :none

        # each Wing carries its own sizing loads
        wing_a, wing_b = TASOPT.structures.Wing(), TASOPT.structures.Wing()
        wing_a.sizing_loads.pullup.load_factor = 2.5
        @test wing_b.sizing_loads.pullup.load_factor == 3.0
    end

    @testset "TOML input" verbose=true begin
        WingMassState = TASOPT.structures.WingMassState
        WingSizingLoads = TASOPT.structures.WingSizingLoads

        @testset "[Wing.SizingLoads] populates the load cases" begin
            ac, _ = read_aircraft_from_default_dict() do d
                pushover = sizing_load_table(d, :pushover)
                pushover["active"] = true
                pushover["load_factor"] = -1.5
                pushover["mass_state"] = "light"
                taxi = sizing_load_table(d, :taxi)
                taxi["active"] = true
                taxi["tail_load_scale"] = 0.8
            end
            loads = ac.wing.sizing_loads

            @test loads.pushover.active
            @test loads.pushover.load_factor == -1.5
            @test loads.pushover.mass_state == WingMassState.Light
            @test loads.taxi.active
            @test loads.taxi.tail_load_scale == 0.8
            @test loads.taxi.load_factor == -2.0
            @test loads.taxi.mass_state == WingMassState.MTOW_MaxFuel

            # cases that the input does not touch keep their defaults
            for name in (:pullup, :pullup_light, :pushover_light, :gust_up_light, :gust_down_light)
                @test case_inputs(loads, name) == case_inputs(WingSizingLoads(), name)
            end

            # missing keys fall back to the defaults
            ac_missing, _ = read_aircraft_from_default_dict() do d
                delete!(sizing_load_table(d, :pushover_light), "load_factor")
                delete!(sizing_load_table(d, :pushover_light), "mass_state")
            end
            @test case_inputs(ac_missing.wing.sizing_loads, :pushover_light) ==
                  case_inputs(WingSizingLoads(), :pushover_light)
        end

        @testset "save and read reproduce the load cases exactly" begin
            ac, _ = read_aircraft_from_default_dict() do d
                sizing_load_table(d, :pullup)["load_factor"] = 2.6
                sizing_load_table(d, :pushover)["active"] = true
                sizing_load_table(d, :pushover)["mass_state"] = "mzfw"
                sizing_load_table(d, :pullup_light)["active"] = true
                sizing_load_table(d, :pullup_light)["load_factor"] = 2.25
                sizing_load_table(d, :taxi)["active"] = true
                sizing_load_table(d, :taxi)["tail_load_scale"] = 0.75
                sizing_load_table(d, :gust_up_light)["mass_state"] = "DesignCruise"
                sizing_load_table(d, :gust_down_light)["load_factor"] = -0.5
            end

            mktempdir() do dir
                path = joinpath(dir, "sizing_loads_roundtrip.toml")
                save_aircraft_model(ac, path)
                saved = TASOPT.TOML.parsefile(path)

                @test !haskey(saved["Mission"], "Nlift")
                @test saved["Wing"]["SizingLoads"] isa AbstractDict
                for name in sizing_load_case_names
                    table = saved["Wing"]["SizingLoads"][string(name)]
                    for key in ("load_factor", "mass_state", "active", "tail_load_scale")
                        @test haskey(table, key)
                    end
                    for key in ("kind", "n_computed", "gust_velocity_EAS", "gust_speed_EAS")
                        @test !haskey(table, key)
                    end
                end

                ac_reread, logs = read_aircraft_from_text(read(path, String))
                for name in sizing_load_case_names
                    @test case_inputs(ac_reread.wing.sizing_loads, name) == case_inputs(ac.wing.sizing_loads, name)
                end
                @test !any(mentions_nlift_deprecation, logs)
            end
        end

        @testset "Mission.Nlift is a deprecated alias of the pull-up load factor" begin
            ac, logs = read_aircraft_from_default_dict() do d
                delete!(sizing_load_table(d, :pullup), "load_factor")
                d["Mission"]["Nlift"] = 2.7
            end
            @test ac.wing.sizing_loads.pullup.load_factor == 2.7
            @test ac.parg[igNlift] == 0.0
            @test any(mentions_nlift_deprecation, logs)

            # the table wins when both are given
            ac_both, _ = read_aircraft_from_default_dict() do d
                sizing_load_table(d, :pullup)["load_factor"] = 2.4
                d["Mission"]["Nlift"] = 2.7
            end
            @test ac_both.wing.sizing_loads.pullup.load_factor == 2.4
        end

        # The user input with the `active` keys removed, so that only the defaults decide
        function without_active_keys!(d)
            for name in sizing_load_case_names
                delete!(sizing_load_table(d, name), "active")
            end
        end
        braced_extra_cases = (:pushover, :pullup_light, :pushover_light, :taxi)
        active_cases(ac) = Set(n for n in sizing_load_case_names if getfield(ac.wing.sizing_loads, n).active)

        @testset "braced wings activate the extra load cases by default" begin
            for config in ("strut_braced", "truss_braced")
                ac, _ = read_aircraft_from_default_dict() do d
                    d["Options"]["wing_config"] = config
                    without_active_keys!(d)
                end
                @test active_cases(ac) == Set([:pullup, braced_extra_cases...])
            end

            # default_tbw.toml has no [Wing.SizingLoads] table
            ac_tbw, _ = read_aircraft_from_text(read(joinpath(TASOPT.__TASOPTroot__,
                                                              "../example/defaults/default_tbw.toml"), String))
            @test active_cases(ac_tbw) == Set([:pullup, braced_extra_cases...])

            # a cantilever wing keeps only the pull-up case
            ac_cantilever, _ = read_aircraft_from_default_dict(without_active_keys!)
            @test active_cases(ac_cantilever) == Set([:pullup])
        end

        @testset "an explicit active key wins over the braced defaults" begin
            ac, _ = read_aircraft_from_default_dict() do d
                d["Options"]["wing_config"] = "strut_braced"
                without_active_keys!(d)
                sizing_load_table(d, :taxi)["active"] = false
                sizing_load_table(d, :gust_up_light)["active"] = true
            end
            @test active_cases(ac) == Set([:pullup, :pushover, :pullup_light, :pushover_light, :gust_up_light])
        end

        @testset "unknown mass_state is an error" begin
            @test_throws Exception read_aircraft_from_default_dict() do d
                sizing_load_table(d, :pushover)["mass_state"] = "HalfFull"
            end
        end
    end
end

# FAR/CS 25.337 positive manoeuvre load factor: n = clamp(2.1 + 24000 / (W_lb + 10000), 2.5, 3.8), W_lb = W / 4.4482216152605
const newtons_per_lbf = 4.4482216152605

@testset "FAR positive manoeuvre load factor" verbose=true begin
    far_manoeuvre_load_factor = TASOPT.structures.far_manoeuvre_load_factor

    @testset "hand calculation" begin
        # 20000 lb: 2.1 + 24000/30000 = 2.9
        @test far_manoeuvre_load_factor(20000 * newtons_per_lbf) ≈ 2.9 rtol = 1e-12
        # 10000 lb: 2.1 + 24000/20000 = 3.3
        @test far_manoeuvre_load_factor(10000 * newtons_per_lbf) ≈ 3.3 rtol = 1e-12
        # 5000 lb: 2.1 + 24000/15000 = 3.7
        @test far_manoeuvre_load_factor(5000 * newtons_per_lbf) ≈ 3.7 rtol = 1e-12
        # 88964.43230521 N = 20000 lb
        @test far_manoeuvre_load_factor(88964.43230521) ≈ 2.9 rtol = 1e-9
        # the load factor reaches 3.8 at W_lb = 24000/1.7 - 10000 = 4117.6
        @test far_manoeuvre_load_factor(4117.6 * newtons_per_lbf) ≈ 3.8 rtol = 1e-5
        # 50000 lb: 2.1 + 24000/60000 = 2.5, the lower bound
        @test far_manoeuvre_load_factor(50000 * newtons_per_lbf) ≈ 2.5 rtol = 1e-12
    end

    @testset "the load factor is limited to 2.5 and 3.8" begin
        # 4000 lb: 2.1 + 24000/14000 = 3.8143 -> 3.8
        @test far_manoeuvre_load_factor(4000 * newtons_per_lbf) == 3.8
        @test far_manoeuvre_load_factor(0.0) == 3.8
        # 5.0e5 N = 112404 lb: 2.1 + 24000/122404 = 2.2961 -> 2.5
        @test far_manoeuvre_load_factor(5.0e5) == 2.5
        # 145000 lb: 2.1 + 24000/155000 = 2.2548 -> 2.5
        @test far_manoeuvre_load_factor(145000 * newtons_per_lbf) == 2.5
        # 200000 lb: 2.1 + 24000/210000 = 2.2143 -> 2.5
        @test far_manoeuvre_load_factor(200000 * newtons_per_lbf) == 2.5
        @test far_manoeuvre_load_factor(1.0e8) == 2.5
    end

    @testset "opt_far_load_factor is off by default" begin
        ac = CoreLogging.with_logger(CoreLogging.NullLogger()) do
            load_default_model()
        end
        @test ac.options.opt_far_load_factor === false
        @test ac.wing.sizing_loads.pullup.load_factor == 3.0
    end

    @testset "far_load_factor key in [Options] sets the option, and save and read reproduce it exactly" begin
        for value in (true, false)
            ac, _ = read_aircraft_from_default_dict(d -> d["Options"]["far_load_factor"] = value)
            @test ac.options.opt_far_load_factor == value

            mktempdir() do dir
                path = joinpath(dir, "far_load_factor_roundtrip.toml")
                save_aircraft_model(ac, path)
                @test TASOPT.TOML.parsefile(path)["Options"]["far_load_factor"] == value
                @test first(read_aircraft_from_text(read(path, String))).options.opt_far_load_factor == value
            end
        end
    end

    @testset "saving keeps the user input load factor of the pull-up cases" begin
        ac, _ = read_aircraft_from_default_dict() do d
            d["Options"]["far_load_factor"] = true
            sizing_load_table(d, :pullup)["load_factor"] = 2.6
            sizing_load_table(d, :pullup_light)["load_factor"] = 2.25
        end
        mktempdir() do dir
            path = joinpath(dir, "far_load_factor_inputs.toml")
            save_aircraft_model(ac, path)
            saved_cases = TASOPT.TOML.parsefile(path)["Wing"]["SizingLoads"]
            @test saved_cases["pullup"]["load_factor"] == 2.6
            @test saved_cases["pullup_light"]["load_factor"] == 2.25
        end
    end

    @testset "Options summary and show report the option" begin
        options = CoreLogging.with_logger(CoreLogging.NullLogger()) do
            load_default_model().options
        end
        summary_text(opt) = mktemp() do path, io
            redirect_stdout(() -> summary(opt), io)
            close(io)
            read(path, String)
        end
        @test occursin("FAR Load Factor: false", summary_text(options))
        @test occursin("FAR Load Factor: false", sprint(show, options))

        options_far = deepcopy(options)
        options_far.opt_far_load_factor = true
        @test occursin("FAR Load Factor: true", summary_text(options_far))
        @test occursin("FAR Load Factor: true", sprint(show, options_far))
    end
end

const braced_min_cap_thickness = 1.5e-3 # [m]
const braced_min_web_thickness = 1.0e-3 # [m]
const wing_config_names = ("cantilever", "strut_braced", "truss_braced")

"Aircraft read from the default input for `config`, with the minimum gauge keys removed from the user input and `edit!` applied to the Dict."
function read_aircraft_with_gauge_inputs(edit!, config)
    return first(read_aircraft_from_default_dict() do d
        d["Options"]["wing_config"] = config
        delete!(d["Wing"], "min_cap_thickness")
        delete!(d["Wing"], "min_web_thickness")
        edit!(d)
    end)
end

read_aircraft_with_gauge_inputs(config) = read_aircraft_with_gauge_inputs(d -> nothing, config)

"`(min_cap_thickness, min_web_thickness)` of a wing section."
gauge_floors(section) = (section.min_cap_thickness, section.min_web_thickness)

@testset "TBW minimum wing gauge" verbose=true begin
    @testset "defaults depend on the wing configuration" begin
        for config in ("strut_braced", "truss_braced")
            ac = read_aircraft_with_gauge_inputs(config)
            @test gauge_floors(ac.wing.inboard)  == (braced_min_cap_thickness, braced_min_web_thickness)
            @test gauge_floors(ac.wing.outboard) == (braced_min_cap_thickness, braced_min_web_thickness)
            # tails have no floor
            for tail in (ac.htail, ac.vtail)
                @test gauge_floors(tail.inboard) == (0.0, 0.0)
                @test gauge_floors(tail.outboard) == (0.0, 0.0)
            end
        end

        ac_cantilever = read_aircraft_with_gauge_inputs("cantilever")
        @test gauge_floors(ac_cantilever.wing.inboard) == (0.0, 0.0)
        @test gauge_floors(ac_cantilever.wing.outboard) == (0.0, 0.0)

        # default_tbw.toml is a braced wing
        ac_tbw, _ = read_aircraft_from_text(read(joinpath(TASOPT.__TASOPTroot__,
                                                          "../example/defaults/default_tbw.toml"), String))
        @test all(>(0), gauge_floors(ac_tbw.wing.inboard))
        @test gauge_floors(ac_tbw.wing.outboard) == gauge_floors(ac_tbw.wing.inboard)
    end

    @testset "explicit [Wing] keys set both sections and win over the defaults" begin
        for config in wing_config_names
            ac = read_aircraft_with_gauge_inputs(config) do d
                d["Wing"]["min_cap_thickness"] = 0.002
                d["Wing"]["min_web_thickness"] = 0.0007
            end
            @test gauge_floors(ac.wing.inboard)  == (0.002, 0.0007)
            @test gauge_floors(ac.wing.outboard) == (0.002, 0.0007)
        end

        # distance strings use the units of the other [Wing] distances
        ac = read_aircraft_with_gauge_inputs("strut_braced") do d
            d["Wing"]["min_cap_thickness"] = "0.1 in"
            d["Wing"]["min_web_thickness"] = "0.05 in"
        end
        @test collect(gauge_floors(ac.wing.inboard))  ≈ [0.1 * 0.0254, 0.05 * 0.0254]
        @test collect(gauge_floors(ac.wing.outboard)) ≈ [0.1 * 0.0254, 0.05 * 0.0254]

        # an explicit zero switches the braced default off; a single key leaves the other at its default
        ac_zero = read_aircraft_with_gauge_inputs("strut_braced") do d
            d["Wing"]["min_cap_thickness"] = 0.0
        end
        @test gauge_floors(ac_zero.wing.inboard)  == (0.0, braced_min_web_thickness)
        @test gauge_floors(ac_zero.wing.outboard) == (0.0, braced_min_web_thickness)

        ac_web = read_aircraft_with_gauge_inputs("truss_braced") do d
            d["Wing"]["min_web_thickness"] = 0.0
        end
        @test gauge_floors(ac_web.wing.inboard)  == (braced_min_cap_thickness, 0.0)
        @test gauge_floors(ac_web.wing.outboard) == (braced_min_cap_thickness, 0.0)
    end

    @testset "save and read reproduce the floors exactly" begin
        cases = (("truss_braced", d -> nothing), # the braced defaults
                 ("strut_braced", d -> (d["Wing"]["min_cap_thickness"] = 0.0021; d["Wing"]["min_web_thickness"] = 0.00073)),
                 ("cantilever",   d -> (d["Wing"]["min_cap_thickness"] = 0.0; d["Wing"]["min_web_thickness"] = 0.0)))
        for (config, edit!) in cases
            ac = read_aircraft_with_gauge_inputs(edit!, config)
            mktempdir() do dir
                path = joinpath(dir, "min_gauge_roundtrip.toml")
                save_aircraft_model(ac, path)
                saved = TASOPT.TOML.parsefile(path)
                @test saved["Wing"]["min_cap_thickness"] == ac.wing.inboard.min_cap_thickness
                @test saved["Wing"]["min_web_thickness"] == ac.wing.inboard.min_web_thickness

                ac_reread, _ = read_aircraft_from_text(read(path, String))
                @test gauge_floors(ac_reread.wing.inboard)  == gauge_floors(ac.wing.inboard)
                @test gauge_floors(ac_reread.wing.outboard) == gauge_floors(ac.wing.outboard)
            end
        end
    end
end

# Stress and safety factors of the structural materials are kept in StructuralAlloy

@testset "structural material factors round trip through save and read" verbose=true begin
    material_getters = (
        caps = ac -> ac.wing.inboard.caps.material,
        webs = ac -> ac.wing.inboard.webs.material,
        strut = ac -> ac.wing.strut.caps.material,
        skin = ac -> ac.fuselage.skin.material,
        cone = ac -> ac.fuselage.cone.material,
        bending = ac -> ac.fuselage.bendingmaterial_h.material,
        floor = ac -> ac.fuselage.floor.material,
    )
    # distinct factors per material, which no combined value could reproduce
    factors = Dict(name => (1.2 + 0.01i, 1.7 + 0.01i) for (i, name) in enumerate(keys(material_getters)))

    ac, _ = read_aircraft_from_default_dict() do d
        for (name, (max_avg, safety)) in factors
            d["Structures"]["$(name)_max_avg_stress"] = max_avg
            d["Structures"]["$(name)_safety_factor"] = safety
        end
    end

    @testset "the input keys set the factors of each material, and save and read reproduce them exactly" begin
        for (name, getter) in pairs(material_getters)
            material = getter(ac)
            max_avg, safety = factors[name]
            @test material.max_avg_stress == max_avg
            @test material.safety_factor == safety
            @test material.σmax ≈ material.YTS / max_avg / safety
        end

        mktempdir() do dir
            path = joinpath(dir, "material_factors_roundtrip.toml")
            save_aircraft_model(ac, path)
            saved = TASOPT.TOML.parsefile(path)["Structures"]
            ac_reread, _ = read_aircraft_from_text(read(path, String))
            for (name, getter) in pairs(material_getters)
                max_avg, safety = factors[name]
                @test saved["$(name)_max_avg_stress"] == max_avg
                @test saved["$(name)_safety_factor"] == safety
                @test getter(ac_reread).max_avg_stress == max_avg
                @test getter(ac_reread).safety_factor == safety
                @test getter(ac_reread).σmax ≈ getter(ac).σmax rtol = 1e-14
            end
        end
    end
end

# Inputs of the fuselage strut-attach reinforcement. The Fuselage weight_frac_strut_attach key is currently unused by the
# weight, but is still read and saved.

"Default input text for `config` without the `weight_frac_strut_attach` key, and with `line` below [Fuselage.Weights] if given."
function toml_with_strut_attach_fraction(config, line = "")
    txt = replace(toml_with_wing_config(config), r"^[ \t]*weight_frac_strut_attach[ \t]*=.*\n"m => "")
    return isempty(line) ? txt : insert_below_header(txt, "Fuselage.Weights", line)
end

@testset "fuselage strut-attach reinforcement inputs" verbose=true begin
    attach_fraction(txt) = first(read_aircraft_from_text(txt)).fuselage.weight_frac_strut_attach

    @testset "TOML weight_frac_strut_attach: the default depends on the wing configuration, an explicit key wins" begin
        @test attach_fraction(toml_with_strut_attach_fraction("cantilever")) == 0.0
        @test attach_fraction(toml_with_strut_attach_fraction("strut_braced")) == 0.05
        @test attach_fraction(toml_with_strut_attach_fraction("truss_braced")) == 0.05

        @test attach_fraction(toml_with_strut_attach_fraction("strut_braced", "weight_frac_strut_attach = 0.2")) == 0.2
        @test attach_fraction(toml_with_strut_attach_fraction("truss_braced", "weight_frac_strut_attach = 0.0")) == 0.0
        @test attach_fraction(toml_with_strut_attach_fraction("cantilever", "weight_frac_strut_attach = 0.1")) == 0.1
    end

    @testset "save and read reproduce the fraction exactly" begin
        for (config, line, expected) in (("strut_braced", "", 0.05),
                                         ("truss_braced", "weight_frac_strut_attach = 0.0123", 0.0123),
                                         ("cantilever", "", 0.0))
            ac = first(read_aircraft_from_text(toml_with_strut_attach_fraction(config, line)))
            @test ac.fuselage.weight_frac_strut_attach == expected
            mktempdir() do dir
                path = joinpath(dir, "strut_attach_roundtrip.toml")
                save_aircraft_model(ac, path)
                saved = TASOPT.TOML.parsefile(path)
                @test saved["Fuselage"]["Weights"]["weight_frac_strut_attach"] == expected
                @test attach_fraction(read(path, String)) == expected
            end
        end
    end

    # The reinforcement is this fraction of the strut weight, like the wing attachment hardware
    @testset "[Wing.Strut] weight_frac_attachments: the default is 0.06, an explicit key wins, save and read are exact" begin
        attachments(txt) = first(read_aircraft_from_text(txt)).wing.strut.weight_frac_attachments
        @test TASOPT.structures.Strut().weight_frac_attachments == 0.06

        for config in ("strut_braced", "truss_braced")
            # the default input has no key, and neither has an input without a [Wing.Strut] table
            @test attachments(toml_with_wing_config(config)) == 0.06
            @test attachments(toml_with_strut_inputs(config)) == 0.06
            @test attachments(toml_with_strut_inputs(config; block = "weight_frac_attachments = 0.08")) == 0.08
            @test attachments(toml_with_strut_inputs(config; block = "weight_frac_attachments = 0.0")) == 0.0

            for (block, expected) in (("", 0.06), ("weight_frac_attachments = 0.0123", 0.0123))
                ac = first(read_aircraft_from_text(toml_with_strut_inputs(config; block)))
                mktempdir() do dir
                    path = joinpath(dir, "strut_attachments_roundtrip.toml")
                    save_aircraft_model(ac, path)
                    @test TASOPT.TOML.parsefile(path)["Wing"]["Strut"]["weight_frac_attachments"] == expected
                    @test attachments(read(path, String)) == expected
                end
            end
        end
    end
end

# Engine spanwise location input

@testset "engine_location input" verbose=true begin
    example_path(name) = joinpath(TASOPT.__TASOPTroot__, "../example", name)
    read_example(name) = first(read_aircraft_from_text(read(example_path(name), String)))

    @testset "default input puts the engine at the span break" begin
        wing = read_example("defaults/default_input.toml").wing
        @test wing.η_engine == 0.285
        @test wing.η_engine == wing.layout.ηs
    end

    @testset "default_tbw.toml puts the engine inboard of the strut attachment" begin
        wing = read_example("defaults/default_tbw.toml").wing
        @test wing.layout.ηs == 0.6
        @test wing.η_engine == 0.285
    end

    @testset "a missing key puts the engine at the span break" begin
        wing = read_example("cryo_input.toml").wing
        @test !haskey(TASOPT.TOML.parsefile(example_path("cryo_input.toml"))["Wing"], "engine_location")
        @test wing.η_engine == wing.layout.ηs
    end

    @testset "the key sets the engine location, and save and read reproduce it exactly" begin
        ac, _ = read_aircraft_from_default_dict() do d
            d["Wing"]["engine_location"] = 0.4321
        end
        @test ac.wing.η_engine == 0.4321
        @test ac.wing.layout.ηs == 0.285
        mktempdir() do dir
            path = joinpath(dir, "engine_location_roundtrip.toml")
            save_aircraft_model(ac, path)
            @test TASOPT.TOML.parsefile(path)["Wing"]["engine_location"] == 0.4321
            @test first(read_aircraft_from_text(read(path, String))).wing.η_engine == 0.4321
        end
    end
end

@testset "summary table" verbose=true begin
    WingCaseLoads = TASOPT.structures.WingCaseLoads
    ac = CoreLogging.with_logger(CoreLogging.NullLogger()) do
        load_default_model()
    end
    sl = ac.wing.sizing_loads

    function summary_text(ac)
        io = IOBuffer()
        TASOPT.wing_sizing_loads_summary(ac; io = io)
        String(take!(io))
    end
    # first line that names the case, so that pullup does not match pullup_light
    case_row(text, name) = findfirst(line -> occursin(Regex("\\b$name\\b"), line), split(text, '\n'))
    row_text(text, name) = split(text, '\n')[case_row(text, name)]
    # the leading digits of a value appear in the text, whatever the number format
    shows_value(row, x) = occursin(first(replace(string(round(abs(x), sigdigits = 3)), r"\D" => ""), 3),
                                   replace(row, r"[.,]" => ""))
    marked(row) = occursin(r"\*|root|break|governing"i, row)

    @testset "no evaluated cases" begin
        @test isempty(sl.results)
        text = summary_text(ac)
        @test occursin(r"no .*cases"i, text)
        @test all(isnothing(case_row(text, name)) for name in fieldnames(typeof(sl))[1:7])
    end

    @testset "table of the evaluated cases" begin
        sl.pushover.n_computed, sl.taxi.n_computed = -1.0, -2.0
        # inserted out of order: the rows follow the order of the case fields
        sl.results[:taxi] = WingCaseLoads(S_root = 1.11e5, M_root = 2.22e6, S_break = 3.33e4, M_break = 4.44e5)
        sl.results[:pushover] = WingCaseLoads(S_root = 5.55e5, M_root = 6.66e6, S_break = 7.77e4, M_break = 8.88e5)
        sl.results[:pullup] = WingCaseLoads(S_root = 1.23e5, M_root = 4.56e6, S_break = 7.89e4, M_break = 9.87e5)
        sl.governing_root, sl.governing_break = :pushover, :taxi
        text = summary_text(ac)
        rows = Dict(name => row_text(text, name) for name in (:pullup, :pushover, :taxi))

        @test isnothing(case_row(text, :pullup_light))
        @test case_row(text, :pullup) < case_row(text, :pushover) < case_row(text, :taxi)
        for (name, loads) in sl.results
            @test all(shows_value(rows[name], getfield(loads, f)) for f in (:S_root, :M_root, :S_break, :M_break))
        end
        @test occursin("-1", rows[:pushover])
        @test occursin("-2", rows[:taxi])
        @test marked(rows[:pushover])
        @test marked(rows[:taxi])
        @test !marked(rows[:pullup])
    end

    @testset "geometry does not print the summary" begin
        io = IOBuffer()
        TASOPT.geometry(ac; io = io)
        @test isnothing(case_row(String(take!(io)), :pushover))
    end
end

@testset "regression-test parameter dump of the jury strut" verbose=true begin
    truss = sized_jury_aircraft("truss_braced", 0.3)

    @testset "the regression-test parameter dump writes the jury weight and moment" begin
        # `reset_regression_test` writes into the test directory, so restore the reference files afterwards
        reference_files = [joinpath(TASOPT.__TASOPTroot__, "../test", f) for f in ("default_sized.jl", "default_structures.jl")]
        originals = read.(reference_files, String)
        dump = try
            CoreLogging.with_logger(CoreLogging.NullLogger()) do
                TASOPT.reset_regression_test(truss)
            end
            read(reference_files[2], String)
        finally
            write.(reference_files, originals)
        end
        # NaN if the dump has no line for `name`
        function value(name)
            line = match(Regex("^" * name * " = (\\S+)", "m"), dump)
            return isnothing(line) ? NaN : parse(Float64, line.captures[1])
        end
        @test value("wing.strut.weight") == truss.wing.strut.weight
        @test value("wing.jury.weight") == truss.wing.jury.weight
        @test value("wing.jury.dxW") == truss.wing.jury.dxW
    end
end

# ============================================================================
# Golden values: hand-calculated reference values
# ============================================================================

@testset "wing lift curve slope" verbose=true begin
    wing_lift_curve_slope = TASOPT.wing_lift_curve_slope

    @testset "default_tbw.toml wing at cruise Mach" begin
        ac, _ = read_aircraft_from_text(read(joinpath(TASOPT.__TASOPTroot__,
                                                      "../example/defaults/default_tbw.toml"), String))
        AR, sweep, Mach = ac.wing.layout.AR, ac.wing.layout.sweep, ac.para[iaMach, ipcruise1, 1]
        @test (AR, sweep, Mach) == (16.0, 20.0, 0.8)

        # β = sqrt(1 - 0.8^2) = 0.6 and tand(20)^2 = 0.13247433
        # CLα = 2π·16 / (2 + sqrt(4 + (16·0.6)^2 · (1 + 0.13247433/0.6^2)))
        #     = 100.53096491 / (2 + sqrt(130.07342885)) = 100.53096491 / (2 + 11.40497386) = 7.4995271 1/rad
        @test wing_lift_curve_slope(AR, sweep, Mach) ≈ 7.499527110770484 rtol = 1e-6
    end

    @testset "limits" begin
        # a thin infinite wing in incompressible flow has the 2D slope 2π
        @test wing_lift_curve_slope(1e6, 0.0, 0.0) ≈ 2π rtol = 1e-5
        # incompressible, the slope stays below 2π and grows with the aspect ratio
        @test issorted([wing_lift_curve_slope(AR, 20.0, 0.0) for AR in (5.0, 10.0, 20.0, 30.0)])
        @test all(wing_lift_curve_slope(AR, 20.0, 0.0) < 2π for AR in (5.0, 10.0, 20.0, 30.0))
        # sweep lowers the slope
        @test wing_lift_curve_slope(10.0, 30.0, 0.8) < wing_lift_curve_slope(10.0, 0.0, 0.8)
    end
end

@testset "gust load increment" verbose=true begin
    gust_load_increment = TASOPT.gust_load_increment
    U_de, V_EAS, CL_alpha, c_mean = 15.24, 154.3, 4.5, 4.0   # m/s EAS (50 ft/s, 300 kt), 1/rad, m
    ρ0 = TASOPT.ρSL   # 1.2255627 kg/m³

    # Pratt–Walker: μ = 2(W/S) / (ρ0·c_mean·CLα·g),  K_g = 0.88μ/(5.3+μ),  Δn = ρ0·U_de·V·CLα·K_g / (2(W/S))
    @testset "hand calculation" begin
        # W/S = 5500 N/m²: μ = 11000 / (1.2255627·4·4.5·9.81) = 11000/216.40986 = 50.829476
        #   K_g = 0.88·50.829476/(5.3+50.829476) = 0.79690640
        #   Δn  = 1.2255627·15.24·154.3·4.5·0.79690640 / 11000 = 12968.7746·0.79690640/11000 = 0.93953632
        @test gust_load_increment(5500.0, c_mean, CL_alpha, U_de, V_EAS) ≈ 0.9395363236301301 rtol = 1e-6
        # W/S = 3500 N/m²: μ = 7000/216.40986 = 32.346030, K_g = 0.75610912
        #   Δn  = 12968.7746·0.75610912/7000 = 1.40082982
        @test gust_load_increment(3500.0, c_mean, CL_alpha, U_de, V_EAS) ≈ 1.4008298164662953 rtol = 1e-6
    end
end

# ============================================================================
# Cantilever non-interference: braced-wing options do not change the cantilever wing
# ============================================================================

# On a cantilever wing the pull-up governs, so the extra sizing cases must not change the sized wing

cryo_input_path() = joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml")
sized_cryo_aircraft() = sized_once(() -> read_aircraft_model(cryo_input_path()), :cryo; iter = 50)

@testset "sized cantilever wing" verbose=true begin
    extra_cases = (:pushover, :pullup_light, :pushover_light, :taxi)

    # The recorded case loads use the panel weights of the previous sizing iteration, so they
    # equal the loads of the sized wing only to the convergence tolerance.
    function test_recorded_pullup(ac)
        sl = ac.wing.sizing_loads
        @test Set(keys(sl.results)) == Set([:pullup])
        @test sl.governing_root == :pullup
        @test sl.governing_break == :pullup
        @test sl.pullup.n_computed == 3.0
        @test sl.governing_strut == :none
        pullup = sl.results[:pullup]
        @test abs(pullup.M_root) ≈ abs(ac.wing.inboard.moment) rtol = 1e-3
        @test abs(pullup.M_break) ≈ abs(ac.wing.outboard.moment) rtol = 1e-3
        @test abs(pullup.S_root) ≈ abs(ac.wing.inboard.max_shear_load) rtol = 1e-3
        @test abs(pullup.S_break) ≈ abs(ac.wing.outboard.max_shear_load) rtol = 1e-3
    end

    # `load_model` builds a fresh aircraft, so that both sizings start from the same state
    function test_extra_cases_do_not_govern(load_model, base, key)
        extra = sized_once(key; iter = 50) do
            ac = load_model()
            ac.wing.sizing_loads.pushover.active = ac.wing.sizing_loads.pullup_light.active =
                ac.wing.sizing_loads.pushover_light.active = ac.wing.sizing_loads.taxi.active = true
            ac
        end
        sl = extra.wing.sizing_loads
        @test Set(keys(sl.results)) == Set((:pullup, extra_cases...))
        @test sl.governing_root == :pullup
        @test sl.governing_break == :pullup
        @test sl.governing_strut == :none
        @test extra.wing.weight == base.wing.weight
        @test extra.parg[igWMTO] == base.parg[igWMTO]
    end

    @testset "fuel in the wing" begin
        base = sized_default_aircraft()
        test_recorded_pullup(base)
        test_extra_cases_do_not_govern(load_default_model, base, :default_extra_cases)
    end

    @testset "fuel in the fuselage" begin
        base = sized_cryo_aircraft()
        test_recorded_pullup(base)
        test_extra_cases_do_not_govern(() -> read_aircraft_model(cryo_input_path()), base, :cryo_extra_cases)
    end

    @testset "has no minimum gauge, strut, jury or strut-attach results" begin
        ac = sized_default_aircraft()
        strut, sl = ac.wing.strut, ac.wing.sizing_loads
        for section in (ac.wing.inboard, ac.wing.outboard)
            @test section.min_cap_thickness == 0
            @test section.min_web_thickness == 0
        end
        @test strut.axial_force_compression == 0
        @test strut.weight == 0
        @test strut.dxW == 0
        @test sl.governing_strut == :none
        @test sl.governing_strut_compression == :none
        @test strut.EI_min == 0
        @test strut.buckling_load == 0
        @test strut.shear_load == 0
        @test strut.bending_moment == 0
        @test strut.sizing_criterion === :none
        for field in (:eta_wing, :length, :chord, :area, :axial_force, :buckling_load, :weight, :dxW)
            @test getfield(ac.wing.jury, field) == 0
        end
        @test !hasproperty(ac.htail, :strut)
        @test !hasproperty(ac.vtail, :strut)
        @test ac.fuselage.weight_frac_strut_attach == 0.0
        @test ac.fuselage.strut_attach.W == 0
        @test ac.fuselage.strut_attach_criterion === :none
        @test strut.weight_frac_attachments == 0.06 # the struct default adds nothing: the strut weight is zero
        @test strut.reaction_vertical == 0
    end

    @testset "sizing responds to the pull-up load factor" begin
        lighter = sized_once(:pullup_2p5) do
            ac = load_default_model()
            ac.wing.sizing_loads.pullup.load_factor = 2.5
            ac
        end
        @test lighter.wing.weight < sized_default_aircraft().wing.weight
    end

    @testset "FAR load factor option" begin
        far_manoeuvre_load_factor = TASOPT.structures.far_manoeuvre_load_factor
        base = sized_default_aircraft()

        # Option off with the pull-up load factor given explicitly: identical to the default sizing
        explicit = sized_once(:far_off_explicit_load_factor; iter = 50) do
            first(read_aircraft_from_default_dict() do d
                d["Options"]["far_load_factor"] = false
                sizing_load_table(d, :pullup)["load_factor"] = 3.0
            end)
        end
        @test !explicit.options.opt_far_load_factor
        @test explicit.parg[igWMTO] == base.parg[igWMTO]
        @test explicit.wing.weight == base.wing.weight
        @test explicit.wing.sizing_loads.pullup.n_computed == 3.0

        # Option on, with the positive and negative manoeuvre cases and the taxi case active
        far = sized_once(:far_load_factor; iter = 50) do
            ac = first(read_aircraft_from_default_dict(d -> d["Options"]["far_load_factor"] = true))
            sl = ac.wing.sizing_loads
            sl.pushover.active = sl.pullup_light.active = sl.pushover_light.active = sl.taxi.active = true
            ac
        end
        sl = far.wing.sizing_loads
        n_far = far_manoeuvre_load_factor(far.parg[igWMTO])

        # 197000 lb: 2.1 + 24000/207000 = 2.216 -> 2.5, below the default 3.0
        @test n_far == 2.5
        @test sl.pullup.n_computed ≈ n_far rtol = 1e-6
        @test sl.pullup_light.n_computed ≈ n_far rtol = 1e-6

        # the user inputs are kept, and the negative and ground cases are unchanged
        @test sl.pullup.load_factor == 3.0
        @test sl.pullup_light.load_factor == 3.0
        @test sl.pushover.n_computed == -1.0
        @test sl.pushover_light.n_computed == -1.0
        @test sl.taxi.n_computed == -2.0

        # the lower load factor gives a lighter wing
        @test far.wing.weight < base.wing.weight
    end
end

# ============================================================================
# Edge cases and guards: infeasible sizing inputs and tiny loads
# ============================================================================

# Smoke test: sizing the default aircraft at high wing aspect ratio must either
# converge, stop with a "max iterations" warning, or fail with a specific diagnostic.
# Bare DomainErrors and other uncontrolled crashes are regressions.

# Messages of the guards that are allowed to stop sizing
const clean_failures = (
    ArgumentError => r"negative tail area",
    ErrorException => r"cap sizing infeasible",
)

"Size the default aircraft at aspect ratio `AR`; return `(outcome, detail)`."
function size_at_aspect_ratio(AR)
    ac = CoreLogging.with_logger(CoreLogging.NullLogger()) do
        load_default_model()
    end
    ac.wing.layout.AR = AR

    logger = Test.TestLogger(min_level = CoreLogging.Warn)
    try
        CoreLogging.with_logger(logger) do
            size_aircraft!(ac, printiter = false)
        end
    catch e
        msg = sprint(showerror, e)
        for (type, pattern) in clean_failures
            e isa type && occursin(pattern, msg) && return (:clean_error, msg)
        end
        return (:uncontrolled_error, msg)
    end
    hit_itermax = any(r -> occursin("max iterations", r.message), logger.logs)
    return (hit_itermax ? :itermax_warning : :converged, "")
end

@testset "TBW aspect ratio smoke" verbose=true begin
    for AR in (10, 12, 14, 16, 18, 20)
        @testset "AR = $AR" begin
            outcome, detail = size_at_aspect_ratio(AR)
            @info "AR = $AR: $outcome" detail
            @test outcome in (:converged, :itermax_warning, :clean_error)
            outcome == :uncontrolled_error && @error "Uncontrolled failure at AR = $AR" detail
        end
    end
end

@testset "size_member_box sets the minimum gauges when the loads are tiny" verbose=true begin
    st = TASOPT.structures
    material = TASOPT.StructuralAlloy("TASOPT-Al")
    caps = st.StructuralMember(material = material)
    webs = st.StructuralMember(material = TASOPT.StructuralAlloy("TASOPT-Al"))
    E, σcap, τweb = material.E, material.σmax, material.τmax
    fu = material.safety_factor
    sigfac = 1.1

    # box of chord 2 m, so that the chord enters the stiffness as c^4 and the loads as c^2
    c, h, w, web_ratio = 2.0, 0.15, 0.5, 0.8
    hw = web_ratio * h
    min_cap, min_web = 1.5e-3, 1.0e-3

    box(; kwargs...) = st.size_member_box(caps, webs;
        chord = c, height_to_chord = h, width_to_chord = w, web_to_box_height = web_ratio,
        tension = 0.0, compression = 0.0, buckling_length = 16.0, K = 1.0, shear_load = 0.0, sigfac = sigfac,
        min_cap_thickness = min_cap, min_web_thickness = min_web, kwargs...)

    EI_xx(tbcap, tbweb) = E * c^4 * w * (h^3 - (h - 2tbcap)^3) / 12 + E * c^4 * tbweb * hw^3 / 6
    EI_yy(tbcap, tbweb) = E * c^4 * tbcap * w^3 / 6 + E * c^4 * tbweb * hw * w^2 / 2

    @testset "the minimum gauges set the box when the loads are tiny" begin
        r = box(tension = 10.0, compression = -10.0, shear_load = 1.0, buckling_length = 1.0)
        @test r.tbcap == min_cap / c
        @test r.tbweb == min_web / c
        @test r.criterion === :gauge
        @test r.area_caps ≈ 2 * (min_cap / c) * w * c^2 rtol = 1e-12
        @test r.area_webs ≈ 2 * (min_web / c) * hw * c^2 rtol = 1e-12
        @test r.EI_min ≈ min(EI_xx(min_cap / c, min_web / c), EI_yy(min_cap / c, min_web / c)) rtol = 1e-10
    end
end

# ============================================================================
# Member and box sizing (unit tests of size_member_box)
# ============================================================================

# Strut and jury strut as a spar box: caps carry the axial force (yield and Euler buckling), webs carry the shear.
# Expected values are written out from the box formulas, not taken from the sizing function.

@testset "size_member_box sizes the caps and webs of a strut box" verbose=true begin
    st = TASOPT.structures
    material = TASOPT.StructuralAlloy("TASOPT-Al")
    caps = st.StructuralMember(material = material)
    webs = st.StructuralMember(material = TASOPT.StructuralAlloy("TASOPT-Al"))
    E, σcap, τweb = material.E, material.σmax, material.τmax
    fu = material.safety_factor
    sigfac = 1.1

    # box of chord 2 m, so that the chord enters the stiffness as c^4 and the loads as c^2
    c, h, w, web_ratio = 2.0, 0.15, 0.5, 0.8
    hw = web_ratio * h
    min_cap, min_web = 1.5e-3, 1.0e-3

    box(; kwargs...) = st.size_member_box(caps, webs;
        chord = c, height_to_chord = h, width_to_chord = w, web_to_box_height = web_ratio,
        tension = 0.0, compression = 0.0, buckling_length = 16.0, K = 1.0, shear_load = 0.0, sigfac = sigfac,
        min_cap_thickness = min_cap, min_web_thickness = min_web, kwargs...)

    EI_xx(tbcap, tbweb) = E * c^4 * w * (h^3 - (h - 2tbcap)^3) / 12 + E * c^4 * tbweb * hw^3 / 6
    EI_yy(tbcap, tbweb) = E * c^4 * tbcap * w^3 / 6 + E * c^4 * tbweb * hw * w^2 / 2

    @testset "the webs are sized by the shear and the areas add up" begin
        shear = 5.0e5
        r = box(tension = 1.0e6, shear_load = shear)
        tbweb = TASOPT.structures.calc_web_thickness(τweb * sigfac, shear, c, hw)
        @test tbweb > min_web / c
        @test r.tbweb ≈ tbweb rtol = 1e-12
        @test r.area_webs ≈ 2 * tbweb * hw * c^2 rtol = 1e-12
        @test r.area ≈ r.area_caps + r.area_webs rtol = 1e-14

        # the webs carry no axial force: the cap thickness does not depend on the shear
        @test box(tension = 1.0e6).tbcap == r.tbcap
    end

    @testset "the caps give the box the Euler buckling stiffness f_u |P_c| (K l)^2 / pi^2 about the weaker axis" begin
        P, l, K = 1.0e6, 16.0, 1.0
        tbweb = min_web / c
        EI_required(P, K, l) = fu * P * (K * l)^2 / pi^2
        EI_req = EI_required(P, K, l)

        # flapwise bending (box height direction) is the weaker axis of this box
        r = box(compression = -P, buckling_length = l, K = K)
        @test r.criterion === :buckling
        @test r.buckling_axis === :xx
        @test r.tbweb == tbweb
        @test r.tbcap > min_cap / c
        @test EI_xx(r.tbcap, tbweb) ≈ EI_req rtol = 1e-10
        @test EI_yy(r.tbcap, tbweb) > EI_req
        @test r.EI_min ≈ EI_req rtol = 1e-10
        @test r.buckling_load ≈ pi^2 * EI_xx(r.tbcap, tbweb) / (K * l)^2 rtol = 1e-10
        @test r.buckling_load ≈ fu * P rtol = 1e-10
        @test r.area_caps ≈ 2 * r.tbcap * w * c^2 rtol = 1e-12
        # the caps are no thicker than needed
        @test EI_xx(0.99 * r.tbcap, tbweb) < EI_req

        # the required stiffness is proportional to P and to (K l)^2
        half_length = box(compression = -P, buckling_length = l / 2, K = K)
        half_stiffness = box(compression = -P, buckling_length = l, K = K / 2)
        half_load = box(compression = -P / 2, buckling_length = l, K = K)
        for r_scaled in (half_length, half_stiffness, half_load)
            @test r_scaled.criterion === :buckling
        end
        @test half_length.EI_min ≈ EI_required(P, K, l / 2) rtol = 1e-10
        @test half_length.EI_min ≈ r.EI_min / 4 rtol = 1e-10
        @test half_stiffness.EI_min ≈ r.EI_min / 4 rtol = 1e-10
        @test half_load.EI_min ≈ r.EI_min / 2 rtol = 1e-10
        @test half_length.tbcap < r.tbcap

        # a tall narrow box: chordwise bending (box width direction) is the weaker axis, and the caps follow in closed form
        h2, w2 = 0.3, 0.2
        hw2 = web_ratio * h2
        narrow = box(compression = -P, buckling_length = l, K = K, height_to_chord = h2, width_to_chord = w2)
        @test narrow.criterion === :buckling
        @test narrow.buckling_axis === :yy
        @test narrow.tbcap ≈ 6 * (EI_req - E * c^4 * tbweb * hw2 * w2^2 / 2) / (E * c^4 * w2^3) rtol = 1e-10
        @test narrow.EI_min ≈ EI_req rtol = 1e-10
    end
end

# ============================================================================
# Strut sizing (sized strut-braced and truss-braced aircraft)
# ============================================================================

@testset "strut box sizing" verbose=true begin
    st = TASOPT.structures

    # Sized strut-braced aircraft with the compatibility strut force, shared between the testsets by `key`.
    size_braced(key = :compatibility, edit! = ac -> nothing) = sized_braced_aircraft(key, edit!)
    wing_chord_at_strut(ac) = ac.wing.layout.root_chord * ac.wing.inboard.λ
    # a strut with an end-fixity factor of 2 buckles before it yields; no strut-attach reinforcement of the fuselage.
    # Its chord is the wing chord: at the default chord ratio no box of this height can reach the buckling stiffness
    function slender!(ac)
        ac.wing.strut.end_fixity_K = 2.0
        ac.wing.strut.chord_to_wing_chord = 1.0
        ac.wing.strut.weight_frac_attachments = 0.0
    end

    @testset "the strut chord is the chord ratio times the wing chord at the strut attachment" begin
        for ac in (size_braced(), sized_jury_aircraft("strut_braced"), sized_jury_aircraft("truss_braced", 0.3),
                   sized_prescribed_pullup_aircraft())
            strut = ac.wing.strut
            @test strut.chord_to_wing_chord == 0.5
            @test strut.chord ≈ strut.chord_to_wing_chord * wing_chord_at_strut(ac) rtol = 1e-12
            @test strut.S ≈ 2 * strut.chord * strut.length rtol = 1e-10
        end
    end

    @testset "the strut weight is both wings of caps and webs with the added weight fraction, and the moment uses it" begin
        no_fraction = sized_without_added_weight_fraction()
        for (ac, frac) in ((size_braced(), 0.60), (sized_jury_aircraft("truss_braced", 0.3), 0.60),
                           (sized_prescribed_pullup_aircraft(), 0.60), (no_fraction, 0.0))
            strut = ac.wing.strut
            box = recomputed_strut_box(ac)
            @test strut.weight_frac_added == frac
            @test strut.area > 0
            @test strut.weight ≈ box_member_weight(strut, box; frac) rtol = 1e-6
            @test strut.dxW ≈ strut.weight * strut_dx_arm(ac) rtol = 1e-6
            # neither the area nor the box contains the added fraction
            @test strut.area ≈ box.area_caps + box.area_webs rtol = 1e-6
        end
        with_fraction = sized_jury_aircraft("truss_braced", 0.3)
        @test with_fraction.wing.strut.weight >
              box_member_weight(with_fraction.wing.strut, recomputed_strut_box(with_fraction); frac = 0.0)
    end

    @testset "the own weight of the strut is a transverse load: shear on the webs and bending moment on the caps" begin
        for ac in (size_braced(), sized_jury_aircraft("strut_braced"), sized_jury_aircraft("truss_braced", 0.3),
                   size_braced(:prescribed, ac -> ac.wing.strut.opt_strut_force = TASOPT.StrutForceModel.Prescribed),
                   sized_prescribed_pullup_aircraft())
            strut, jury = ac.wing.strut, ac.wing.jury
            loads = strut_own_weight_loads(ac)
            @test strut.bending_moment > 0
            @test strut.bending_moment ≈ loads.M rtol = 1e-6
            @test strut.shear_load > 0
            @test strut.shear_load ≈ abs(jury.axial_force) + loads.V rtol = 1e-6
            @test strut.shear_load >= abs(jury.axial_force)
        end
    end

    # Probe values of the slender strut-braced aircraft (K = 2, chord ratio 1): tension area of the caps 6.3e-3 m^2,
    # total box area 2.2e-2 m^2 (the buckling cap thickness is proportional to K^2 for thin caps)
    @testset "a slender strut is governed by buckling" begin
        ac = size_braced(:slender, slender!)
        strut = ac.wing.strut
        box = recomputed_strut_box(ac)
        @test strut.sizing_criterion === :buckling
        @test box.criterion === :buckling
        @test box.area_caps * strut.caps.material.σmax * ac.parg[igsigfac] >
              1.1 * max(strut.axial_force, abs(strut.axial_force_compression))
        @test strut.buckling_load ≈ strut.caps.material.safety_factor * abs(strut.axial_force_compression) rtol = 1e-6

        # heavier caps make a heavier strut
        stocky = size_braced(:stocky, ac -> ac.wing.strut.end_fixity_K = 0.1).wing.strut
        @test strut.weight > stocky.weight
    end

    @testset "a stocky strut is governed by tension" begin
        ac = size_braced(:stocky, ac -> ac.wing.strut.end_fixity_K = 0.1)
        strut = ac.wing.strut
        box = recomputed_strut_box(ac)
        @test strut.sizing_criterion === :tension
        # the caps carry at least the tension (the bending moment of the own weight adds to it)
        @test box.area_caps ≥ strut.axial_force / (strut.caps.material.σmax * ac.parg[igsigfac]) * (1 - 1e-6)
        @test strut.buckling_load > strut.caps.material.safety_factor * abs(strut.axial_force_compression)
    end

    @testset "the compressive yield sets the caps when the compression exceeds the tension" begin
        # |P_c| is about 0.35 times the tension for the default sizing cases, so a heavier push-over case
        # (load factor -4) is used to make the compression larger
        ac = size_braced(:heavy_pushover, ac -> (ac.wing.strut.end_fixity_K = 0.1;
                                                 ac.wing.sizing_loads.pushover.load_factor = -4.0))
        strut = ac.wing.strut
        box = recomputed_strut_box(ac)
        σ_allowable = strut.caps.material.σmax * ac.parg[igsigfac]
        @test abs(strut.axial_force_compression) > strut.axial_force
        @test strut.sizing_criterion === :tension
        # the caps carry at least the compression (the bending moment of the own weight adds to it)
        @test box.area_caps ≥ abs(strut.axial_force_compression) / σ_allowable * (1 - 1e-6)
        @test box.area_caps > strut.axial_force / σ_allowable
    end

    @testset "Prescribed strut force is sized by tension only" begin
        pullup_only = sized_prescribed_pullup_aircraft()
        all_cases = size_braced(:prescribed, ac -> ac.wing.strut.opt_strut_force = TASOPT.StrutForceModel.Prescribed)
        for (ac, cases) in ((all_cases, [:pullup, :pushover, :pullup_light, :pushover_light, :taxi]),
                            (pullup_only, [:pullup]))
            strut, sl = ac.wing.strut, ac.wing.sizing_loads
            box = recomputed_strut_box(ac)
            @test Set(keys(sl.results)) == Set(cases)
            @test strut.axial_force > 0
            @test strut.axial_force_compression == 0
            @test strut.sizing_criterion in (:tension, :gauge)   # there is no buckling requirement without compression
            @test strut.area ≈ box.area rtol = 1e-6
            @test box.area_caps * strut.caps.material.σmax * ac.parg[igsigfac] >= strut.axial_force * (1 - 1e-6)
        end

        # the push-over, light and taxi cases are active and compressive, but they do not enter the strut
        @test all_cases.wing.sizing_loads.results[:pushover].strut_force < 0

        # with the pull-up as the only case: finite weights, the pull-up governs the strut and sets the own-weight load factor
        ac = pullup_only
        wing, strut, sl, layout = ac.wing, ac.wing.strut, ac.wing.sizing_loads, ac.wing.layout
        @test isfinite(wing.weight) && wing.weight > 0
        @test isfinite(ac.parg[igWMTO]) && ac.parg[igWMTO] > 0
        @test strut.weight > 0
        @test strut.area > 0
        @test strut.chord > 0
        @test sl.governing_strut == :pullup
        @test sl.results[:pullup].strut_force > 0
        @test sl.results[:pullup].P_axial_inboard < 0
        @test sl.pullup.n_computed ≈ sl.pullup.load_factor

        # the strut length is the diagonal over the height of the attachment and the inboard bay along the sweep
        length_3d = sqrt(strut.z^2 + (0.5 * layout.span * (layout.ηs - layout.ηo) / cosd(layout.sweep))^2)
        @test strut.length > 0
        @test strut.length ≈ length_3d rtol = 1e-6
    end

    @testset "minimum gauges set the strut box when they are large" begin
        ac = size_braced(:large_gauges, ac -> (ac.wing.strut.min_cap_thickness = 2.0e-2;
                                               ac.wing.strut.min_web_thickness = 2.0e-2))
        strut = ac.wing.strut
        hw = strut.web_to_box_height * strut.thickness_to_chord
        @test strut.sizing_criterion === :gauge
        # the caps and the webs are at their gauge: area = 2 t w c + 2 t hw c
        @test strut.area ≈ 2 * 2.0e-2 * strut.chord * (strut.box_width_to_chord + hw) rtol = 1e-6
        @test strut.area ≈ recomputed_strut_box(ac).area rtol = 1e-6
        # the default gauges leave a lighter strut
        @test strut.weight > size_braced().wing.strut.weight
    end

    @testset "a larger strut chord gives a lighter strut when buckling governs" begin
        # the same slender strut (K = 2, governed by buckling at the wing chord) with twice the wing chord
        wide = size_braced(:wide_chord, ac -> (slender!(ac); ac.wing.strut.chord_to_wing_chord = 2.0))
        slender = size_braced(:slender, slender!)
        @test slender.wing.strut.sizing_criterion === :buckling
        @test wide.wing.strut.chord ≈ 2.0 * wing_chord_at_strut(wide) rtol = 1e-12
        @test wide.wing.strut.chord > 1.5 * slender.wing.strut.chord
        @test wide.wing.strut.S > slender.wing.strut.S
        @test wide.wing.strut.area ≈ recomputed_strut_box(wide).area rtol = 1e-6
        # the second moment of the same box scales with the chord to the fourth power, so the caps can be much thinner
        @test wide.wing.strut.weight < slender.wing.strut.weight
    end
end

# ============================================================================
# Jury strut sizing
# ============================================================================

# Strut and jury weights are the weight of the caps and webs of both wings times (1 + weight_frac_added); the areas stay bare
@testset "jury strut sizing" verbose=true begin
    st = TASOPT.structures
    jury_load_fraction = 0.02   # jury load as a fraction of the strut compression
    jury_default_fields = (:eta_wing, :length, :chord, :area, :axial_force, :buckling_load, :weight, :dxW)

    @testset "strut_buckling_length is the longer part of the strut cut at the jury" begin
        L = 15.0
        @test st.strut_buckling_length(L, 0.5) ≈ 0.5 * L rtol = 1e-12
        @test st.strut_buckling_length(L, 0.3) ≈ 0.7 * L rtol = 1e-12
        @test st.strut_buckling_length(L, 0.8) ≈ 0.8 * L rtol = 1e-12
        @test st.strut_buckling_length(L, nothing) == L
    end

    truss = Dict(η => sized_jury_aircraft("truss_braced", η) for η in (0.3, 0.8))
    braced = sized_jury_aircraft("strut_braced")

    @testset "the main strut buckles over the longer part, whatever the jury position" begin
        for (η, ac) in truss
            strut = ac.wing.strut
            l_eff = max(η, 1 - η) * strut.length
            box = recomputed_strut_box(ac)
            @test strut.buckling_load ≈ pi^2 * strut.EI_min / (strut.end_fixity_K * l_eff)^2 rtol = 1e-6
            @test strut.EI_min ≈ box.EI_min rtol = 1e-6
            @test strut.buckling_load ≈ box.buckling_load rtol = 1e-6
            @test strut.buckling_load >= strut.caps.material.safety_factor * abs(strut.axial_force_compression) * (1 - 1e-6)
        end
    end

    @testset "the jury follows the strut and is sized as a pinned box under a fraction of the strut compression" begin
        for (η, ac) in truss
            strut, jury, layout = ac.wing.strut, ac.wing.jury, ac.wing.layout
            @test jury.eta_strut == η
            @test jury.length ≈ strut.z * (1 - η) rtol = 1e-8
            @test jury.eta_wing ≈ layout.ηo + η * (layout.ηs - layout.ηo) rtol = 1e-8
            @test jury.axial_force ≈ -jury_load_fraction * abs(strut.axial_force_compression) rtol = 1e-8

            # the jury chord is the chord ratio times the strut chord; the box ratios are those of the strut
            @test jury.chord_to_strut_chord == 0.5
            @test jury.chord ≈ jury.chord_to_strut_chord * strut.chord rtol = 1e-12

            box = recomputed_jury_box(ac)
            @test jury.area > 0
            @test jury.area ≈ box.area rtol = 1e-6
            @test jury.buckling_load ≈ box.buckling_load rtol = 1e-6
            @test jury.buckling_load >= jury.caps.material.safety_factor * abs(jury.axial_force) * (1 - 1e-6)
            if box.criterion === :buckling
                @test jury.buckling_load ≈ jury.caps.material.safety_factor * abs(jury.axial_force) rtol = 1e-6
            end
            @test jury.weight ≈ box_member_weight(jury, box) rtol = 1e-6
        end
    end

    @testset "jury weight is 2 g l (rho_cap A_caps + rho_web A_webs) (1 + weight_frac_added), and the moment uses the total weight" begin
        for (ac, frac) in ((sized_jury_aircraft("truss_braced", 0.3), 0.60), (sized_without_added_weight_fraction(), 0.0))
            jury = ac.wing.jury
            box = recomputed_jury_box(ac)
            @test jury.weight_frac_added == frac
            @test jury.area > 0
            @test jury.area ≈ box.area rtol = 1e-6
            @test jury.weight ≈ box_member_weight(jury, box; frac) rtol = 1e-6
            @test jury.dxW ≈ jury.weight * jury_dx_arm(ac) rtol = 1e-6
        end
    end

    @testset "size_jury! scales the weight and the moment exactly by (1 + weight_frac_added) and leaves the area alone" begin
        ac = sized_jury_aircraft("truss_braced", 0.3)
        etao, etas = ac.wing.layout.ηo, ac.wing.layout.ηs
        sinL, cosL = sind(ac.wing.layout.sweep), cosd(ac.wing.layout.sweep)
        function sized_jury(frac)
            wing = deepcopy(ac.wing)
            wing.jury.weight_frac_added = frac
            st.size_jury!(wing, etao, etas, ac.parg[igsigfac], sinL, cosL)
            return wing.jury
        end
        j0, j1, j2 = sized_jury(0.0), sized_jury(0.3), sized_jury(1.0)
        @test j0.weight ≈ box_member_weight(j0, recomputed_jury_box(ac); frac = 0.0) rtol = 1e-6
        @test j1.weight ≈ 1.3 * j0.weight rtol = 1e-12
        @test j2.weight ≈ 2.0 * j0.weight rtol = 1e-12
        @test j1.dxW ≈ 1.3 * j0.dxW rtol = 1e-12
        @test j2.dxW ≈ 2.0 * j0.dxW rtol = 1e-12
        @test j1.area == j0.area == j2.area
        @test j1.length == j0.length
        @test j1.chord == j0.chord
        @test j1.buckling_load == j0.buckling_load
    end

    @testset "a jury at the fraction eta_strut cuts the stiffness the strut needs against buckling to the square of the longer part" begin
        material = TASOPT.StructuralAlloy("TASOPT-Al")
        caps = st.StructuralMember(material = material)
        webs = st.StructuralMember(material = TASOPT.StructuralAlloy("TASOPT-Al"))
        P, L, K = 3.0e5, 30.0, 1.0
        strut_box(buckling_length) = st.size_member_box(caps, webs;
            chord = 2.0, height_to_chord = 0.15, width_to_chord = 0.5, web_to_box_height = 0.8,
            tension = 0.0, compression = -P, buckling_length, K, shear_load = 0.0, sigfac = 1.1,
            min_cap_thickness = 1.5e-3, min_web_thickness = 1.0e-3)
        EI_required(l) = material.safety_factor * P * (K * l)^2 / pi^2

        without_jury = strut_box(L)
        @test without_jury.criterion === :buckling
        @test without_jury.EI_min ≈ EI_required(L) rtol = 1e-10
        for (eta_strut, l_eff) in ((0.5, 0.5 * L), (0.3, 0.7 * L), (0.8, 0.8 * L))
            with_jury = strut_box(st.strut_buckling_length(L, eta_strut))
            @test with_jury.criterion === :buckling
            @test with_jury.EI_min ≈ EI_required(l_eff) rtol = 1e-10
            @test with_jury.EI_min / without_jury.EI_min ≈ (l_eff / L)^2 rtol = 1e-10
            @test with_jury.tbcap < without_jury.tbcap
            @test with_jury.area < without_jury.area
        end
        # a jury at mid-length quarters the required stiffness
        @test strut_box(st.strut_buckling_length(L, 0.5)).EI_min ≈ without_jury.EI_min / 4 rtol = 1e-10
    end

    @testset "a jury shortens the buckling length and so cannot make the strut heavier" begin
        # the same loads and the same box with the buckling length of the full strut instead of its longer part
        for (η, ac) in truss
            strut = ac.wing.strut
            with_jury = recomputed_strut_box(ac)
            without_jury = st.size_member_box(strut.caps, strut.webs;
                chord = strut.chord, height_to_chord = strut.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
                web_to_box_height = strut.web_to_box_height, tension = strut.axial_force,
                compression = strut.axial_force_compression, buckling_length = strut.length, K = strut.end_fixity_K,
                shear_load = strut.shear_load, bending_moment = strut.bending_moment, sigfac = ac.parg[igsigfac],
                min_cap_thickness = strut.min_cap_thickness, min_web_thickness = strut.min_web_thickness)
            @test with_jury.area <= without_jury.area * (1 + 1e-12)
            @test with_jury.EI_min <= without_jury.EI_min * (1 + 1e-12)
        end
        # a strut-braced wing has no jury
        jury = braced.wing.jury
        @test jury.eta_strut == 0.5
        for field in jury_default_fields
            @test getfield(jury, field) == 0
        end
        @test braced.wing.strut.buckling_load ≈ pi^2 * braced.wing.strut.EI_min /
              (braced.wing.strut.end_fixity_K * braced.wing.strut.length)^2 rtol = 1e-10
    end
end

# ============================================================================
# Wing sizing load cases (unit tests: mass states, root loading, envelope, braced section loads)
# ============================================================================

@testset "wing sizing load cases" verbose=true begin
    # Wing of the default aircraft, whose span is set here because it is only set when the aircraft is sized
    function default_wing_with_span()
        wing = CoreLogging.with_logger(CoreLogging.NullLogger()) do
            load_default_model().wing
        end
        wing.layout.span = 60.0
        return wing
    end
    # Section weights and moments are zero until the wing is sized
    function default_wing_with_section_weights()
        wing = default_wing_with_span()
        wing.inboard.weight, wing.inboard.dyW = 4.0e4, 6.0e4
        wing.outboard.weight, wing.outboard.dyW = 2.0e4, 8.0e4
        return wing
    end
    masses = (W_cruise = 5.0e5, W_MTO = 5.5e5, W_fuel = 1.0e5, W_reserve = 2.0e4,
              W_empty = 3.0e5, W_fuel_cap = 1.5e5, W_fuel_at_MTO = 8.0e4, r_fuel_design = 0.5)
    Lhtail = 2.0e4
    gammat, gammas, We, dyeinn = 0.5, 0.75, 1.0e4, 5.0

    @testset "mass states of the cases" verbose=true begin
        WingMassState = TASOPT.structures.WingMassState
        wing_mass_state = TASOPT.structures.wing_mass_state

        weights = (W_cruise = 900e3, W_MTO = 1000e3, W_fuel = 200e3, W_reserve = 30e3,
                   W_empty = 500e3, W_fuel_cap = 300e3, W_fuel_at_MTO = 150e3, r_fuel_design = 0.6)

        # expected (W_lift, r_fuel_wing, W_fuel_fuse) with fuel in the fuselage
        fuselage_fuel = (
            (WingMassState.DesignCruise, 900e3, 0.0, 200e3),
            (WingMassState.MTOW_MaxFuel, 1000e3, 0.0, 150e3),
            (WingMassState.MZFW, 830e3, 0.0, 30e3),
            (WingMassState.Light, 530e3, 0.0, 30e3),
        )
        # expected with fuel in the wing
        wing_fuel = (
            (WingMassState.DesignCruise, 900e3, 0.6, 0.0),
            (WingMassState.MTOW_MaxFuel, 1000e3, 0.5, 0.0),
            (WingMassState.MZFW, 830e3, 0.1, 0.0),
            (WingMassState.Light, 530e3, 0.1, 0.0),
        )

        @testset "fuel in the wing" begin
            for (state, W_lift, r_fuel_wing, W_fuel_fuse) in wing_fuel
                result = wing_mass_state(state, true; weights...)
                @test result.W_lift ≈ W_lift
                @test result.r_fuel_wing ≈ r_fuel_wing
                @test result.W_fuel_fuse == W_fuel_fuse
            end
        end

        @testset "fuel in the fuselage" begin
            for (state, W_lift, r_fuel_wing, W_fuel_fuse) in fuselage_fuel
                result = wing_mass_state(state, false; weights...)
                @test result.W_lift ≈ W_lift
                @test result.r_fuel_wing == r_fuel_wing
                @test result.W_fuel_fuse ≈ W_fuel_fuse
            end
        end

        @testset "no wing tank capacity without fuel in the wing, and the MTOW wing fuel fraction is capped at one" begin
            for state in instances(WingMassState.T)
                result = wing_mass_state(state, false; weights..., W_fuel_cap = 0.0)
                @test all(isfinite, values(result))
                @test result.r_fuel_wing == 0.0
            end
            result = wing_mass_state(WingMassState.MTOW_MaxFuel, true; weights..., W_fuel_at_MTO = 400e3)
            @test result.r_fuel_wing == 1.0
            @test result.W_fuel_fuse == 0.0
        end
    end

    @testset "masses of the cases" verbose=true begin
        # 1.7e5 = W_pay - W_paymax + W_fuel is the fuel that fits at MTOW with maximum payload
        sizing_masses(; W_fuel_cap = 3.0e5, has_wing_fuel = true) = TASOPT.wing_sizing_masses(;
            W_cruise = 9.0e5, W_MTO = 1.0e6, W_pay = 1.5e5, W_paymax = 2.0e5, W_fuel = 2.2e5,
            freserve = 0.1, W_fuel_cap, has_wing_fuel)

        @testset "fuel in the wing" begin
            m = sizing_masses()
            @test m.W_cruise == 9.0e5
            @test m.W_MTO == 1.0e6
            @test m.W_fuel == 2.2e5
            @test m.W_fuel_cap == 3.0e5
            @test m.W_reserve ≈ 2.0e4
            @test m.W_empty ≈ 6.3e5
            @test m.W_fuel_at_MTO ≈ 1.7e5
            @test m.r_fuel_design ≈ 1.7e5 / 3.0e5
        end

        @testset "the tank capacity caps the fuel at MTOW" begin
            m = sizing_masses(W_fuel_cap = 1.0e5)
            @test m.W_fuel_at_MTO == 1.0e5
            @test m.r_fuel_design ≈ 1.7
        end

        @testset "fuel in the fuselage" begin
            m = sizing_masses(W_fuel_cap = 0.0, has_wing_fuel = false)
            @test m.W_fuel_at_MTO ≈ 1.7e5
            @test m.r_fuel_design == 0.0
            @test m.W_fuel_cap == 0.0
            @test all(isfinite, values(m))
        end
    end

    @testset "root loading of the cases" verbose=true begin
        wing_case_root_loading = TASOPT.wing_case_root_loading
        wing_loading = TASOPT.wing_loading

        wing = default_wing_with_span()
        cases = wing.sizing_loads
        rclt, rcls = 1.0, 1.0
        W_lift, n_pullup = 5.0e5, 3.0
        root_loading(case) = wing_case_root_loading(case, wing, rclt, rcls, W_lift, Lhtail, n_pullup)

        @testset "root loading of the pull-up, push-over, ground and gust cases" begin
            # pull-up: wing_loading at the pull-up load factor
            loads = root_loading(cases.pullup)
            @test loads.Nload == 3.0
            @test loads.po == wing_loading(wing, rclt, rcls, n_pullup, W_lift, Lhtail)

            # push-over: the tail load scales with the load factor, and with tail_load_scale
            loads = root_loading(cases.pushover)
            @test loads.Nload == -1.0
            @test loads.po == wing_loading(wing, rclt, rcls, -1.0, W_lift, (-1.0 / n_pullup) * Lhtail)
            half = deepcopy(cases.pushover)
            half.tail_load_scale = 0.5
            @test root_loading(half).po ==
                  wing_loading(wing, rclt, rcls, -1.0, W_lift, 0.5 * (-1.0 / n_pullup) * Lhtail)
            @test root_loading(half).po != loads.po

            # ground: no wing lift
            loads = root_loading(cases.taxi)
            @test loads.Nload == -2.0
            @test loads.po == 0.0

            @test_throws ArgumentError root_loading(cases.gust_up_light)   # gust cases are not supported yet
        end
    end

    @testset "envelope of the cases" verbose=true begin
        Cantilever = TASOPT.structures.Cantilever
        wing_case_loads = TASOPT.wing_case_loads
        wing_case_root_loading = TASOPT.wing_case_root_loading

        base_wing = default_wing_with_section_weights()
        # sentinel section loads, which the evaluation must leave untouched
        base_wing.inboard.max_shear_load, base_wing.inboard.moment = 111.0, 222.0
        base_wing.outboard.max_shear_load, base_wing.outboard.moment = 333.0, 444.0

        rclt, rcls = 1.0, 1.0
        # lift weight of each mass state used below, from `masses`
        W_lift_cruise, W_lift_MTO, W_lift_light = 5.0e5, 5.5e5, 3.2e5

        evaluate(wing) = TASOPT.evaluate_wing_sizing_cases!(wing, masses, true, rclt, rcls, Lhtail,
                                                            gammat, gammas, We, dyeinn)
        # loads of a case computed directly: one inboard engine per wing, on a copy of the wing
        direct_loads(wing, Nload, po) = wing_case_loads(Cantilever(), deepcopy(wing), po, gammat, gammas,
                                                        Nload, We, 0, 0.0, 1, dyeinn)
        maneuver_loads(wing, case, W_lift) = begin
            Nload, po = wing_case_root_loading(case, wing, rclt, rcls, W_lift, Lhtail, 3.0)
            direct_loads(wing, Nload, po)
        end
        section_loads(wing) = (wing.inboard.max_shear_load, wing.inboard.moment,
                               wing.outboard.max_shear_load, wing.outboard.moment)
        same_loads(a, b) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))

        @testset "default: only the pull-up case" begin
            wing = deepcopy(base_wing)
            envelope = evaluate(wing)
            sl = wing.sizing_loads
            expected = maneuver_loads(wing, sl.pullup, W_lift_cruise)

            @test Set(keys(sl.results)) == Set([:pullup])
            @test same_loads(sl.results[:pullup], expected)
            @test sl.governing_root == :pullup
            @test sl.governing_break == :pullup
            @test envelope == (S_root = abs(expected.S_root), M_root = abs(expected.M_root),
                               S_break = abs(expected.S_break), M_break = abs(expected.M_break))
            @test sl.pullup.n_computed == 3.0
            @test section_loads(wing) == section_loads(base_wing)
        end

        @testset "maneuver and ground cases" begin
            wing = deepcopy(base_wing)
            sl = wing.sizing_loads
            for name in (:pushover, :pullup_light, :pushover_light, :taxi)
                getfield(sl, name).active = true
            end
            envelope = evaluate(wing)
            active = (:pullup, :pushover, :pullup_light, :pushover_light, :taxi)

            @test Set(keys(sl.results)) == Set(active)
            @test same_loads(sl.results[:pullup], maneuver_loads(wing, sl.pullup, W_lift_cruise))
            @test same_loads(sl.results[:pushover], maneuver_loads(wing, sl.pushover, W_lift_MTO))
            @test same_loads(sl.results[:pullup_light], maneuver_loads(wing, sl.pullup_light, W_lift_light))
            @test same_loads(sl.results[:pushover_light], maneuver_loads(wing, sl.pushover_light, W_lift_light))
            @test same_loads(sl.results[:taxi], direct_loads(wing, -2.0, 0.0))
            @test [getfield(sl, name).n_computed for name in active] == [3.0, -1.0, 3.0, -1.0, -2.0]
            @test section_loads(wing) == section_loads(base_wing)

            # the cases differ, so the governing cases are not trivially the first one
            M_root = [abs(sl.results[name].M_root) for name in active]
            M_break = [abs(sl.results[name].M_break) for name in active]
            @test length(unique(M_root)) == length(active)
            @test sl.governing_root == active[argmax(M_root)]
            @test sl.governing_break == active[argmax(M_break)]
            governing_root, governing_break = sl.results[sl.governing_root], sl.results[sl.governing_break]
            @test envelope == (S_root = abs(governing_root.S_root), M_root = abs(governing_root.M_root),
                               S_break = abs(governing_break.S_break), M_break = abs(governing_break.M_break))

            # evaluating again gives the same results
            results = deepcopy(sl.results)
            @test evaluate(wing) == envelope
            @test Set(keys(sl.results)) == Set(keys(results))
            @test all(same_loads(sl.results[name], results[name]) for name in active)
        end

        @testset "gust cases are skipped, stale results are removed, and no active case gives no envelope" begin
            wing = deepcopy(base_wing)
            sl = wing.sizing_loads
            sl.gust_up_light.active = true
            sl.results[:taxi] = TASOPT.structures.WingCaseLoads(M_root = 1.0e12)
            envelope = evaluate(wing)
            @test Set(keys(sl.results)) == Set([:pullup])
            @test sl.gust_up_light.n_computed == 0.0
            @test envelope.M_root == abs(sl.results[:pullup].M_root)

            sl.pullup.active = false
            @test evaluate(wing) === nothing
            @test isempty(sl.results)
        end
    end

    @testset "braced wing section loads" verbose=true begin
        st = TASOPT.structures
        Braced, Cantilever = st.Braced, st.Cantilever
        StrutForceModel = TASOPT.StrutForceModel
        wing_case_loads = TASOPT.wing_case_loads
        wing_case_root_loading = TASOPT.wing_case_root_loading

        base_wing = default_wing_with_section_weights()
        base_wing.strut.z = 154 * in_to_m
        base_wing.strut.area = 5.0e-3
        # stiffness of the previous sizing iteration
        base_wing.inboard.EI[1], base_wing.outboard.EI[1] = 5.0e9, 2.5e9

        span, ηs, ηo = base_wing.layout.span, base_wing.layout.ηs, base_wing.layout.ηo
        z = base_wing.strut.z
        L = 0.5 * span * (ηs - ηo)    # length of the inboard bay
        ls = sqrt(z^2 + L^2)          # strut length in the plane of the strut
        W_lift = (pullup = 5.0e5, pushover = 5.5e5, pullup_light = 3.2e5, pushover_light = 3.2e5, taxi = 5.5e5)

        "A copy of the base wing with the strut force model `model` and `edit!` applied."
        function wing_with(model; edit! = w -> nothing)
            wing = deepcopy(base_wing)
            wing.strut.opt_strut_force = model
            edit!(wing)
            return wing
        end
        "Root loading `(Nload, po)` of the named case"
        root_loading(wing, name) = wing_case_root_loading(getfield(wing.sizing_loads, name), wing,
                                                          1.0, 1.0, W_lift[name], Lhtail, 3.0)
        "Braced loads of the named case, computed on a copy of `wing`"
        function braced_loads(wing, name; neinn = 1)
            Nload, po = root_loading(wing, name)
            return wing_case_loads(Braced(), deepcopy(wing), po, gammat, gammas, Nload, We, 0, 0.0, neinn, dyeinn)
        end
        case_names = (:pullup, :pushover, :pullup_light, :pushover_light, :taxi)

        # A wing without inboard weight and engines, loaded by the lift only
        lift_only(model) = wing_with(model; edit! = w -> (w.inboard.weight = 0.0; w.inboard.dyW = 0.0))

        @testset "prescribed strut force makes the root moment equal the moment at the strut attachment" begin
            wing = lift_only(StrutForceModel.Prescribed)
            Nload, po = root_loading(wing, :pullup)
            loads = braced_loads(wing, :pullup; neinn = 0)

            R = (po * span / 12) * (ηs - ηo) * (1 + 2gammas) + loads.S_break
            @test loads.strut_force ≈ R * ls / z rtol = 1e-12
            @test loads.P_axial_inboard ≈ -loads.strut_force * L / ls rtol = 1e-12
            # the prescribed force makes the root moment equal the moment at the strut attachment
            @test loads.M_root ≈ loads.M_break rtol = 1e-9

            # outboard loads are those of the cantilever
            cantilever = wing_case_loads(Cantilever(), deepcopy(wing), po, gammat, gammas, Nload, We, 0, 0.0, 0, dyeinn)
            @test loads.S_break == cantilever.S_break
            @test loads.M_break == cantilever.M_break
        end

        @testset "compatibility strut force follows the strut_force_compatibility formula" begin
            wing = lift_only(StrutForceModel.Compatibility)
            Nload, po = root_loading(wing, :pullup)
            loads = braced_loads(wing, :pullup; neinn = 0)

            # unbraced bay: clamped at the root, linear spanload from po to po*gammas, outboard loads at the break
            outboard = wing_case_loads(Cantilever(), deepcopy(wing), po, gammat, gammas, Nload, We, 0, 0.0, 0, dyeinn)
            function M0(x)
                u = L - x
                a, c = po * (1 - (1 - gammas) * x / L), po * (1 - gammas) / L
                return outboard.M_break + outboard.S_break * u + a * u^2 / 2 - c * u^3 / 3
            end
            EI(x) = base_wing.inboard.EI[1] + (base_wing.outboard.EI[1] - base_wing.inboard.EI[1]) * x / L
            sinθ = z / ls
            f_strut = ls / (base_wing.strut.caps.material.E * base_wing.strut.area * sinθ^2)
            R = st.strut_force_compatibility(M0, EI, L, f_strut)

            @test loads.strut_force ≈ R * ls / z rtol = 1e-9
            @test loads.P_axial_inboard ≈ -loads.strut_force * L / ls rtol = 1e-9
            @test loads.M_root ≈ M0(0) - R * L rtol = 1e-9
            @test loads.S_root ≈ outboard.S_break + po * L * (1 + gammas) / 2 - R rtol = 1e-9
        end

        @testset "strut stiffness sets the share of the load the strut carries" begin
            pullup(area) = braced_loads(wing_with(StrutForceModel.Compatibility;
                                                  edit! = w -> w.strut.area = area), :pullup)
            flexible, nominal, stiff = pullup(1e-9), pullup(5e-3), pullup(1.0)

            @test 0 < abs(flexible.strut_force) < 1e-3 * nominal.strut_force
            @test 0 < nominal.strut_force < stiff.strut_force
            # M_root = M0(0) - R L falls as the strut takes more of the load
            @test stiff.M_root < nominal.M_root < flexible.M_root

            # a strut that carries nothing leaves the root loads of the cantilever
            wing = wing_with(StrutForceModel.Compatibility)
            Nload, po = root_loading(wing, :pullup)
            cantilever = wing_case_loads(Cantilever(), deepcopy(wing), po, gammat, gammas, Nload, We, 0, 0.0, 1, dyeinn)
            @test flexible.M_root ≈ cantilever.M_root rtol = 1e-4
            @test flexible.S_root ≈ cantilever.S_root rtol = 1e-4
        end

        @testset "evaluate_wing_sizing_cases! records the strut loads of a braced wing" begin
            evaluate(wing; kwargs...) = TASOPT.evaluate_wing_sizing_cases!(wing, masses, true, 1.0, 1.0, Lhtail,
                                                                           gammat, gammas, We, dyeinn; kwargs...)
            section_loads(w) = (w.inboard.max_shear_load, w.inboard.moment, w.outboard.max_shear_load, w.outboard.moment)

            for model in (StrutForceModel.Compatibility, StrutForceModel.Prescribed)
                wing = wing_with(model)
                sl = wing.sizing_loads
                for name in (:pushover, :pullup_light, :pushover_light, :taxi)
                    getfield(sl, name).active = true
                end
                sections_before = section_loads(wing)
                envelope = evaluate(wing; has_strut = true)
                results = sl.results

                @test Set(keys(results)) == Set(case_names)
                @test section_loads(wing) == sections_before
                forces = Dict(name => results[name].strut_force for name in case_names)
                # a pull-up loads the strut in tension, a push-over in compression
                @test forces[:pullup] > 0 && forces[:pushover] < 0
                @test results[:pullup].P_axial_inboard < 0 && results[:pushover].P_axial_inboard > 0
                @test all(abs(results[name].M_inboard_extreme) >= abs(results[name].M_root) for name in case_names)
                @test results[:pullup].P_axial_inboard ≈ -forces[:pullup] * L / ls rtol = 1e-12

                # the case with the largest tension governs the strut
                @test sl.governing_strut in case_names
                @test sl.governing_strut == case_names[argmax([forces[name] for name in case_names])]
                @test envelope.T_strut ≈ maximum(values(forces)) rtol = 1e-12

                # the case with the most negative strut force governs the compression of the strut
                @test envelope.T_strut_compression ≈ minimum(values(forces)) rtol = 1e-12
                @test envelope.T_strut_compression < 0
                @test sl.governing_strut_compression in case_names
                @test sl.governing_strut_compression == case_names[argmin([forces[name] for name in case_names])]
                @test abs(envelope.P_axial) ≈ maximum(abs(results[name].P_axial_inboard) for name in case_names) rtol = 1e-12

                # the inboard section is sized for the largest moment and shear over the cases
                @test envelope.M_root ≈ maximum(max(abs(results[name].M_root), abs(results[name].M_inboard_extreme))
                                                for name in case_names)
                @test envelope.S_root ≈ maximum(max(abs(results[name].S_root), abs(results[name].S_break))
                                                for name in case_names)
            end
        end

        @testset "a cantilever wing has no strut loads" begin
            wing = deepcopy(base_wing)
            envelope = TASOPT.evaluate_wing_sizing_cases!(wing, masses, true, 1.0, 1.0, Lhtail, gammat, gammas, We, dyeinn)
            @test wing.sizing_loads.governing_strut == :none
            @test wing.sizing_loads.results[:pullup].strut_force == 0.0
            @test wing.sizing_loads.results[:pullup].P_axial_inboard == 0.0
            @test !haskey(envelope, :T_strut)
            @test !haskey(envelope, :T_strut_compression)
            @test wing.sizing_loads.governing_strut_compression == :none
        end
    end
end

# ============================================================================
# Weight build-up and CG: jury strut and strut-attach reinforcement
# ============================================================================

@testset "jury strut weight in the aircraft weight build-up" verbose=true begin
    truss = sized_jury_aircraft("truss_braced", 0.3)

    # Copy of `ac` with `extra_weight` [N] added to the jury and its moment set for the same arm
    function with_extra_jury_weight(ac, extra_weight)
        heavier = deepcopy(ac)
        heavier.wing.jury.weight += extra_weight
        heavier.wing.jury.dxW = heavier.wing.jury.weight * jury_dx_arm(heavier)
        return heavier
    end

    @testset "the jury moment is its weight times the arm to its spanwise station on the swept wing" begin
        jury = truss.wing.jury
        @test jury.weight > 0
        @test jury.dxW ≈ jury.weight * jury_dx_arm(truss) rtol = 1e-6
        for ac in (sized_jury_aircraft("strut_braced"), sized_default_aircraft())
            @test ac.wing.jury.weight == 0
            @test ac.wing.jury.dxW == 0
        end
    end

    @testset "update_weights! scales the jury weight with WMTO like the strut weight" begin
        extra = 5e4
        heavier = with_extra_jury_weight(truss, extra)
        reference = deepcopy(truss)
        WMTO0 = truss.parg[igWMTO]
        TASOPT.update_weights!(reference, 1.0)
        TASOPT.update_weights!(heavier, 1.0)

        # the jury is part of the weight fraction closure: (Wpay + Wfuse) / WMTO = 1 - sum of fractions
        payload_and_fuselage = truss.parg[igWpay] + truss.fuselage.weight
        fraction_sum = 1 - payload_and_fuselage / reference.parg[igWMTO]
        @test heavier.parg[igWMTO] ≈ payload_and_fuselage / (1 - fraction_sum - extra / WMTO0) rtol = 1e-8
        @test heavier.parg[igWMTO] > reference.parg[igWMTO]

        # each component keeps its fraction of WMTO
        @test heavier.wing.jury.weight / heavier.parg[igWMTO] ≈ (truss.wing.jury.weight + extra) / WMTO0 rtol = 1e-8
        @test heavier.wing.strut.weight / heavier.parg[igWMTO] ≈ truss.wing.strut.weight / WMTO0 rtol = 1e-8
    end

    @testset "a heavier jury at its station moves the CG limits towards it, and size_htail accounts for it" begin
        extra = 1e5
        # jury at the middle of the strut: its station lies between the CG limits (at 0.3 it is forward of the forward limit)
        truss = sized_jury_aircraft("truss_braced", 0.5)
        heavier = with_extra_jury_weight(truss, extra)
        _, _, _, _, xCG_fwd, xCG_aft = TASOPT.CG_limits(truss)
        _, _, _, _, xCG_fwd_heavier, xCG_aft_heavier = TASOPT.CG_limits(heavier)

        # the extra weight acts at the wing box plus its moment arm, between the forward and aft limits
        x_jury = heavier.wing.layout.box_x + jury_dx_arm(heavier)
        @test xCG_fwd < x_jury < xCG_aft
        @test xCG_fwd < xCG_fwd_heavier < x_jury
        @test x_jury < xCG_aft_heavier < xCG_aft

        function size_tail!(ac)
            para = view(ac.para, :, :, 1)
            TASOPT.size_htail(ac, view(para, :, ipdescentn), view(para, :, ipcruise1), view(para, :, ipcruise1))
            return ac
        end
        base = size_tail!(deepcopy(truss))
        sized_tail = size_tail!(heavier)
        @test sized_tail.parg[igxCGfwd] ≈ xCG_fwd_heavier rtol = 1e-10
        @test sized_tail.parg[igxCGaft] ≈ xCG_aft_heavier rtol = 1e-10
        @test sized_tail.parg[igxCGfwd] > base.parg[igxCGfwd]
        @test abs(sized_tail.htail.layout.S - base.htail.layout.S) > 1e-3 * base.htail.layout.S
    end

    @testset "update_WMTO! adds the jury weight to the weight sum" begin
        extra = 5e4
        heavier = with_extra_jury_weight(truss, extra)
        reference = deepcopy(truss)
        TASOPT.update_WMTO!(reference, 1.0)
        TASOPT.update_WMTO!(heavier, 1.0)
        @test heavier.parg[igWMTO] - reference.parg[igWMTO] ≈ extra / (1 - truss.fuselage.HPE_sys.W) rtol = 1e-8
    end

    @testset "initialize_sizing_loop! starts the jury at zero weight" begin
        ac = deepcopy(truss)
        @test ac.wing.jury.weight > 0
        TASOPT.initialize_sizing_loop!(ac)
        @test ac.wing.jury.weight == 0
    end

    @testset "balance_aircraft! puts the jury weight and moment into the CG" begin
        extra = 1e5
        heavier = with_extra_jury_weight(truss, extra)

        # the main gear sits at the aft CG limit plus a fixed distance, and the limit moves with the jury;
        # keep the gear at the same x so that only the jury weight differs
        xcgB(ac) = TASOPT.CG_limits(ac)[6]
        heavier.landing_gear.main_gear.distance_CG_to_landing_gear += xcgB(truss) - xcgB(heavier)

        function cruise_xCG(ac)
            ac = deepcopy(ac)
            TASOPT.balance_aircraft!(ac, 1, ipcruise1, 0.5, 1.0, 0.0, TrimVar.None)
            return ac.para[iaxCG, ipcruise1, 1]
        end
        xCG, xCG_heavier = cruise_xCG(truss), cruise_xCG(heavier)

        # the extra weight acts at the wing box plus its moment arm, so the CG moves toward that point without passing it
        x_jury = heavier.wing.layout.box_x + jury_dx_arm(heavier)
        @test xCG_heavier != xCG
        @test sign(xCG_heavier - xCG) == sign(x_jury - xCG)
        @test abs(xCG_heavier - x_jury) < abs(xCG - x_jury)
    end
end

# Fuselage reinforcement at the strut attachment: a weight fraction of the strut weight of the previous sizing iteration
# (like the attachment hardware of the main wing), added to the fuselage weight and moment and lumped at the wing-box x.
# The vertical strut reaction and Fuselage weight_frac_strut_attach are currently unused: they do not change it.

@testset "fuselage strut-attach reinforcement" verbose=true begin
    @testset "sized braced aircraft: the reinforcement is the fraction of the strut weight, at the wing box" begin
        for ac in (sized_braced_aircraft(), sized_jury_aircraft("truss_braced"))
            reinforcement = ac.fuselage.strut_attach
            strut = ac.wing.strut
            @test strut.weight_frac_attachments == 0.06
            @test strut.weight > 0
            @test reinforcement.W > 0
            # the reinforcement is sized with the strut weight of the previous sizing iteration
            @test reinforcement.W ≈ strut.weight_frac_attachments * strut.weight rtol = 1e-3
            @test ac.fuselage.strut_attach_criterion === :strut_weight
            @test reinforcement.x ≈ ac.wing.layout.box_x rtol = 1e-3
            # the default weight_frac_strut_attach of 0.05 times the vertical reaction would be larger, and is not used
            @test ac.fuselage.weight_frac_strut_attach == 0.05
            @test strut.reaction_vertical > 0
            @test ac.fuselage.weight_frac_strut_attach * strut.reaction_vertical > reinforcement.W
        end
    end

    @testset "the reinforcement adds to the fuselage weight" begin
        with_fraction = sized_braced_aircraft()
        # the same aircraft without the reinforcement: no attachment fraction of the strut weight
        without = sized_braced_aircraft(:strut_attach_zero, ac -> (ac.wing.strut.weight_frac_attachments = 0.0))

        @test without.wing.strut.weight_frac_attachments == 0.0
        @test without.fuselage.strut_attach.W == 0
        @test without.fuselage.strut_attach_criterion === :none
        @test without.wing.strut.reaction_vertical > 0

        # the fuselage load-path feedback of the sizing loop changes the rest of the fuselage a little,
        # so the difference equals the reinforcement only roughly
        W_attach = with_fraction.fuselage.strut_attach.W
        @test W_attach > 0
        @test with_fraction.fuselage.weight > without.fuselage.weight
        @test with_fraction.fuselage.weight - without.fuselage.weight ≈ W_attach rtol = 0.5
        @test with_fraction.fuselage.moment > without.fuselage.moment
        @test with_fraction.parg[igWMTO] > without.parg[igWMTO]
    end

    @testset "weight_frac_strut_attach does not change the reinforcement or the sized aircraft" begin
        default = sized_braced_aircraft()
        unused = sized_braced_aircraft(:strut_attach_no_reaction_fraction, ac -> (ac.fuselage.weight_frac_strut_attach = 0.0))

        @test default.fuselage.weight_frac_strut_attach == 0.05
        @test unused.fuselage.weight_frac_strut_attach == 0.0
        @test unused.wing.strut.reaction_vertical > 0
        @test unused.fuselage.strut_attach_criterion === :strut_weight
        @test default.fuselage.strut_attach_criterion === :strut_weight
        @test unused.fuselage.strut_attach.W ≈ 0.06 * unused.wing.strut.weight rtol = 1e-3
        @test unused.fuselage.strut_attach.x ≈ unused.wing.layout.box_x rtol = 1e-3
        @test unused.fuselage.strut_attach.W ≈ default.fuselage.strut_attach.W rtol = 1e-10
        @test unused.fuselage.weight ≈ default.fuselage.weight rtol = 1e-10
        @test unused.fuselage.moment ≈ default.fuselage.moment rtol = 1e-10
        @test unused.parg[igWMTO] ≈ default.parg[igWMTO] rtol = 1e-10
    end
end

# ============================================================================
# Sized aircraft: governing sizing case and trend over the strut attachment station
# ============================================================================

# Switching sizing cases on and off changes the governing cases of a sized wing

@testset "governing case follows the active cases" verbose=true begin
    # The truss-braced default_tbw.toml has no cantilever load cases, and at AR = 16 a cantilever
    # wing cannot be sized, so the reference wing is made a cantilever with AR = 13.
    path = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw.toml")
    txt = replace(read(path, String), "wing_config = \"truss_braced\"" => "wing_config = \"cantilever\"",
                  "AR    = 16.0" => "AR    = 13.0")
    ac, _ = read_aircraft_from_text(txt)
    @test !ac.options.has_strut
    @test ac.wing.layout.AR == 13.0

    sl = ac.wing.sizing_loads
    function size_with_active(active)
        for name in (:pullup, :pushover)
            getfield(sl, name).active = name in active
        end
        CoreLogging.with_logger(CoreLogging.NullLogger()) do
            size_aircraft!(ac, iter = 50, printiter = false)
        end
        return (Set(keys(sl.results)), sl.governing_root, sl.governing_break)
    end

    @test size_with_active((:pullup,)) == (Set([:pullup]), :pullup, :pullup)
    # the push-over case (-1 g) has smaller loads than the pull-up (3 g), so it governs once alone
    @test size_with_active((:pushover,)) == (Set([:pushover]), :pushover, :pushover)
    @test size_with_active((:pullup, :pushover)) == (Set([:pullup, :pushover]), :pullup, :pullup)
end

# Structural model of the strut-braced wing over the strut attachment station ηs: the total structural weight has an
# interior optimum, the minimum gauge floors hold, and a slender strut is sized by buckling.

"Strut-braced aircraft with the strut attached at `ηs`, sized once and shared with the other strut tests."
function sized_strut_station_aircraft(ηs)
    ηs == 0.7 && return sized_jury_aircraft("strut_braced")
    ηs == 0.9 && return sized_braced_aircraft()
    return sized_once((:strut_station, ηs)) do
        ac = read_aircraft_with_wing_config("strut_braced")
        ac.wing.layout.ηs = ηs
        ac
    end
end

"Structural weight of the wing, strut and jury strut [N]."
structural_weight(ac) = ac.wing.weight + ac.wing.strut.weight + ac.wing.jury.weight

"Physical cap and web thickness [m] of a wing section, from the thickness normalised by the perpendicular chord."
function wing_section_gauges(wing, section)
    perpendicular_chord = section.co * cosd(wing.layout.sweep)
    return (cap = section.caps.thickness * perpendicular_chord,
            web = section.webs.thickness * perpendicular_chord)
end

const strut_stations = (0.5, 0.7, 0.9)

@testset "strut-braced wing over the strut attachment station" verbose=true begin
    low, middle, high = strut_stations
    aircraft = Dict(ηs => sized_strut_station_aircraft(ηs) for ηs in strut_stations)

    @testset "total structural weight has an interior optimum" begin
        weight = Dict(ηs => structural_weight(ac) for (ηs, ac) in aircraft)
        @test all(w -> isfinite(w) && w > 0, values(weight))
        @test weight[middle] < weight[low]
        @test weight[middle] < weight[high]
    end

    @testset "wing sections respect the minimum gauge" begin
        for ac in values(aircraft), section in (ac.wing.inboard, ac.wing.outboard)
            @test section.min_cap_thickness > 0
            @test section.min_web_thickness > 0
            gauges = wing_section_gauges(ac.wing, section)
            @test gauges.cap ≥ section.min_cap_thickness * (1 - 1e-9)
            @test gauges.web ≥ section.min_web_thickness * (1 - 1e-9)
        end
    end

    @testset "the slender strut is sized by buckling, the short one by tension" begin
        @test aircraft[low].wing.strut.sizing_criterion === :tension
        @test aircraft[middle].wing.strut.sizing_criterion === :buckling
        @test aircraft[high].wing.strut.sizing_criterion === :buckling
    end
end

# ============================================================================
# Weight build-up output: the empty-weight subtotal closes with WMTO - Wfuel - Wpay
# ============================================================================

"Text printed by `weight_buildup` for `ac`."
weight_buildup_text(ac) = sprint(io -> TASOPT.weight_buildup(ac; io = io))

"Weight [N] printed on the line of the weight build-up `text` that starts with `label` and the sign `+` or `=`, or `nothing` if there is no such line."
function printed_weight(text, label, sign = '+')
    m = match(Regex("^" * label * "\\s+\\" * sign * "\\s+([0-9.]+) N", "m"), text)
    return m === nothing ? nothing : parse(Float64, m.captures[1])
end

@testset "Weight build-up empty weight" verbose=true begin
    # Printed values are rounded to 0.1 N; a missing line is compared as NaN and so fails
    printing_tolerance = 0.1

    for (name, ac, has_strut, has_jury) in (
            ("cantilever",   sized_default_aircraft(),                  false, false),
            ("strut_braced", sized_jury_aircraft("strut_braced"),       true,  false),
            ("truss_braced", sized_jury_aircraft("truss_braced", 0.5),  true,  true))
        @testset "$name" begin
            text = weight_buildup_text(ac)
            parg = ac.parg

            # The subtotal printed after the components equals the empty weight of the sized aircraft
            Wempty_sized = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
            @test something(printed_weight(text, "Wempty", '='), NaN) ≈ Wempty_sized atol = printing_tolerance

            @test something(printed_weight(text, "Weng"), NaN) ≈ parg[igWeng] atol = printing_tolerance

            # The strut and jury lines are printed only for the wing configurations that have these members
            if has_strut
                @test something(printed_weight(text, "Wstrut"), NaN) ≈ ac.wing.strut.weight atol = printing_tolerance
            else
                @test printed_weight(text, "Wstrut") === nothing
            end
            if has_jury
                @test something(printed_weight(text, "Wjury"), NaN) ≈ ac.wing.jury.weight atol = printing_tolerance
            else
                @test printed_weight(text, "Wjury") === nothing
            end
        end
    end
end

# ============================================================================
# Inboard shear and moment distribution of a sizing load case, with the jump of the strut reaction
# ============================================================================

@testset "inboard load distribution of the sizing load cases" verbose=true begin
    @testset "cantilever wing stores no inboard distribution" begin
        sl = sized_default_aircraft().wing.sizing_loads
        @test !isempty(sl.results)
        for r in values(sl.results)
            @test isempty(r.inboard_stations)
            @test isempty(r.inboard_shear)
            @test isempty(r.inboard_moment)
            @test r.strut_reaction == 0.0
        end
    end

    for (name, ac) in (("strut_braced", sized_jury_aircraft("strut_braced")),
                       ("truss_braced", sized_jury_aircraft("truss_braced", 0.3)))
        @testset "$name governing strut case" begin
            wing = ac.wing
            sl = wing.sizing_loads
            @test sl.governing_strut != :none
            r = sl.results[sl.governing_strut]

            L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
            z = wing.strut.z
            ls = sqrt(z^2 + L^2)

            @test length(r.inboard_stations) == 11
            @test length(r.inboard_shear) == 11
            @test length(r.inboard_moment) == 11
            @test r.inboard_stations[1] == 0.0
            @test r.inboard_stations[end] ≈ L rtol = 1e-9

            # the distribution starts at the sized root loads
            @test r.inboard_shear[1] ≈ r.S_root rtol = 1e-6
            @test r.inboard_moment[1] ≈ r.M_root rtol = 1e-6

            # the moment is continuous at the strut attachment, the shear jumps there by the strut reaction
            @test r.inboard_moment[end] ≈ r.M_break rtol = 1e-6
            @test r.strut_reaction ≈ r.strut_force * z / ls rtol = 1e-9
            @test r.strut_reaction != 0.0
            @test r.inboard_shear[end] + r.strut_reaction ≈ r.S_break rtol = 1e-6

            # the stored extreme moment is the moment of largest magnitude of the distribution
            @test r.M_inboard_extreme == r.inboard_moment[argmax(abs.(r.inboard_moment))]
        end
    end
end
