const CoreLogging = Base.CoreLogging

# --- Helpers ---
"Default input text without the wing configuration lines."
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

"Aircraft built by `make()` and sized on the first request for `key`, then shared (do not modify; `deepcopy` first)."
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

"Box of the strut recomputed with `size_member_box` from the loads, geometry and inputs recorded on the sized aircraft `ac`."
function recomputed_strut_box(ac)
    wing, strut = ac.wing, ac.wing.strut
    has_jury = ac.options.opt_wing_config == WingConfig.TrussBraced
    return TASOPT.structures.size_member_box(strut.caps, strut.webs;
        chord = strut.chord_to_wing_chord * wing.layout.root_chord * wing.inboard.λ * strut.cos_lambda,
        height_to_chord = strut.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
        web_to_box_height = strut.web_to_box_height,
        tension = strut.axial_force, compression = strut.axial_force_compression,
        buckling_length = TASOPT.structures.strut_buckling_length(strut.length, has_jury ? wing.jury.eta_strut : nothing),
        K = strut.end_fixity_K, shear_load = strut.shear_load, bending_moment = strut.bending_moment,
        sigfac = ac.parg[igsigfac],
        min_cap_thickness = strut.min_cap_thickness, min_web_thickness = strut.min_web_thickness)
end

"Box of the jury strut recomputed with `size_member_box`: a pinned member under its recorded compression."
function recomputed_jury_box(ac)
    strut, jury = ac.wing.strut, ac.wing.jury
    return TASOPT.structures.size_member_box(jury.caps, jury.webs;
        chord = jury.chord_to_strut_chord * (strut.root_chord + (strut.chord - strut.root_chord) * jury.eta_strut),
        height_to_chord = jury.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
        web_to_box_height = strut.web_to_box_height,
        tension = 0.0, compression = jury.axial_force, buckling_length = jury.length, K = 1.0,
        shear_load = 0.0, sigfac = ac.parg[igsigfac],
        min_cap_thickness = jury.min_cap_thickness, min_web_thickness = jury.min_web_thickness)
end

"""
Shear `V` [N] and midspan moment `M` [N m] of the net transverse load of the strut of `ac`: in each case the strut weight
times the load factor and the strut lift act oppositely, uniformly along the strut; the largest net load per length acts
on a pinned beam over the buckling length (`V = q l_b / 2`, `M = q l_b^2 / 8`).
"""
function strut_net_transverse_loads(ac)
    wing, strut = ac.wing, ac.wing.strut
    sl = wing.sizing_loads
    L, z = strut.length, strut.z
    q_net(name, r) = abs(r.strut_lift - getfield(sl, name).n_computed * strut.weight / 2) / L * sqrt(L^2 - z^2) / L
    q_perp = maximum(q_net(name, r) for (name, r) in sl.results)
    has_jury = ac.options.opt_wing_config == WingConfig.TrussBraced
    l_b = TASOPT.structures.strut_buckling_length(L, has_jury ? wing.jury.eta_strut : nothing)
    return (V = q_perp * l_b / 2, M = q_perp * l_b^2 / 8)
end

"Weight [N] of the strut or jury strut `member` for both wings from the areas of its `box`, with the added weight fraction `frac`."
box_member_weight(member, box; frac = member.weight_frac_added) =
    2 * TASOPT.gee * member.length * (member.caps.material.ρ * box.area_caps + member.webs.material.ρ * box.area_webs) * (1 + frac)

# Strut-braced wing with ηs = 0.9 (where the Prescribed strut force sizes a nearly weightless wing), sized once per `key`; `edit!` changes it first.
function sized_braced_aircraft(key = :compatibility, edit! = ac -> nothing)
    return sized_once(key) do
        ac = read_aircraft_with_wing_config("strut_braced")
        ac.wing.layout.ηs = 0.9
        ac.wing.η_engine = 0.9 # the engine stays at the break
        edit!(ac)
        ac
    end
end

"Pitching-moment arm of the strut about the wing-box x."
strut_dx_arm(ac) = 0.25 * ac.wing.layout.span * (ac.wing.layout.ηs - ac.wing.layout.ηo) * tand(ac.wing.layout.sweep)

"Pitching-moment arm of the jury about the wing-box x."
jury_dx_arm(ac) = 0.5 * ac.wing.layout.span * (ac.wing.jury.eta_wing - ac.wing.layout.ηo) * tand(ac.wing.layout.sweep)

# Truss-braced wing like `sized_jury_aircraft("truss_braced", 0.3)` with a strut tapered to `taper_ratio`
function sized_tapered_strut_aircraft(taper_ratio)
    return sized_once((:tapered_strut, taper_ratio)) do
        ac = read_aircraft_with_wing_config("truss_braced")
        ac.wing.layout.ηs = 0.7
        ac.wing.jury.eta_strut = 0.3
        ac.wing.strut.taper_ratio = taper_ratio
        ac
    end
end

# Wing with the strut at ηs = 0.7 (the default ηs does not converge with a strut), sized once; `eta_strut` sets the jury
function sized_jury_aircraft(config, eta_strut = 0.5)
    return sized_once((:jury, config, eta_strut)) do
        ac = read_aircraft_with_wing_config(config)
        ac.wing.layout.ηs = 0.7
        config == "truss_braced" && (ac.wing.jury.eta_strut = eta_strut)
        ac
    end
end

# The sized `sized_jury_aircraft` (or `size_braced`) with the strut attach moved 3 m aft of its wing box, sized again
function sized_attach_aft_aircraft(config)
    return sized_once((:attach_aft, config)) do
        base = config == "strut_braced" ? sized_braced_aircraft() : sized_jury_aircraft(config, 0.3)
        ac = deepcopy(base)
        ac.wing.strut.x_strut_attach = base.wing.layout.box_x + 3.0
        ac
    end
end

# --- TOML input, defaults, and save and read round trips ---
@testset "TBW wing configuration option" verbose=true begin
    ac_default = CoreLogging.with_logger(CoreLogging.NullLogger()) do
        load_default_model()
    end

    @testset "wing_config key selects the wing configuration" begin
        @test ac_default.options.opt_wing_config == WingConfig.Cantilever

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

        both_txt = insert_below_header(has_strut_txt, "Options", "wing_config = \"cantilever\"") # wing_config wins
        @test first(read_aircraft_from_text(both_txt)).options.opt_wing_config == WingConfig.Cantilever
    end
end

# [Wing.Strut] input keys and the `ac.wing.strut` fields they populate
const strut_block_fields = (
    z_fuselage_attach = :z,
    thickness_to_chord = :thickness_to_chord,
    local_velocity_ratio = :local_velocity_ratio,
    box_width_to_chord = :box_width_to_chord,
    web_to_box_height = :web_to_box_height,
    end_fixity_K = :end_fixity_K,
    weight_frac_added = :weight_frac_added,
    weight_frac_attachments = :weight_frac_attachments,
    opt_strut_force = :opt_strut_force,
    strut_taper_ratio = :taper_ratio,
    lift_fraction_strut_to_total = :lift_fraction_strut_to_total,
)

# [Wing.Strut] body with values that differ from the defaults; distances carry units
const custom_strut_block = """
    z_fuselage_attach = "170 in"
    thickness_to_chord = 0.12
    local_velocity_ratio = 1.1
    box_width_to_chord = 0.45
    web_to_box_height = 0.7
    end_fixity_K = 0.8
    weight_frac_added = 0.5
    weight_frac_attachments = 0.08
    opt_strut_force = "prescribed"
    strut_taper_ratio = 0.6
    lift_fraction_strut_to_total = 0.4
"""

"Default input text for `config` with the [Wing.Strut] table replaced by `block` (no table if empty)."
function toml_with_strut_inputs(config; block = "")
    txt = replace(toml_with_wing_config(config),
        r"^[ \t]*(z_strut|strut_toc|strut_local_velocity_ratio)[ \t]*=.*\n"m => "",
        r"^\[Wing\.Strut\]\n(?:(?!\[).*\n)*"m => "")
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
        thickness_to_chord = 0.12,
        local_velocity_ratio = 1.1,
        box_width_to_chord = 0.45,
        web_to_box_height = 0.7,
        end_fixity_K = 0.8,
        weight_frac_added = 0.5,
        weight_frac_attachments = 0.08,
        opt_strut_force = TASOPT.StrutForceModel.Prescribed,
        strut_taper_ratio = 0.6,
        lift_fraction_strut_to_total = 0.4,
    )

    @testset "[Wing.Strut] populates the strut fields" begin
        ac = first(read_aircraft_from_text(toml_with_strut_inputs("truss_braced"; block = custom_strut_block)))
        for (key, field) in pairs(strut_block_fields)
            @test getfield(ac.wing.strut, field) == custom_values[key]
        end
        @test ac.wing.strut.z ≈ 170 * 0.0254
    end

    @testset "default input gives the default strut values" begin
        strut = first(read_aircraft_from_text(toml_with_wing_config("strut_braced"))).wing.strut
        @test strut.z ≈ 154 * 0.0254
        @test (strut.thickness_to_chord, strut.local_velocity_ratio, strut.box_width_to_chord) == (0.15, 1.0, 0.5)
        @test (strut.web_to_box_height, strut.end_fixity_K, strut.taper_ratio) == (0.75, 1.0, 1.0)
        @test (strut.weight_frac_added, strut.weight_frac_attachments) == (0.60, 0.06)
        @test strut.lift_fraction_strut_to_total == 0.25
        @test first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced"))).wing.jury.weight_frac_added == 0.60
    end

    @testset "save and read reproduce the strut inputs exactly" begin
        ac = first(read_aircraft_from_text(toml_with_strut_inputs("strut_braced"; block = custom_strut_block)))
        mktempdir() do dir
            path = joinpath(dir, "strut_roundtrip.toml")
            save_aircraft_model(ac, path)
            saved = TASOPT.TOML.parsefile(path)
            @test saved["Wing"]["Strut"]["opt_strut_force"] == "prescribed"
            for key in keys(strut_block_fields)
                @test haskey(saved["Wing"]["Strut"], string(key))
            end
            @test !haskey(saved["Wing"], "z_strut")

            ac_reread = first(read_aircraft_from_text(read(path, String)))
            for field in values(strut_block_fields)
                @test getfield(ac_reread.wing.strut, field) == getfield(ac.wing.strut, field)
            end
        end
    end

    @testset "spanload_exponent must be positive; the retired y_fuselage_attach key is read and ignored" begin
        @test_throws ErrorException read_aircraft_from_text(
            toml_with_strut_inputs("strut_braced"; block = "spanload_exponent = 0.0"))
        ac, logs = read_aircraft_from_text(toml_with_strut_inputs("strut_braced"; block = "y_fuselage_attach = \"1.5 ft\""))
        @test ac.wing.strut.spanload_exponent == 0.5
        @test any(r -> r.level == CoreLogging.Info && occursin("y_fuselage_attach", r.message), logs)
    end

    @testset "opt_strut_force defaults to compatibility, and an unknown model is rejected" begin
        @test first(read_aircraft_from_text(toml_with_wing_config("strut_braced"))).wing.strut.opt_strut_force ==
              TASOPT.StrutForceModel.Compatibility
        @test_throws Exception read_aircraft_from_text(
            toml_with_strut_inputs("strut_braced"; block = "opt_strut_force = \"rigid\""))
    end

    @testset "lift_fraction_strut_to_total must be in [0, 1), and the old mean_cl_ratio_to_wing key is rejected" begin
        read_strut_block(block) = read_aircraft_from_text(toml_with_strut_inputs("strut_braced"; block))
        for fraction in ("0.0", "0.999")
            @test first(read_strut_block("lift_fraction_strut_to_total = $fraction")).wing.strut.lift_fraction_strut_to_total ==
                  parse(Float64, fraction)
        end
        for fraction in ("1.0", "-0.1")
            @test_throws ErrorException read_strut_block("lift_fraction_strut_to_total = $fraction")
        end
        @test_throws ErrorException read_strut_block("mean_cl_ratio_to_wing = 0.5")
    end

    al_default = StructuralAlloy("TASOPT-Al", max_avg_stress = 1.1, safety_factor = 1.5)
    al_7075 = StructuralAlloy("Al-7075", max_avg_stress = 1.2, safety_factor = 2.0)
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
        @test jury_box_inputs(ac_strut.wing.jury) == default_jury_box # a strut-braced wing does not read [Wing.Jury]

        ac_truss = read_aircraft_with_box_inputs(set_box_inputs!, "truss_braced")
        @test strut_box_inputs(ac_truss.wing.strut) == explicit_strut_box
        @test jury_box_inputs(ac_truss.wing.jury) == explicit_jury_box
        @test alloy_signature(ac_truss.wing.jury.caps.material) == alloy_signature(ac_truss.wing.strut.caps.material)

        @test alloy_signature(ac_truss.wing.inboard.caps.material) == alloy_signature(al_default) # wing untouched

        # a single gauge key leaves the other at its default
        ac_one = read_aircraft_with_box_inputs("truss_braced") do d
            wing_subtable(d, "Strut")["min_cap_thickness"] = 0.002
            wing_subtable(d, "Jury")["min_web_thickness"] = 0.0004
        end
        @test (ac_one.wing.strut.min_cap_thickness, ac_one.wing.strut.min_web_thickness) == (0.002, 1.0e-3)
        @test (ac_one.wing.jury.min_cap_thickness, ac_one.wing.jury.min_web_thickness) == (1.5e-3, 0.0004)
    end

    @testset "the web keys fall back to the strut keys one by one" begin
        # no web keys: webs equal the caps, also for the jury
        ac = read_aircraft_with_box_inputs("truss_braced") do d
            d["Structures"]["strut_material"] = "Al-7075"
            d["Structures"]["strut_max_avg_stress"] = 1.2
            d["Structures"]["strut_safety_factor"] = 2.0
            for key in ("strut_web_material", "strut_web_max_avg_stress", "strut_web_safety_factor")
                delete!(d["Structures"], key)
            end
        end
        @test alloy_signature(ac.wing.strut.webs.material) == alloy_signature(al_7075)
        @test alloy_signature(ac.wing.jury.webs.material) == alloy_signature(al_7075)
    end

    @testset "[Wing.Jury] sets eta_strut of a truss-braced wing, and save and read reproduce it exactly" begin
        @test first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced"))).wing.jury.eta_strut == 0.5 # default
        @test first(read_aircraft_from_text(toml_with_jury_inputs("strut_braced", "eta_strut = 0.3"))).wing.jury.eta_strut == 0.5 # not read

        ac = first(read_aircraft_from_text(toml_with_jury_inputs("truss_braced", "eta_strut = 0.3")))
        mktempdir() do dir
            path = joinpath(dir, "jury_roundtrip.toml")
            save_aircraft_model(ac, path)
            saved = TASOPT.TOML.parsefile(path)
            @test saved["Wing"]["Jury"]["eta_strut"] == 0.3
            @test first(read_aircraft_from_text(read(path, String))).wing.jury.eta_strut == 0.3
        end
    end

    @testset "[Wing.Jury] thickness_to_chord: default, reading, validation, and save and read round trip" begin
        read_jury_block(config, block) = first(read_aircraft_from_text(toml_with_jury_inputs(config, block)))

        @test read_jury_block("truss_braced", "").wing.jury.thickness_to_chord == 0.10 # key absent
        @test read_jury_block("truss_braced", "eta_strut = 0.3").wing.jury.thickness_to_chord == 0.10
        @test read_jury_block("truss_braced", "thickness_to_chord = 0.12").wing.jury.thickness_to_chord == 0.12
        @test read_jury_block("strut_braced", "thickness_to_chord = 0.12").wing.jury.thickness_to_chord == 0.10 # not read

        # independent of the strut thickness
        ac = read_jury_block("truss_braced", "thickness_to_chord = 0.12")
        @test ac.wing.strut.thickness_to_chord == 0.15

        for thickness_to_chord in ("0.0", "-0.1", "1.0", "1.5")
            err = try
                read_jury_block("truss_braced", "thickness_to_chord = $thickness_to_chord")
                nothing
            catch e
                e
            end
            @test err isa ErrorException
            @test occursin("Input error", sprint(showerror, err))
        end

        mktempdir() do dir
            path = joinpath(dir, "jury_thickness_roundtrip.toml")
            save_aircraft_model(ac, path)
            @test TASOPT.TOML.parsefile(path)["Wing"]["Jury"]["thickness_to_chord"] == 0.12
            @test first(read_aircraft_from_text(read(path, String))).wing.jury.thickness_to_chord == 0.12
        end
    end
end

@testset "fuel in the strut: option and TOML input" verbose=true begin
    read_strut_fuel_aircraft(config, in_strut) = read_aircraft_with_box_inputs(d -> d["Fuel"]["fuel_in_strut"] = in_strut, config)
    error_message(f) = try
        f()
        ""
    catch err
        sprint(showerror, err)
    end

    @testset "default: no fuel in the strut" begin
        for config in ("cantilever", "strut_braced", "truss_braced")
            ac = read_aircraft_with_wing_config(config)
            @test ac.options.has_strut_fuel === false
            @test ac.wing.strut.fuel_capacity == 0.0
        end
    end

    @testset "fuel_in_strut key is read exactly, and save and read reproduce it" begin
        for in_strut in (true, false)
            ac = read_strut_fuel_aircraft("strut_braced", in_strut)
            @test ac.options.has_strut_fuel === in_strut
            @test ac.wing.strut.fuel_capacity == 0.0 # the capacity is computed when sizing
            mktempdir() do dir
                path = joinpath(dir, "strut_fuel_roundtrip.toml")
                save_aircraft_model(ac, path)
                @test TASOPT.TOML.parsefile(path)["Fuel"]["fuel_in_strut"] === in_strut
                @test first(read_aircraft_from_text(read(path, String))).options.has_strut_fuel === in_strut
            end
        end
    end

    @testset "fuel_in_strut = true needs wing fuel and a braced wing" begin
        # the fuselage-tank truss-braced input: the only default with a braced wing and `fuel_in_wing = false`
        lh2 = TASOPT.TOML.parsefile(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw_lh2.toml"))
        @test lh2["Fuel"]["fuel_in_wing"] == false
        lh2["Fuel"]["fuel_in_strut"] = true
        @test occursin("fuel_in_strut", error_message(() -> read_aircraft_from_text(sprint(TASOPT.TOML.print, lh2))))
        @test occursin("fuel_in_strut", error_message(() -> read_strut_fuel_aircraft("cantilever", true)))
    end
end

const sizing_load_case_names = (:pullup, :pushover, :pullup_light, :pushover_light,
                                :taxi, :gust_up_light, :gust_down_light)

"The [Wing.SizingLoads.`case`] table of the parsed input `d`, created if missing."
sizing_load_table(d, case) =
    get!(get!(d["Wing"], "SizingLoads", Dict{String,Any}()), string(case), Dict{String,Any}())

"The inputs of load case `name` that the TOML block carries."
case_inputs(loads, name) = (c = getfield(loads, name); (c.load_factor, c.mass_state, c.active, c.tail_load_scale))

mentions_nlift_deprecation(log) =
    log.level == CoreLogging.Info && occursin("Nlift", log.message) && occursin("deprecat", log.message)

@testset "wing sizing load case inputs" verbose=true begin
    @testset "defaults" begin
        WingMassState, WingLoadKind = TASOPT.structures.WingMassState, TASOPT.structures.WingLoadKind
        loads = TASOPT.structures.WingSizingLoads()

        expected = (pullup          = (WingLoadKind.Maneuver,  3.0, WingMassState.DesignCruise),
                    pushover        = (WingLoadKind.Maneuver, -1.0, WingMassState.MTOW_MaxFuel),
                    pullup_light    = (WingLoadKind.Maneuver,  3.0, WingMassState.Light),
                    pushover_light  = (WingLoadKind.Maneuver, -1.0, WingMassState.Light),
                    taxi            = (WingLoadKind.Ground,   -2.0, WingMassState.MTOW_MaxFuel),
                    gust_up_light   = (WingLoadKind.Gust,      1.0, WingMassState.Light),
                    gust_down_light = (WingLoadKind.Gust,     -1.0, WingMassState.Light))
        for (name, (kind, load_factor, mass_state)) in pairs(expected)
            case = getfield(loads, name)
            @test (case.kind, case.load_factor, case.mass_state) == (kind, load_factor, mass_state)
            @test case.active == (name == :pullup) # only the pull-up case is active
        end
        @test isempty(loads.results)
        @test (loads.governing_root, loads.governing_break, loads.governing_strut) == (:pullup, :pullup, :none)
    end

    @testset "TOML input" verbose=true begin
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
                for name in sizing_load_case_names
                    table = saved["Wing"]["SizingLoads"][string(name)]
                    @test all(haskey(table, key) for key in ("load_factor", "mass_state", "active", "tail_load_scale"))
                    @test !any(haskey(table, key) for key in ("kind", "n_computed"))
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
            @test any(mentions_nlift_deprecation, logs)
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
            ac, _ = read_aircraft_from_default_dict() do d
                d["Options"]["wing_config"] = "strut_braced"
                without_active_keys!(d)
            end
            @test active_cases(ac) == Set([:pullup, braced_extra_cases...])

            # a cantilever wing keeps only the pull-up case
            ac_cantilever, _ = read_aircraft_from_default_dict(without_active_keys!)
            @test active_cases(ac_cantilever) == Set([:pullup])
        end
    end
end

# FAR/CS 25.337: n = clamp(2.1 + 24000 / (W_lb + 10000), 2.5, 3.8), W_lb = W / 4.4482216152605
@testset "FAR positive manoeuvre load factor" begin
    n_far(W_lb) = TASOPT.structures.far_manoeuvre_load_factor(W_lb * 4.4482216152605)
    @test n_far(20000) ≈ 2.9 rtol = 1e-12 # 2.1 + 24000/30000
    @test n_far(4000) ≈ 3.8 # upper bound
    @test n_far(50000) ≈ 2.5 rtol = 1e-12 # 2.1 + 24000/60000, the lower bound
end

"Aircraft read from the default input for `config`, with the [Wing] minimum gauge keys removed and `edit!` applied to the Dict."
function read_aircraft_with_gauge_inputs(edit!, config)
    return first(read_aircraft_from_default_dict() do d
        d["Options"]["wing_config"] = config
        delete!(d["Wing"], "min_cap_thickness")
        delete!(d["Wing"], "min_web_thickness")
        edit!(d)
    end)
end

gauge_floors(section) = (section.min_cap_thickness, section.min_web_thickness)

@testset "TBW minimum wing gauge" verbose=true begin
    @testset "defaults depend on the wing configuration" begin
        for config in ("strut_braced", "truss_braced")
            ac = read_aircraft_with_gauge_inputs(d -> nothing, config)
            @test gauge_floors(ac.wing.inboard) == gauge_floors(ac.wing.outboard) == (1.5e-3, 1.0e-3)
            @test gauge_floors(ac.htail.inboard) == (0.0, 0.0) # tails have no floor
        end
        ac_cantilever = read_aircraft_with_gauge_inputs(d -> nothing, "cantilever")
        @test gauge_floors(ac_cantilever.wing.inboard) == gauge_floors(ac_cantilever.wing.outboard) == (0.0, 0.0)
    end

    @testset "explicit [Wing] keys set both sections and win over the defaults" begin
        ac = read_aircraft_with_gauge_inputs("truss_braced") do d
            d["Wing"]["min_cap_thickness"] = "0.1 in" # distance strings are allowed
            d["Wing"]["min_web_thickness"] = 0.0007
        end
        @test collect(gauge_floors(ac.wing.inboard)) ≈ [0.1 * 0.0254, 0.0007]
        @test collect(gauge_floors(ac.wing.outboard)) ≈ [0.1 * 0.0254, 0.0007]

        # an explicit zero switches the braced default off; the other key keeps its default
        ac_zero = read_aircraft_with_gauge_inputs("strut_braced") do d
            d["Wing"]["min_cap_thickness"] = 0.0
        end
        @test gauge_floors(ac_zero.wing.inboard) == gauge_floors(ac_zero.wing.outboard) == (0.0, 1.0e-3)
    end
end

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
            max_avg, safety = factors[name]
            @test (getter(ac).max_avg_stress, getter(ac).safety_factor) == (max_avg, safety)
            @test getter(ac).σmax ≈ getter(ac).YTS / max_avg / safety
        end

        mktempdir() do dir
            path = joinpath(dir, "material_factors_roundtrip.toml")
            save_aircraft_model(ac, path)
            saved = TASOPT.TOML.parsefile(path)["Structures"]
            ac_reread, _ = read_aircraft_from_text(read(path, String))
            for (name, getter) in pairs(material_getters)
                max_avg, safety = factors[name]
                @test (saved["$(name)_max_avg_stress"], saved["$(name)_safety_factor"]) == (max_avg, safety)
                @test (getter(ac_reread).max_avg_stress, getter(ac_reread).safety_factor) == (max_avg, safety)
            end
        end
    end
end

# Strut loads in the fuselage sizing: the side-wall patch, belly tie and pitch-bending relief of `fusew!`.
# `fusew!` arguments of the default fuselage (the values of the fuselage sizing test in unit_test_structures.jl)
const fusew_cbox = 3.119635043953234
const fusew_xwbox = 15.688410884477133
const fusew_xwing = 18.465770971525686
const fusew_args = (6.0, 219964.5779, 0.0, 0, 0.0, 0.0, 0.0, 0.0, "", 56016.1756316446, 14400.815502104786,
    9645.22120844135, 0.4, 0.7, 1.2379754637499652e6, 1.0469704911343623e6, 7.958844386073947,
    0.3, 1.0, 36.212475021199154, 35.05057234350367, fusew_xwing, fusew_xwbox, fusew_cbox, 14.164410884477133)

const fusew_default_fuselage = load_default_model().fuselage

"Copy of the default fuselage, sized by `fusew!` with the keywords `kwargs`."
function fuselage_sized_with(; kwargs...)
    fuse = deepcopy(fusew_default_fuselage)
    TASOPT.fusew!(fuse, fusew_args...; kwargs...)
    return fuse
end

# Strut geometry used to size the fuselage in the tests below
const strut_patch_z = 154 * 0.0254 # [m]
const strut_patch_y_root = 1.2 # [m]
const strut_patch_box_width = 1.0 # [m]

"""
Weight, x and governing criterion of the strut load-transfer patch (side walls and belly tie), written out independently of `fusew!`.
Each side wall is a trapezoid whose load-path width runs from `box_width` at the strut attach to `cbox` at the wing box (height `z`).
"""
function expected_strut_patch(fuse; force_vertical, force_horizontal, box_width = strut_patch_box_width,
                              x_attach = fusew_xwbox, xwbox = fusew_xwbox, cbox = fusew_cbox,
                              z = strut_patch_z, y_root = strut_patch_y_root)
    σ, ρ, E, ν = fuse.skin.σ, fuse.skin.ρ, fuse.material.E, 0.3
    Vw = force_vertical / 2
    n_intervals = 20
    dh = z / n_intervals
    load_width_yield = load_width_buckling = load_width_governing = 0.0
    for i in 1:n_intervals
        h = (i - 0.5) * dh
        l = box_width + (cbox - box_width) * h / z
        t_yield = Vw / (l * σ)
        t_buckle = (Vw * l * 12 * (1 - ν^2) / (4π^2 * E))^(1 / 3)
        load_width_yield += l * t_yield * dh
        load_width_buckling += l * t_buckle * dh
        load_width_governing += l * max(t_yield, t_buckle) * dh
    end
    W_walls = 2 * TASOPT.gee * ρ * load_width_governing
    W_tie = TASOPT.gee * ρ * force_horizontal / σ * 2 * y_root
    x = (W_walls * 0.5 * (xwbox + x_attach) + W_tie * x_attach) / (W_walls + W_tie)
    return (W = W_walls + W_tie, W_walls = W_walls, W_tie = W_tie, x = x,
            criterion = load_width_buckling > load_width_yield ? :wall_buckling : :wall_yield)
end

@testset "strut loads in the fuselage sizing" verbose=true begin
    strut_geometry = (strut_z = strut_patch_z, strut_y_root = strut_patch_y_root, strut_box_width = strut_patch_box_width)
    # With a strut box width of 1 m, yield governs the walls of `yield_loads` and buckling those of `buckling_loads`.
    yield_loads = (strut_force_vertical = 6.0e7, strut_force_horizontal = 1.5e7)
    buckling_loads = (strut_force_vertical = 1.0e6, strut_force_horizontal = 5.0e5)
    base = fuselage_sized_with()

    expected_patch(fuse, loads; kwargs...) = expected_strut_patch(fuse; force_vertical = loads.strut_force_vertical,
        force_horizontal = loads.strut_force_horizontal, kwargs...)

    @testset "without strut forces there is no patch" begin
        @test base.strut_patch.W == 0
        @test base.strut_patch_criterion === :none
        @test fuselage_sized_with(; strut_geometry...).weight == base.weight
    end

    @testset "the strut forces add a patch of the weight of the side walls and the belly tie" begin
        for (loads, criterion) in ((yield_loads, :wall_yield), (buckling_loads, :wall_buckling))
            fuse = fuselage_sized_with(; loads..., strut_geometry...)
            expected = expected_patch(fuse, loads)
            @test expected.criterion === criterion # the loads are chosen to be governed by the criterion
            @test fuse.strut_patch.W ≈ expected.W rtol = 1e-10
            @test fuse.strut_patch_criterion === criterion
            @test fuse.strut_patch.x ≈ fusew_xwbox rtol = 1e-12 # the strut attach follows the wing box by default
            @test fuse.weight ≈ base.weight + fuse.strut_patch.W rtol = 1e-12
            @test fuse.moment ≈ base.moment + fuse.strut_patch.W * fusew_xwbox rtol = 1e-12
        end

        # without the belly tie the wall weight scales with the force to the power 1 (yield) or 1/3 (buckling)
        walls_only(loads, factor) = fuselage_sized_with(; strut_force_vertical = factor * loads.strut_force_vertical,
            strut_force_horizontal = 0.0, strut_geometry...).strut_patch.W
        @test walls_only(yield_loads, 2) ≈ 2 * walls_only(yield_loads, 1) rtol = 1e-10
        @test walls_only(buckling_loads, 2) ≈ 2^(1 / 3) * walls_only(buckling_loads, 1) rtol = 1e-10
    end

    @testset "the patch is in addition to the strut-attachment hardware, and sits at the strut attach x" begin
        with_hardware = fuselage_sized_with(; strut_weight = 40e3, strut_weight_frac_attachments = 0.03, yield_loads..., strut_geometry...)
        without_hardware = fuselage_sized_with(; yield_loads..., strut_geometry...)
        @test with_hardware.strut_patch.W == without_hardware.strut_patch.W
        @test with_hardware.weight ≈ without_hardware.weight + 0.03 * 40e3 rtol = 1e-12

        # the walls sit midway between the wing box and the strut attach, the tie at the strut attach: the patch is at their weighted mean
        x_attach = fusew_xwbox + 2.0
        shifted = fuselage_sized_with(; yield_loads..., strut_geometry..., strut_x_attach = x_attach)
        expected = expected_patch(shifted, yield_loads; x_attach)
        @test shifted.strut_patch.W == without_hardware.strut_patch.W
        @test shifted.strut_patch.x ≈ expected.x rtol = 1e-12
        @test shifted.moment ≈ base.moment + shifted.strut_patch.W * expected.x rtol = 1e-12
    end

    @testset "the reaction share relieves the horizontal bending material only for a strut attach aft of the wing box" begin
        bending(; share, x_attach) = fuselage_sized_with(; yield_loads..., strut_geometry...,
            strut_reaction_share = share, strut_x_attach = x_attach).bendingmaterial_h.weight.W
        reference = base.bendingmaterial_h.weight.W
        @test bending(share = 0.5, x_attach = fusew_xwbox) == reference # attach at the wing box: bit-identical
        @test bending(share = 0.0, x_attach = fusew_xwbox + 3.0) == reference # no share: bit-identical
        @test bending(share = 0.5, x_attach = fusew_xwbox + 3.0) < reference # an aft attach relieves the bending material
    end

    @testset "sized strut-braced aircraft: the computed strut loads follow the strut force and geometry" begin
        ac = sized_braced_aircraft()
        wing, strut, sl = ac.wing, ac.wing.strut, ac.wing.sizing_loads
        # In-plane strut force (the axial forces are scaled to the total strut length by 1 / cos_lambda)
        L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
        ls = sqrt(strut.z^2 + L^2)
        T_plane = max(strut.axial_force, abs(strut.axial_force_compression)) * strut.cos_lambda
        @test strut.force_vertical ≈ 2 * T_plane * strut.z / ls rtol = 1e-9
        @test strut.force_horizontal ≈ T_plane * L / ls rtol = 1e-9

        # The share of the wing support taken by the strut, in the case with the largest strut force
        results = sl.results
        worst = argmax(name -> abs(results[name].strut_force), collect(keys(results)))
        share = results[worst].strut_reaction / (results[worst].S_root + results[worst].strut_reaction)
        @test 0 < strut.reaction_share < 1
        @test strut.reaction_share ≈ share rtol = 1e-2
        @test ac.fuselage.strut_patch.W > 0 # the fuselage is sized with the patch
        @test ac.fuselage.strut_patch_criterion in (:wall_yield, :wall_buckling)
    end

    # The sized fuselage is compared with the patch recomputed from the converged strut and wing of the aircraft
    sized_patch(ac; kwargs...) = expected_strut_patch(ac.fuselage;
        force_vertical = ac.wing.strut.force_vertical, force_horizontal = ac.wing.strut.force_horizontal,
        xwbox = ac.wing.layout.box_x, cbox = ac.wing.layout.root_chord * ac.wing.inboard.cross_section.width_to_chord,
        z = ac.wing.strut.z, y_root = ac.wing.layout.root_span / 2, kwargs...)

    @testset "sized tapered strut: the strut box width at the fuselage attach is the chord at the fuselage end" begin
        ac = sized_tapered_strut_aircraft(0.5)
        strut = ac.wing.strut
        expected = sized_patch(ac; box_width = strut.root_chord * strut.box_width_to_chord, x_attach = ac.wing.layout.box_x)
        @test ac.fuselage.strut_patch.W ≈ expected.W rtol = 1e-4
    end

    @testset "[Wing.Strut] x_strut_attach: absent is NaN, an explicit key is read exactly, save and read round-trip" begin
        x_attach(txt) = first(read_aircraft_from_text(txt)).wing.strut.x_strut_attach
        @test isnan(x_attach(toml_with_wing_config("strut_braced")))
        ac = first(read_aircraft_from_text(toml_with_strut_inputs("strut_braced"; block = "x_strut_attach = 21.234")))
        @test ac.wing.strut.x_strut_attach == 21.234
        mktempdir() do dir
            path = joinpath(dir, "strut_x_attach_roundtrip.toml")
            save_aircraft_model(ac, path)
            @test TASOPT.TOML.parsefile(path)["Wing"]["Strut"]["x_strut_attach"] == 21.234
            @test x_attach(read(path, String)) == 21.234
        end
    end
end

@testset "engine_location input: the key sets the engine location, and save and read reproduce it exactly" begin
    ac, _ = read_aircraft_from_default_dict() do d
        d["Wing"]["engine_location"] = 0.4321
    end
    @test ac.wing.η_engine == 0.4321
    mktempdir() do dir
        path = joinpath(dir, "engine_location_roundtrip.toml")
        save_aircraft_model(ac, path)
        @test TASOPT.TOML.parsefile(path)["Wing"]["engine_location"] == 0.4321
        @test first(read_aircraft_from_text(read(path, String))).wing.η_engine == 0.4321
    end
end

@testset "summary table of the evaluated sizing cases" begin
    WingCaseLoads = TASOPT.structures.WingCaseLoads
    ac = CoreLogging.with_logger(CoreLogging.NullLogger()) do
        load_default_model()
    end
    sl = ac.wing.sizing_loads
    io = IOBuffer()
    lines(text) = split(text, '\n')
    # first line that names the case, so that pullup does not match pullup_light
    case_row(text, name) = findfirst(line -> occursin(Regex("\\b$name\\b"), line), lines(text))
    # the leading digits of a value appear in the row, whatever the number format
    shows_value(row, x) = occursin(first(replace(string(round(abs(x), sigdigits = 3)), r"\D" => ""), 3),
                                   replace(row, r"[.,]" => ""))
    marked(row) = occursin(r"\*|root|break|governing"i, row)

    sl.pushover.n_computed, sl.taxi.n_computed = -1.0, -2.0
    # inserted out of order: the rows follow the order of the case fields
    sl.results[:taxi] = WingCaseLoads(S_root = 1.11e5, M_root = 2.22e6, S_break = 3.33e4, M_break = 4.44e5)
    sl.results[:pushover] = WingCaseLoads(S_root = 5.55e5, M_root = 6.66e6, S_break = 7.77e4, M_break = 8.88e5)
    sl.results[:pullup] = WingCaseLoads(S_root = 1.23e5, M_root = 4.56e6, S_break = 7.89e4, M_break = 9.87e5)
    sl.governing_root, sl.governing_break = :pushover, :taxi
    TASOPT.wing_sizing_loads_summary(ac; io = io)
    text = String(take!(io))
    rows = Dict(name => lines(text)[case_row(text, name)] for name in (:pullup, :pushover, :taxi))

    @test isnothing(case_row(text, :pullup_light))
    @test case_row(text, :pullup) < case_row(text, :pushover) < case_row(text, :taxi)
    for (name, loads) in sl.results
        @test all(shows_value(rows[name], getfield(loads, f)) for f in (:S_root, :M_root, :S_break, :M_break))
    end
    @test marked(rows[:pushover]) && marked(rows[:taxi]) && !marked(rows[:pullup])
end

@testset "geometry output of the strut tank" verbose=true begin
    geometry_text(ac) = sprint(io -> TASOPT.geometry(ac; io = io))
    printed(label, text) = (m = match(Regex("^" * label * raw"\s*=\s*([-0-9.]+)", "m"), text);
                            m === nothing ? NaN : parse(Float64, m.captures[1]))

    for config in ("strut_braced", "truss_braced")
        @testset "$config" begin
            ac = read_aircraft_with_wing_config(config)
            ac.options.has_strut_fuel = true
            ac.wing.strut.fuel_capacity = 2.0e5 # strut tank capacity [N]
            ac.parg[igWfmax] = 1.0e6            # total tank capacity [N]
            text_on = geometry_text(ac)
            ac.options.has_strut_fuel = false
            text_off = geometry_text(ac)

            @test printed("fuel capacity", text_on) ≈ 2.0e5 atol = 0.1
            @test printed("strut share of tank volume", text_on) ≈ 0.2 atol = 1e-3

            # switched off, the strut block is unchanged whatever the capacity field holds
            @test occursin("Strut Layout:", text_off)
            @test !occursin("fuel capacity", text_off)
            @test !occursin("strut share of tank volume", text_off)
            @test text_off == replace(text_on, r"^(fuel capacity|strut share of tank volume).*\n"m => "")
        end
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

# --- Golden values: hand-calculated reference values ---
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
end

@testset "gust load increment" verbose=true begin
    gust_load_increment = TASOPT.gust_load_increment
    U_de, V_EAS, CL_alpha, c_mean = 17.0688, 154.3, 4.5, 4.0   # m/s EAS (56 ft/s, 300 kt), 1/rad, m
    ρ0 = TASOPT.ρSL   # 1.2255627 kg/m³

    # Pratt–Walker: μ = 2(W/S) / (ρ0·c_mean·CLα·g),  K_g = 0.88μ/(5.3+μ),  Δn = ρ0·U_de·V·CLα·K_g / (2(W/S))
    @testset "hand calculation" begin
        # W/S = 5500 N/m²: μ = 11000 / (1.2255627·4·4.5·9.81) = 11000/216.40986 = 50.829476
        #   K_g = 0.88·50.829476/(5.3+50.829476) = 0.79690640
        #   Δn  = 1.2255627·17.0688·154.3·4.5·0.79690640 / 11000 = 14525.0275·0.79690640/11000 = 1.05228068
        @test gust_load_increment(5500.0, c_mean, CL_alpha, U_de, V_EAS) ≈ 1.0522806792964052 rtol = 1e-6
        # W/S = 3500 N/m²: μ = 7000/216.40986 = 32.346030, K_g = 0.75610912
        #   Δn  = 14525.0275·0.75610912/7000 = 1.56892939
        @test gust_load_increment(3500.0, c_mean, CL_alpha, U_de, V_EAS) ≈ 1.5689293899587433 rtol = 1e-6
    end
end

# --- Cantilever non-interference: braced-wing options do not change the cantilever wing ---
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
        test_recorded_pullup(sized_cryo_aircraft())
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
        @test strut.shear_load == 0
        @test strut.sizing_criterion === :none
        for field in (:eta_wing, :length, :chord, :area, :axial_force, :buckling_load, :weight, :dxW)
            @test getfield(ac.wing.jury, field) == 0
        end
        @test ac.fuselage.strut_attach.W == 0
        @test ac.fuselage.strut_attach_criterion === :none
    end

    @testset "FAR load factor option" begin
        far_manoeuvre_load_factor = TASOPT.structures.far_manoeuvre_load_factor
        base = sized_default_aircraft()

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

        # the user input is kept, and the negative and ground cases are unchanged
        @test sl.pullup.load_factor == 3.0
        @test sl.pushover.n_computed == -1.0
        @test sl.taxi.n_computed == -2.0

        # the lower load factor gives a lighter wing
        @test far.wing.weight < base.wing.weight
    end
end

# --- Edge cases and guards: infeasible sizing inputs and tiny loads ---
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
    for AR in (12, 16, 20)
        @testset "AR = $AR" begin
            outcome, detail = size_at_aspect_ratio(AR)
            @info "AR = $AR: $outcome" detail
            @test outcome in (:converged, :itermax_warning, :clean_error)
            outcome == :uncontrolled_error && @error "Uncontrolled failure at AR = $AR" detail
        end
    end
end

# --- Member and box sizing (unit tests of size_member_box) ---
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
        @test r.buckling_load ≈ fu * P rtol = 1e-10
        @test r.area_caps ≈ 2 * r.tbcap * w * c^2 rtol = 1e-12

        # the required stiffness is proportional to P and to (K l)^2
        half_length = box(compression = -P, buckling_length = l / 2, K = K)
        half_load = box(compression = -P / 2, buckling_length = l, K = K)
        @test half_length.EI_min ≈ r.EI_min / 4 rtol = 1e-10
        @test half_load.EI_min ≈ r.EI_min / 2 rtol = 1e-10

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

# --- Strut sizing (sized strut-braced and truss-braced aircraft) ---
@testset "strut box sizing" verbose=true begin
    st = TASOPT.structures

    # Sized strut-braced aircraft with the compatibility strut force, shared between the testsets by `key`.
    size_braced(key = :compatibility, edit! = ac -> nothing) = sized_braced_aircraft(key, edit!)
    wing_chord_at_strut(ac) = ac.wing.layout.root_chord * ac.wing.inboard.λ
    # a strut with an end-fixity factor of 2 buckles before it yields; no strut-attach reinforcement of the fuselage.
    # Its chord is the wing chord: at the default chord ratio no box of this height can reach the buckling stiffness
    function slender!(ac)
        ac.wing.strut.lift_fraction_strut_to_total = 0.0 # these tests are about the box sizing, not the strut lift
        ac.wing.strut.end_fixity_K = 2.0
        ac.wing.strut.chord_to_wing_chord = 1.0
        ac.wing.strut.weight_frac_attachments = 0.0
    end

    # a strut with a low end-fixity factor yields before it buckles; no strut lift (box sizing, not strut lift)
    stocky!(ac) = (ac.wing.strut.lift_fraction_strut_to_total = 0.0; ac.wing.strut.end_fixity_K = 0.1)

    @testset "the strut chord is the chord ratio times the wing chord at the strut attachment" begin
        for ac in (size_braced(), sized_jury_aircraft("truss_braced", 0.3))
            strut = ac.wing.strut
            @test strut.chord_to_wing_chord == 0.5
            @test strut.chord ≈ strut.chord_to_wing_chord * wing_chord_at_strut(ac) rtol = 1e-12
            @test strut.S ≈ 2 * strut.chord * strut.length * strut.cos_lambda rtol = 1e-10
            # without taper the chord at the fuselage end is the chord at the wing end
            @test strut.taper_ratio == 1.0
            @test strut.root_chord == strut.chord
        end
    end

    @testset "a tapered strut is sized with the wing-end chord, and weighs the mean chord over the wing-end chord more" begin
        ac = sized_tapered_strut_aircraft(0.5)
        strut, box = ac.wing.strut, recomputed_strut_box(ac)
        @test strut.area ≈ box.area_caps + box.area_webs rtol = 1e-6
        @test strut.weight ≈ (1 + 1 / 0.5) / 2 * box_member_weight(strut, box) rtol = 1e-6
    end

    @testset "the strut weight is both wings of caps and webs with the added weight fraction, and the moment uses it" begin
        ac = size_braced()
        strut, box = ac.wing.strut, recomputed_strut_box(ac)
        @test strut.weight ≈ box_member_weight(strut, box) rtol = 1e-6
        @test strut.dxW ≈ strut.weight * strut_dx_arm(ac) rtol = 1e-6
        @test strut.area ≈ box.area_caps + box.area_webs rtol = 1e-6 # the area does not contain the added fraction
    end

    @testset "the net of the own weight and the lift of the strut is a transverse load: shear on the webs and bending moment on the caps" begin
        for ac in (size_braced(), sized_jury_aircraft("truss_braced", 0.3))
            strut, jury = ac.wing.strut, ac.wing.jury
            loads = strut_net_transverse_loads(ac)
            @test strut.bending_moment ≈ loads.M rtol = 1e-6
            @test strut.shear_load ≈ abs(jury.axial_force) + loads.V rtol = 1e-6
        end
    end

    @testset "the strut transverse shear and moment are sized for the net of own weight, lift and fuel inertia" begin
        ac = size_braced()
        wing, strut = ac.wing, ac.wing.strut
        results = wing.sizing_loads.results
        n_wings = 2.0

        # Strut fuel inertia load of one wing per case: a tenth of the largest strut lift, with the sign of the load factor
        f_scale = 0.1 * maximum(abs(r.strut_lift) for r in values(results))
        fuel_load(name) = sign(getfield(wing.sizing_loads, name).n_computed) * f_scale

        # Sizing envelope of the recorded cases, as `evaluate_wing_sizing_cases!` returns it
        cases = collect(values(results))
        forces = [c.strut_force for c in cases]
        envelope = (S_root = maximum(max(abs(c.S_root), abs(c.S_break)) for c in cases),
                    M_root = maximum(max(abs(c.M_root), abs(c.M_inboard_extreme)) for c in cases),
                    S_break = maximum(abs(c.S_break) for c in cases), M_break = maximum(abs(c.M_break) for c in cases),
                    T_strut = max(0.0, maximum(forces)), T_strut_compression = min(0.0, minimum(forces)),
                    P_axial = cases[argmax([abs(c.P_axial_inboard) for c in cases])].P_axial_inboard)

        # Transverse shear and midspan moment of the strut after `wing_weights!` with the strut fuel loads `f(name)` in the cases
        function strut_transverse_loads(f)
            w = deepcopy(wing)
            for (name, r) in w.sizing_loads.results
                r.strut_fuel_load = f(name)
            end
            Weng1, (; neinn, dyeinn, neout, dyeout) = ac.parg[igWeng] / ac.parg[igneng], TASOPT.wing_engine_arms(w)
            TASOPT.wing_weights!(w, 1.0e5, 0.3, 0.8, 3.0, Weng1, neout, dyeout, neinn, dyeinn, ac.parg[igsigfac], ac.parg[igrhofuel];
                                 n_wings, has_strut = true, section_loads_envelope = envelope)
            return (V = w.strut.shear_load, M = w.strut.bending_moment)
        end
        # Hand calculation: the net of lift, own weight and fuel acts on the strut, its component perpendicular to the strut over the length
        function expected_transverse_loads(f)
            L, z = strut.length, strut.z
            q_net(name, r) = abs(r.strut_lift - getfield(wing.sizing_loads, name).n_computed * strut.weight / n_wings - f(name)) /
                             L * sqrt(L^2 - z^2) / L
            q_perp = maximum(q_net(name, r) for (name, r) in results)
            return (V = q_perp * L / 2, M = q_perp * L^2 / 8)
        end

        without, with_fuel = strut_transverse_loads(name -> 0.0), strut_transverse_loads(fuel_load)
        @test without.M ≈ expected_transverse_loads(name -> 0.0).M rtol = 1e-9
        @test with_fuel.V ≈ expected_transverse_loads(fuel_load).V rtol = 1e-9
        @test with_fuel.M ≈ expected_transverse_loads(fuel_load).M rtol = 1e-9
        @test with_fuel.V < without.V && with_fuel.M < without.M # the pull-up governs and the fuel inertia opposes its lift
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
    end

    @testset "a stocky strut is governed by tension" begin
        ac = size_braced(:stocky, stocky!)
        strut = ac.wing.strut
        box = recomputed_strut_box(ac)
        @test strut.sizing_criterion === :tension
        # the caps carry at least the tension (the bending moment of the own weight adds to it)
        @test box.area_caps ≥ strut.axial_force / (strut.caps.material.σmax * ac.parg[igsigfac]) * (1 - 1e-6)
        @test strut.buckling_load > strut.caps.material.safety_factor * abs(strut.axial_force_compression)
    end

    @testset "Prescribed strut force is sized by tension only" begin
        ac = size_braced(:prescribed, ac -> ac.wing.strut.opt_strut_force = TASOPT.StrutForceModel.Prescribed)
        strut, sl, layout = ac.wing.strut, ac.wing.sizing_loads, ac.wing.layout
        box = recomputed_strut_box(ac)
        @test Set(keys(sl.results)) == Set([:pullup, :pushover, :pullup_light, :pushover_light, :taxi])
        @test strut.axial_force > 0
        @test strut.axial_force_compression == 0 # the compressive push-over case does not enter the strut
        @test sl.results[:pushover].strut_force < 0
        @test strut.sizing_criterion in (:tension, :gauge) # there is no buckling requirement without compression
        @test strut.area ≈ box.area rtol = 1e-6
        @test box.area_caps * strut.caps.material.σmax * ac.parg[igsigfac] >= strut.axial_force * (1 - 1e-6)

        # without an envelope, wing_weights! takes the strut force of the pull-up, scaled from the length ls to the swept strut length
        ls = sqrt(strut.z^2 + (0.5 * layout.span * (layout.ηs - layout.ηo))^2)
        @test sl.results[:pullup].strut_lift > 0
        @test strut.axial_force * ls / strut.length ≈ sl.results[:pullup].strut_force rtol = 1e-6

        # the strut length is the diagonal over the height of the attachment and the inboard bay along the sweep
        length_3d = sqrt(strut.z^2 + (0.5 * layout.span * (layout.ηs - layout.ηo) / cosd(layout.sweep))^2)
        @test strut.length ≈ length_3d rtol = 1e-6
    end

    @testset "a strut attach aft of the wing box lengthens the strut run and reduces its planform sweep" begin
        for config in ("strut_braced", "truss_braced")
            ac = sized_attach_aft_aircraft(config)
            wing, strut = ac.wing, ac.wing.strut
            layout = wing.layout
            dy = 0.5 * layout.span * (layout.ηs - layout.ηo)
            dx = strut.x_strut_attach - layout.box_x # attach offset, final box_x
            strut_run = dy * tand(layout.sweep) - dx # streamwise run of the strut
            ls = sqrt(strut.z^2 + dy^2)
            @test dx > 0
            @test strut.length ≈ sqrt(ls^2 + strut_run^2) rtol = 1e-4
            @test strut.cos_lambda ≈ ls / strut.length rtol = 1e-4
            @test st.strut_sweep(wing) ≈ atand(strut_run / dy) rtol = 1e-4
            @test st.strut_sweep(wing) < 0.5 * layout.sweep
            @test strut.dxW ≈ strut.weight * (dx + dy * tand(layout.sweep)) / 2 rtol = 1e-4
        end
    end

    @testset "with the strut attach at the wing box the strut sweep is the wing sweep" begin
        ac = size_braced()
        @test isnan(ac.wing.strut.x_strut_attach)
        @test st.strut_sweep(ac.wing) ≈ ac.wing.layout.sweep rtol = 1e-6
    end

    @testset "minimum gauges set the strut box when they are large" begin
        ac = size_braced(:large_gauges, ac -> (ac.wing.strut.lift_fraction_strut_to_total = 0.0; # box sizing, not strut lift
                                               ac.wing.strut.min_cap_thickness = 3.0e-2;
                                               ac.wing.strut.min_web_thickness = 3.0e-2))
        strut = ac.wing.strut
        hw = strut.web_to_box_height * strut.thickness_to_chord
        @test strut.sizing_criterion === :gauge
        # the caps and the webs are at their gauge: area = 2 t w c + 2 t hw c
        @test strut.area ≈ 2 * 3.0e-2 * strut.chord * strut.cos_lambda * (strut.box_width_to_chord + hw) rtol = 1e-6
        @test strut.area ≈ recomputed_strut_box(ac).area rtol = 1e-6
    end
end

# --- Fuel stored in the strut box (sized strut-braced aircraft) ---
@testset "fuel in the strut tank" verbose=true begin
    st = TASOPT.structures
    sized_with_strut_fuel() = sized_braced_aircraft(:strut_fuel, ac -> ac.options.has_strut_fuel = true)

    # Capacity [N] and moment about the wing box [N m] of the wing tanks (both wings) from the recorded cap and web
    # thicknesses, as in `wing_weights!`; `inner_panel_fuel` and `inner_panel_structure` are the weights of one inner panel.
    function wing_tank_capacity_and_moment(ac)
        wing = ac.wing
        c, b, ηo, ηs, λi, λo = wing.layout.root_chord, wing.layout.span, wing.layout.ηo, wing.layout.ηs,
                               wing.inboard.λ, wing.outboard.λ
        sinL, cosL = sind(wing.layout.sweep), cosd(wing.layout.sweep)
        width = wing.inboard.cross_section.width_to_chord
        area(section) = st.calc_sparbox_internal_area(width, first(st.get_average_sparbox_heights(section.cross_section)),
                                                      section.caps.thickness, section.webs.thickness)
        A_root, A_break = area(wing.inboard), area(wing.outboard)
        A_inner = (A_root + A_break * λi^2) / (1 + λi^2)
        V_center = c^2 * b * ηo / 2
        V_inner = c^2 * b * (ηs - ηo) * (1 + λi + λi^2) / 6 * cosL
        V_outer = c^2 * b * (1 - ηs) * (λi^2 + λi * λo + λo^2) / 6 * cosL
        dxV_inner = c^2 * b^2 * (ηs - ηo)^2 * (1 + 2λi + 3λi^2) / 48 * sinL
        dxV_outer = c^2 * b^2 * (1 - ηs)^2 * (λi^2 + 2λi * λo + 3λo^2) / 48 * sinL +
                    c^2 * b^2 * (ηs - ηo) * (1 - ηs) * (λi^2 + λi * λo + λo^2) / 12 * sinL
        weight_per_area = ac.parg[igrhofuel] * TASOPT.gee
        capacity = 2 * weight_per_area * ((ac.options.has_centerbox_fuel ? A_root * V_center : 0.0) +
                                          A_inner * V_inner + A_break * V_outer)
        moment = 2 * weight_per_area * (A_inner * dxV_inner + A_break * dxV_outer)
        web_height = wing.inboard.cross_section.web_to_box_height * wing.inboard.cross_section.thickness_to_chord
        structure_area(section) = section.caps.ρ * 2 * section.caps.thickness * width +
                                  section.webs.ρ * 2 * section.webs.thickness * web_height
        structure_inner = (structure_area(wing.inboard) + structure_area(wing.outboard) * λi^2) / (1 + λi^2)
        return (capacity = capacity, moment = moment, inner_panel_fuel = weight_per_area * A_inner * V_inner,
                inner_panel_structure = TASOPT.gee * structure_inner * V_inner * (1 + TASOPT.wing_additional_weight(wing)))
    end

    @testset "without strut fuel the strut holds none, and the tank capacity and fuel moment are those of the wing" begin
        ac = sized_braced_aircraft()
        wing_tank = wing_tank_capacity_and_moment(ac)
        parg = ac.parg
        @test !ac.options.has_strut_fuel
        @test ac.wing.strut.fuel_capacity == 0
        @test parg[igWfmax] ≈ wing_tank.capacity rtol = 1e-6
        @test parg[igdxWfuel] ≈ wing_tank.moment * parg[igWfuel] / parg[igWfmax] rtol = 1e-6
    end

    @testset "the strut tank capacity is the internal volume of the strut box times the fuel weight density, for both wings" begin
        ac = sized_with_strut_fuel()
        strut, box = ac.wing.strut, recomputed_strut_box(ac)
        area = st.calc_sparbox_internal_area(strut.box_width_to_chord, strut.thickness_to_chord, box.tbcap, box.tbweb)
        mean_chord_squared = (strut.chord^2 + strut.chord * strut.root_chord + strut.root_chord^2) / 3 * strut.cos_lambda^2 # linear taper, perpendicular to the strut axis
        @test ac.options.has_strut_fuel
        @test strut.fuel_capacity > 0
        @test strut.fuel_capacity ≈ 2 * ac.parg[igrhofuel] * TASOPT.gee * area * mean_chord_squared * strut.length rtol = 1e-6
    end

    @testset "the wing tanks and the strut tanks form one tank capacity, filled to one fraction" begin
        ac = sized_with_strut_fuel()
        wing_tank = wing_tank_capacity_and_moment(ac)
        strut, parg = ac.wing.strut, ac.parg
        Wfmax = wing_tank.capacity + strut.fuel_capacity
        @test parg[igWfmax] ≈ Wfmax rtol = 1e-6

        # the fuel moment is that of the tanks at the one fill fraction; the strut tank acts at the strut centroid
        strut_moment = strut.fuel_capacity * strut_dx_arm(ac)
        @test parg[igdxWfuel] ≈ (wing_tank.moment + strut_moment) * parg[igWfuel] / Wfmax rtol = 1e-6

        # the wing panels carry the fuel at that fraction of the combined capacity (not of the wing capacity alone)
        fill = (parg[igWpay] - parg[igWpaymax] + parg[igWfuel]) / Wfmax
        @test ac.wing.inboard.weight ≈ wing_tank.inner_panel_structure + fill * wing_tank.inner_panel_fuel rtol = 1e-3
    end
end

# --- Jury strut sizing ---
@testset "jury strut sizing" verbose=true begin
    st = TASOPT.structures
    truss = sized_jury_aircraft("truss_braced", 0.3)

    @testset "strut_buckling_length is the longer part of the strut cut at the jury" begin
        @test st.strut_buckling_length(15.0, 0.5) ≈ 7.5 rtol = 1e-12
        @test st.strut_buckling_length(15.0, 0.3) ≈ 10.5 rtol = 1e-12
        @test st.strut_buckling_length(15.0, nothing) == 15.0
    end

    @testset "the jury follows the strut and is sized as a pinned box under a fraction of the strut compression" begin
        strut, jury, layout = truss.wing.strut, truss.wing.jury, truss.wing.layout
        box = recomputed_jury_box(truss)
        @test jury.length ≈ strut.z * (1 - 0.3) rtol = 1e-8
        @test jury.eta_wing ≈ layout.ηo + 0.3 * (layout.ηs - layout.ηo) rtol = 1e-8
        @test jury.axial_force ≈ -0.02 * abs(strut.axial_force_compression) rtol = 1e-8 # 2 % of the strut compression
        @test jury.area ≈ box.area rtol = 1e-6
        @test jury.buckling_load ≈ box.buckling_load rtol = 1e-6
        @test jury.buckling_load >= jury.caps.material.safety_factor * abs(jury.axial_force) * (1 - 1e-6)
    end

    @testset "the jury chord is the chord ratio times the strut chord at the jury attachment" begin
        for ac in (truss, sized_tapered_strut_aircraft(0.5))
            strut, jury = ac.wing.strut, ac.wing.jury
            local_strut_chord = strut.root_chord + (strut.chord - strut.root_chord) * jury.eta_strut
            @test jury.chord ≈ jury.chord_to_strut_chord * local_strut_chord rtol = 1e-12
            @test jury.area ≈ recomputed_jury_box(ac).area rtol = 1e-6
        end
    end

    @testset "jury weight is the weight of the caps and webs of both wings times (1 + weight_frac_added), and the moment uses it" begin
        jury, box = truss.wing.jury, recomputed_jury_box(truss)
        @test jury.weight ≈ box_member_weight(jury, box) rtol = 1e-6
        @test jury.dxW ≈ jury.weight * jury_dx_arm(truss) rtol = 1e-6
    end

    @testset "with the strut attach aft of the wing box the jury spans the attach offset as well as the height" begin
        ac = sized_attach_aft_aircraft("truss_braced")
        wing, strut, jury = ac.wing, ac.wing.strut, ac.wing.jury
        layout = wing.layout
        dy = 0.5 * layout.span * (layout.ηs - layout.ηo)
        dx = strut.x_strut_attach - layout.box_x
        @test jury.length ≈ (1 - jury.eta_strut) * sqrt(strut.z^2 + dx^2) rtol = 1e-4
        @test jury.dxW ≈ jury.weight * (jury.eta_strut * dy * tand(layout.sweep) + 0.5 * dx * (1 - jury.eta_strut)) rtol = 1e-4
    end
end

# --- Wing sizing load cases (mass states, root loading, envelope, braced section loads) ---
@testset "wing sizing load cases" verbose=true begin
    # The default wing; its span and section weights are only set when the aircraft is sized
    function default_wing_with_span()
        wing = CoreLogging.with_logger(CoreLogging.NullLogger()) do
            load_default_model().wing
        end
        wing.layout.span = 60.0
        return wing
    end
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

        # expected (W_lift, r_fuel_wing, W_fuel_fuse) with fuel in the wing, and with fuel in the fuselage
        wing_fuel = (
            (WingMassState.DesignCruise, 900e3, 0.6, 0.0),
            (WingMassState.MTOW_MaxFuel, 1000e3, 0.5, 0.0),
            (WingMassState.MZFW, 830e3, 0.1, 0.0),
            (WingMassState.Light, 530e3, 0.1, 0.0),
        )
        fuselage_fuel = (
            (WingMassState.MTOW_MaxFuel, 1000e3, 0.0, 150e3),
            (WingMassState.Light, 530e3, 0.0, 30e3),
        )

        for (has_wing_fuel, expected) in ((true, wing_fuel), (false, fuselage_fuel)),
            (state, W_lift, r_fuel_wing, W_fuel_fuse) in expected
            result = wing_mass_state(state, has_wing_fuel; weights...)
            @test result.W_lift ≈ W_lift
            @test result.r_fuel_wing ≈ r_fuel_wing
            @test result.W_fuel_fuse ≈ W_fuel_fuse
        end
    end

    @testset "masses of the cases" begin
        # 1.7e5 = W_pay - W_paymax + W_fuel is the fuel that fits at MTOW with maximum payload
        m = TASOPT.wing_sizing_masses(; W_cruise = 9.0e5, W_MTO = 1.0e6, W_pay = 1.5e5, W_paymax = 2.0e5, W_fuel = 2.2e5,
                                      freserve = 0.1, W_fuel_cap = 3.0e5, has_wing_fuel = true)
        @test (m.W_cruise, m.W_MTO, m.W_fuel, m.W_fuel_cap) == (9.0e5, 1.0e6, 2.2e5, 3.0e5)
        @test m.W_reserve ≈ 2.0e4
        @test m.W_empty ≈ 6.3e5
        @test m.W_fuel_at_MTO ≈ 1.7e5
        @test m.r_fuel_design ≈ 1.7e5 / 3.0e5
    end

    @testset "load factor of the cases" verbose=true begin
        wing_case_load_factor = TASOPT.wing_case_load_factor
        far_manoeuvre_load_factor = TASOPT.structures.far_manoeuvre_load_factor
        cases = default_wing_with_span().sizing_loads
        W_MTO = 5.5e5

        @testset "with FAR only the positive manoeuvre cases follow the FAR load factor of W_MTO" begin
            n_far = far_manoeuvre_load_factor(W_MTO)
            @test n_far != cases.pullup.load_factor
            for name in (:pullup, :pullup_light)
                @test wing_case_load_factor(name, getfield(cases, name), W_MTO, true) == n_far
            end
            for name in (:pushover, :pushover_light, :taxi)
                @test wing_case_load_factor(name, getfield(cases, name), W_MTO, true) == getfield(cases, name).load_factor
            end
        end

        @testset "gust cases" verbose=true begin
            wing = default_wing_with_span()
            wing.layout.S = wing.layout.span^2 / wing.layout.AR
            W_lift, Mach_cruise, V_cruise_EAS = 3.2e5, 0.8, 154.3   # N, -, m/s EAS (300 kt)
            gust_up, gust_down = deepcopy(cases.gust_up_light), deepcopy(cases.gust_down_light)
            gust_factor(name, case; kwargs...) = wing_case_load_factor(name, case, W_MTO, false;
                wing, W_lift, Mach_cruise, V_cruise_EAS, kwargs...)
            function expected_increment(case, V)
                S, span = wing.layout.S, wing.layout.span
                CL_alpha = TASOPT.wing_lift_curve_slope(wing.layout.AR, wing.layout.sweep, Mach_cruise)
                return TASOPT.gust_load_increment(W_lift / S, S / span, CL_alpha, case.gust_velocity_EAS, V)
            end

            @testset "the load factor is 1 ± Δn with the sign of the load factor of the case" begin
                Δn = expected_increment(gust_up, V_cruise_EAS)
                @test Δn > 0
                @test gust_factor(:gust_up_light, gust_up) ≈ 1 + Δn rtol = 1e-12
                @test gust_factor(:gust_down_light, gust_down) ≈ 1 - Δn rtol = 1e-12
            end

            @testset "the speed of the case overrides the cruise speed, which is used if the case has no speed" begin
                case = deepcopy(gust_up)
                case.gust_speed_EAS = 100.0
                @test gust_factor(:gust_up_light, case; V_cruise_EAS = 200.0) ≈ 1 + expected_increment(case, 100.0) rtol = 1e-12
                case.gust_speed_EAS = 0.0
                @test gust_factor(:gust_up_light, case; V_cruise_EAS = 100.0) ≈ 1 + expected_increment(case, 100.0) rtol = 1e-12
            end

            @testset "a gust case needs a positive speed, a subsonic Mach number and the gust context" begin
                @test_throws ArgumentError gust_factor(:gust_up_light, gust_up; V_cruise_EAS = 0.0)
                @test_throws ArgumentError gust_factor(:gust_up_light, gust_up; Mach_cruise = 1.0)
                @test_throws Exception wing_case_load_factor(:gust_up_light, gust_up, W_MTO, false) # no gust context
            end
        end
    end

    @testset "root loading of the cases" begin
        wing_loading = TASOPT.wing_loading
        wing = default_wing_with_span()
        cases = wing.sizing_loads
        W_lift, n_pullup = 5.0e5, 3.0
        root_loading(case) = TASOPT.wing_case_root_loading(case, wing, 1.0, 1.0, W_lift, Lhtail, n_pullup)

        loads = root_loading(cases.pullup) # wing_loading at the pull-up load factor
        @test loads.Nload == 3.0
        @test loads.po == wing_loading(wing, 1.0, 1.0, n_pullup, W_lift, Lhtail)

        loads = root_loading(cases.pushover) # the tail load scales with the load factor and with tail_load_scale
        @test loads.Nload == -1.0
        @test loads.po == wing_loading(wing, 1.0, 1.0, -1.0, W_lift, (-1.0 / n_pullup) * Lhtail)
        half = deepcopy(cases.pushover)
        half.tail_load_scale = 0.5
        @test root_loading(half).po == wing_loading(wing, 1.0, 1.0, -1.0, W_lift, 0.5 * (-1.0 / n_pullup) * Lhtail)

        loads = root_loading(cases.taxi) # ground: no wing lift
        @test (loads.Nload, loads.po) == (-2.0, 0.0)
    end

    @testset "envelope of the cases" verbose=true begin
        Cantilever = TASOPT.structures.Cantilever
        wing_case_loads = TASOPT.wing_case_loads
        wing_case_root_loading = TASOPT.wing_case_root_loading

        base_wing = default_wing_with_section_weights()
        base_wing.inboard.max_shear_load, base_wing.inboard.moment = 111.0, 222.0 # sentinels, left untouched by the evaluation
        base_wing.outboard.max_shear_load, base_wing.outboard.moment = 333.0, 444.0

        rclt, rcls = 1.0, 1.0
        W_lift_cruise, W_lift_MTO, W_lift_light = 5.0e5, 5.5e5, 3.2e5 # lift weights of the mass states, from `masses`

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

            M_root = [abs(sl.results[name].M_root) for name in active]
            M_break = [abs(sl.results[name].M_break) for name in active]
            @test sl.governing_root == active[argmax(M_root)]
            @test sl.governing_break == active[argmax(M_break)]
            governing_root, governing_break = sl.results[sl.governing_root], sl.results[sl.governing_break]
            @test envelope == (S_root = abs(governing_root.S_root), M_root = abs(governing_root.M_root),
                               S_break = abs(governing_break.S_break), M_break = abs(governing_break.M_break))

            @test evaluate(wing) == envelope # evaluating again gives the same results
        end

        @testset "stale results are removed, and no active case gives no envelope" begin
            wing = deepcopy(base_wing)
            sl = wing.sizing_loads
            sl.results[:taxi] = TASOPT.structures.WingCaseLoads(M_root = 1.0e12)
            envelope = evaluate(wing)
            @test Set(keys(sl.results)) == Set([:pullup])
            @test envelope.M_root == abs(sl.results[:pullup].M_root)

            sl.pullup.active = false
            @test evaluate(wing) === nothing
            @test isempty(sl.results)
        end

        @testset "gust cases enter the results with opposite loads" begin
            wing = deepcopy(base_wing)
            wing.layout.S = wing.layout.span^2 / wing.layout.AR # the default wing has no area until it is sized
            sl = wing.sizing_loads
            sl.gust_up_light.active = sl.gust_down_light.active = true
            TASOPT.evaluate_wing_sizing_cases!(wing, masses, true, rclt, rcls, Lhtail, gammat, gammas, We, dyeinn;
                                               Mach_cruise = 0.8, V_cruise_EAS = 154.3)
            @test sl.gust_up_light.n_computed > 1 > sl.gust_down_light.n_computed
            @test Set(keys(sl.results)) == Set([:pullup, :gust_up_light, :gust_down_light])
            @test sl.results[:gust_up_light].M_root * sl.results[:gust_down_light].M_root < 0 # opposite directions
        end
    end

    @testset "braced wing section loads" verbose=true begin
        st = TASOPT.structures
        Cantilever = st.Cantilever
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
        case_names = (:pullup, :pushover, :pullup_light, :pushover_light, :taxi)

        # A wing without inboard weight and engines, loaded by the lift only
        lift_only(model) = wing_with(model; edit! = w -> (w.inboard.weight = 0.0; w.inboard.dyW = 0.0))

        @testset "evaluate_wing_sizing_cases! records the strut loads of a braced wing" begin
            section_loads(w) = (w.inboard.max_shear_load, w.inboard.moment, w.outboard.max_shear_load, w.outboard.moment)

            for model in (StrutForceModel.Compatibility, StrutForceModel.Prescribed)
                wing = wing_with(model)
                sl = wing.sizing_loads
                for name in (:pushover, :pullup_light, :pushover_light, :taxi)
                    getfield(sl, name).active = true
                end
                sections_before = section_loads(wing)
                envelope = TASOPT.evaluate_wing_sizing_cases!(wing, masses, true, 1.0, 1.0, Lhtail,
                                                              gammat, gammas, We, dyeinn; has_strut = true)
                results = sl.results

                @test Set(keys(results)) == Set(case_names)
                @test section_loads(wing) == sections_before
                forces = Dict(name => results[name].strut_force for name in case_names)
                # a pull-up loads the strut in tension, a push-over in compression
                @test forces[:pullup] > 0 && forces[:pushover] < 0
                @test results[:pullup].P_axial_inboard ≈ -forces[:pullup] * L / ls rtol = 1e-12

                # the largest tension and the most negative strut force govern the strut
                @test sl.governing_strut == case_names[argmax([forces[name] for name in case_names])]
                @test envelope.T_strut ≈ maximum(values(forces)) rtol = 1e-12
                @test sl.governing_strut_compression == case_names[argmin([forces[name] for name in case_names])]
                @test envelope.T_strut_compression ≈ minimum(values(forces)) rtol = 1e-12

                # the inboard section is sized for the largest moment over the cases
                @test envelope.M_root ≈ maximum(max(abs(results[name].M_root), abs(results[name].M_inboard_extreme))
                                                for name in case_names)
            end
        end

        # Strut lift: the struts carry the share k = lift_fraction_strut_to_total of the lift;
        # the wing carries (1 - k), and the strut lift of one wing is k (N W - Lh) / 2.
        n_pullup = 3.0
        case_load_factor = (pullup = 3.0, pushover = -1.0)

        "Set the strut geometry of `wing` and the `lift_fraction_strut_to_total` `k`, the share of the lift the struts carry."
        function with_strut_lift!(wing, k; p = 0.5)
            strut = wing.strut
            wing.layout.S = wing.layout.span^2 / wing.layout.AR
            strut.S, strut.dihedral, strut.local_velocity_ratio = 1.6 * ls, atan(z, L), 1.05
            strut.spanload_exponent = p
            strut.lift_fraction_strut_to_total = k
            return wing
        end

        same_loads(a, b) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))
        evaluate_braced(wing; kwargs...) = TASOPT.evaluate_wing_sizing_cases!(wing, masses, true, 1.0, 1.0, Lhtail,
            gammat, gammas, We, dyeinn; has_strut = true, kwargs...)

        "Closed-form loads of case `name` of a `lift_only` wing for the strut lift share `k`, with the fraction `f_wing` of the strut lift passed to the wing end."
        function expected_lift_case(wing, name, k, f_wing)
            N = case_load_factor[name]
            Lh = getfield(wing.sizing_loads, name).tail_load_scale * (N / n_pullup) * Lhtail
            po_unreduced = TASOPT.wing_loading(wing, 1.0, 1.0, N, W_lift[name], Lh)
            po = (1 - k) * po_unreduced
            F_strut = k * (N * W_lift[name] - Lh) / 2
            F_wing = f_wing * F_strut
            outboard = wing_case_loads(Cantilever(), deepcopy(wing), po, gammat, gammas, N, We, 0, 0.0, 0, dyeinn)
            function M0(x)
                u = L - x
                a, c = po * (1 - (1 - gammas) * x / L), po * (1 - gammas) / L
                return outboard.M_break + outboard.S_break * u + a * u^2 / 2 - c * u^3 / 3 + F_wing * u
            end
            if wing.strut.opt_strut_force == StrutForceModel.Compatibility
                EI(x) = wing.inboard.EI[1] + (wing.outboard.EI[1] - wing.inboard.EI[1]) * x / L
                f_strut = ls / (wing.strut.caps.material.E * wing.strut.area * (z / ls)^2)
                R = st.strut_force_compatibility(M0, EI, L, f_strut)
            else
                # the prescribed rule M_root = M_break, which the point load F_wing at the attachment adds to
                R = (po * span / 12) * (ηs - ηo) * (1 + 2gammas) + outboard.S_break + F_wing
            end
            return (F_strut = F_strut, strut_force = R * ls / z + (1 - f_wing) * F_strut * z / ls,
                    S_break = outboard.S_break, M_break = outboard.M_break,
                    S_root = outboard.S_break + po * L * (1 + gammas) / 2 + F_wing - R, M_root = M0(0) - R * L)
        end

        @testset "a strut that carries no lift leaves every load unchanged" begin
            for model in (StrutForceModel.Compatibility, StrutForceModel.Prescribed)
                wing = with_strut_lift!(wing_with(model), 0.0; p = 0.5)
                reference = wing_with(model)
                for w in (wing, reference), name in (:pushover, :pullup_light, :pushover_light, :taxi)
                    getfield(w.sizing_loads, name).active = true
                end
                envelope, reference_envelope = evaluate_braced(wing), evaluate_braced(reference)

                @test envelope == reference_envelope
                for name in case_names
                    r = wing.sizing_loads.results[name]
                    @test r.strut_lift == 0.0
                    @test same_loads(r, reference.sizing_loads.results[name])
                end
            end
        end

        @testset "closed form of the strut lift loads, with the prescribed and the compatibility strut force" begin
            k = 0.2
            for model in (StrutForceModel.Compatibility, StrutForceModel.Prescribed)
                lifting, without = with_strut_lift!(lift_only(model), k), with_strut_lift!(lift_only(model), 0.0)
                lifting.sizing_loads.pushover.active = without.sizing_loads.pushover.active = true
                envelope = evaluate_braced(lifting; neinn = 0)
                evaluate_braced(without; neinn = 0)

                for name in (:pullup, :pushover)
                    r, r_without = lifting.sizing_loads.results[name], without.sizing_loads.results[name]
                    expected = expected_lift_case(lifting, name, k, 0.5)   # a uniform strut loading puts half on each end
                    @test r.strut_lift ≈ expected.F_strut rtol = 1e-12
                    @test r.S_break ≈ expected.S_break rtol = 1e-9
                    @test r.M_break ≈ expected.M_break rtol = 1e-9
                    @test r.strut_force ≈ expected.strut_force rtol = 1e-9
                    @test r.S_root ≈ expected.S_root rtol = 1e-9
                    @test r.M_root ≈ expected.M_root rtol = 1e-9
                    @test sign(r.strut_lift) == sign(case_load_factor[name]) # the strut lift has the sign of the load factor
                    @test r_without.strut_lift == 0.0
                end

                # the prescribed force makes the root moment equal the moment at the strut attachment
                model == StrutForceModel.Prescribed &&
                    @test without.sizing_loads.results[:pullup].M_root ≈ without.sizing_loads.results[:pullup].M_break rtol = 1e-9

                # the strut lift lowers the strut force of the pull-up, which is the largest tension
                r, r_without = lifting.sizing_loads.results[:pullup], without.sizing_loads.results[:pullup]
                @test 0 < r.strut_force < r_without.strut_force
                @test envelope.T_strut ≈ r.strut_force rtol = 1e-12
            end
        end

        @testset "a cantilever wing ignores the strut lift inputs" begin
            lifting, plain = with_strut_lift!(deepcopy(base_wing), 0.3), deepcopy(base_wing)
            for w in (lifting, plain), name in (:pushover, :taxi)
                getfield(w.sizing_loads, name).active = true
            end
            evaluate_cantilever(w) = TASOPT.evaluate_wing_sizing_cases!(w, masses, true, 1.0, 1.0, Lhtail, gammat, gammas, We, dyeinn)

            @test evaluate_cantilever(lifting) == evaluate_cantilever(plain)
            for name in (:pullup, :pushover, :taxi)
                @test same_loads(lifting.sizing_loads.results[name], plain.sizing_loads.results[name])
            end
        end

        # Strut fuel: the inertia load of the fuel in one strut, positive downward, adds to the strut lift load case by case.
        # Half of it goes to the wing end of the strut and half to the fuselage end.
        fuel_load_wing(model) = (w = lift_only(model); w.strut.spanload_exponent = 0.5; w)
        po_fuel, Nload_fuel, F_strut_fuel, f_fuel = 6.0e4, 3.0, 5.0e4, 2.0e4

        "Closed-form loads for the strut lift `F` and the strut fuel load `f`: the wing end of the strut carries `ξcp F - f/2` as a point load at the attachment."
        function expected_fuel_case(wing, F, f)
            ξcp = TASOPT.strut_lift_centroid_fraction(wing.strut.spanload_exponent)
            F_end = ξcp * F - f / 2
            outboard = wing_case_loads(Cantilever(), deepcopy(wing), po_fuel, gammat, gammas, Nload_fuel, We, 0, 0.0, 0, dyeinn)
            function M0(x)
                u = L - x
                a, c = po_fuel * (1 - (1 - gammas) * x / L), po_fuel * (1 - gammas) / L
                return outboard.M_break + outboard.S_break * u + a * u^2 / 2 - c * u^3 / 3 + F_end * u
            end
            if wing.strut.opt_strut_force == StrutForceModel.Compatibility
                EI(x) = wing.inboard.EI[1] + (wing.outboard.EI[1] - wing.inboard.EI[1]) * x / L
                f_strut = ls / (wing.strut.caps.material.E * wing.strut.area * (z / ls)^2)
                R = st.strut_force_compatibility(M0, EI, L, f_strut)
            else
                R = (po_fuel * span / 12) * (ηs - ηo) * (1 + 2gammas) + outboard.S_break + F_end
            end
            return (R = R, S_root = outboard.S_break + po_fuel * L * (1 + gammas) / 2 + F_end - R, M_root = M0(0) - R * L,
                    strut_force = R * ls / z + ((F - f) - F_end) * z / ls)
        end

        @testset "wing_case_loads of a braced wing: the strut fuel load adds to the wing-end point load and the fuselage-end force" begin
            braced_loads(wing; kwargs...) = wing_case_loads(st.Braced(), deepcopy(wing), po_fuel, gammat, gammas, Nload_fuel,
                                                            We, 0, 0.0, 0, dyeinn; kwargs...)
            for model in (StrutForceModel.Compatibility, StrutForceModel.Prescribed)
                wing = fuel_load_wing(model)
                without = braced_loads(wing; strut_lift = F_strut_fuel)
                with_fuel = braced_loads(wing; strut_lift = F_strut_fuel, strut_fuel_load = f_fuel)
                expected = expected_fuel_case(wing, F_strut_fuel, f_fuel)

                @test same_loads(braced_loads(wing; strut_lift = F_strut_fuel, strut_fuel_load = 0.0), without)
                @test with_fuel.strut_fuel_load == f_fuel
                @test without.strut_fuel_load == 0.0
                @test with_fuel.strut_lift == F_strut_fuel # the fuel load is not strut lift
                # the point load sits at the attachment: the outboard section does not see it
                @test (with_fuel.S_break, with_fuel.M_break) == (without.S_break, without.M_break)
                @test with_fuel.strut_reaction ≈ expected.R rtol = 1e-9
                @test with_fuel.S_root ≈ expected.S_root rtol = 1e-9
                @test with_fuel.M_root ≈ expected.M_root rtol = 1e-9
                @test with_fuel.strut_force ≈ expected.strut_force rtol = 1e-9
                @test with_fuel.strut_force < without.strut_force # the fuel inertia relieves the strut lift
                if model == StrutForceModel.Prescribed
                    # the prescribed reaction takes the whole point load, so the root loads do not change
                    @test with_fuel.S_root ≈ without.S_root rtol = 1e-9
                    @test with_fuel.M_root ≈ without.M_root rtol = 1e-9
                else
                    @test with_fuel.M_root != without.M_root
                end
            end
        end

        @testset "evaluate_wing_sizing_cases! records the strut fuel inertia load n r_fuel_wing capacity / 2 of each case" begin
            fuel_capacity = 4.0e4 # fuel that fills the struts of both wings [N]
            r_fuel_wing = (pullup = 0.5, pushover = 8.0e4 / 1.5e5, pullup_light = 2.0e4 / 1.5e5, taxi = 8.0e4 / 1.5e5)
            n_case = (pullup = 3.0, pushover = -1.0, pullup_light = 3.0, taxi = -2.0)
            function wing_with_fuel(capacity)
                wing = with_strut_lift!(fuel_load_wing(StrutForceModel.Compatibility), 0.2; p = 0.5)
                wing.strut.fuel_capacity = capacity
                for name in (:pushover, :pullup_light, :taxi)
                    getfield(wing.sizing_loads, name).active = true
                end
                return wing
            end
            wing, dry = wing_with_fuel(fuel_capacity), wing_with_fuel(0.0)
            evaluate_braced(wing; neinn = 0)
            evaluate_braced(dry; neinn = 0)

            for name in keys(n_case)
                r, r_dry = wing.sizing_loads.results[name], dry.sizing_loads.results[name]
                expected_load = n_case[name] * r_fuel_wing[name] * fuel_capacity / 2
                @test r.strut_fuel_load ≈ expected_load rtol = 1e-12
                @test sign(r.strut_fuel_load) == sign(n_case[name]) # the inertia load is upward for n < 0
                @test r_dry.strut_fuel_load == 0.0
                @test r.strut_lift == r_dry.strut_lift

                # the case equals the direct evaluation with the fuel load of the case
                (Nload, po) = wing_case_root_loading(getfield(wing.sizing_loads, name), wing, 1.0, 1.0, W_lift[name],
                                                     Lhtail, n_pullup; strut_share = TASOPT.strut_lift_share(wing))
                direct = wing_case_loads(st.Braced(), deepcopy(wing), po, gammat, gammas, Nload, We, 0, 0.0, 0, dyeinn;
                                         strut_lift = r.strut_lift, strut_fuel_load = expected_load)
                @test same_loads(r, direct)
                # the fuel inertia opposes the lift of the strut: it relieves the strut force of the pull-up (n > 0)
                n_case[name] > 0 && r.strut_lift != 0 && @test r.strut_force < r_dry.strut_force
            end
        end
    end
end

# --- Weight build-up and CG: jury strut and strut-attach reinforcement ---
@testset "jury strut weight in the aircraft weight build-up" verbose=true begin
    truss = sized_jury_aircraft("truss_braced", 0.3)

    # Copy of `ac` with `extra_weight` [N] added to the jury and its moment set for the same arm
    function with_extra_jury_weight(ac, extra_weight)
        heavier = deepcopy(ac)
        heavier.wing.jury.weight += extra_weight
        heavier.wing.jury.dxW = heavier.wing.jury.weight * jury_dx_arm(heavier)
        return heavier
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

# Fuselage reinforcement at the strut attachment: a fraction of the strut weight of the previous sizing iteration, lumped at the wing-box x.
@testset "fuselage strut-attach reinforcement" verbose=true begin
    ac = sized_braced_aircraft()
    reinforcement, strut = ac.fuselage.strut_attach, ac.wing.strut
    @test reinforcement.W > 0
    @test reinforcement.W ≈ strut.weight_frac_attachments * strut.weight rtol = 1e-3
    @test ac.fuselage.strut_attach_criterion === :strut_weight
    @test reinforcement.x ≈ ac.wing.layout.box_x rtol = 1e-3
end

"Physical cap and web thickness [m] of a wing section, from the thickness normalised by the perpendicular chord."
function wing_section_gauges(wing, section)
    perpendicular_chord = section.co * cosd(wing.layout.sweep)
    return (cap = section.caps.thickness * perpendicular_chord,
            web = section.webs.thickness * perpendicular_chord)
end

@testset "sized braced wings respect the minimum gauge" verbose=true begin
    for ac in (sized_jury_aircraft("strut_braced"), sized_jury_aircraft("truss_braced", 0.3)),
        section in (ac.wing.inboard, ac.wing.outboard)
        @test section.min_cap_thickness > 0
        @test section.min_web_thickness > 0
        gauges = wing_section_gauges(ac.wing, section)
        @test gauges.cap ≥ section.min_cap_thickness * (1 - 1e-9)
        @test gauges.web ≥ section.min_web_thickness * (1 - 1e-9)
    end
end

# --- Weight build-up output: the empty-weight subtotal closes with WMTO - Wfuel - Wpay ---
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
            ("truss_braced", sized_jury_aircraft("truss_braced", 0.3),  true,  true))
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

# --- Inboard shear and moment distribution of a sizing load case, with the jump of the strut reaction ---
@testset "inboard load distribution of the sizing load cases" verbose=true begin
    @testset "cantilever wing stores no inboard distribution" begin
        sl = sized_default_aircraft().wing.sizing_loads
        @test !isempty(sl.results)
        for r in values(sl.results)
            @test isempty(r.inboard_stations)
            @test r.strut_reaction == 0.0
        end
    end

    for (name, ac) in (("strut_braced", sized_jury_aircraft("strut_braced")),)
        @testset "$name governing strut case" begin
            wing = ac.wing
            sl = wing.sizing_loads
            @test sl.governing_strut != :none
            r = sl.results[sl.governing_strut]

            L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
            z = wing.strut.z
            ls = sqrt(z^2 + L^2)

            @test length(r.inboard_stations) == length(r.inboard_shear) == length(r.inboard_moment) == 11
            @test r.inboard_stations[1] == 0.0
            @test r.inboard_stations[end] ≈ L rtol = 1e-9

            # the distribution starts at the sized root loads
            @test r.inboard_shear[1] ≈ r.S_root rtol = 1e-6
            @test r.inboard_moment[1] ≈ r.M_root rtol = 1e-6

            # the moment is continuous at the strut attachment, the shear jumps there by the strut reaction
            @test r.inboard_moment[end] ≈ r.M_break rtol = 1e-6
            # the strut force is the one at the fuselage end: R ls / z + (1 - ξ_cp) F sin(dihedral)
            ξ_cp = TASOPT.strut_lift_centroid_fraction(wing.strut.spanload_exponent)
            @test r.strut_reaction ≈ (r.strut_force - (1 - ξ_cp) * r.strut_lift * z / ls) * z / ls rtol = 1e-9
            @test r.strut_reaction != 0.0
            @test r.inboard_shear[end] + r.strut_reaction ≈ r.S_break rtol = 1e-6

            # the stored extreme moment is the moment of largest magnitude of the distribution
            @test r.M_inboard_extreme == r.inboard_moment[argmax(abs.(r.inboard_moment))]
        end
    end
end
