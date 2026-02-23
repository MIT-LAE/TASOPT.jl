"""
Enums for aircraft-level configuration options.
"""

using EnumX

"""
    EngineLocation

Engine mounting position on the airframe.

- `Wing`: engines mounted under or on the wing
- `Fuselage`: engines mounted on the aft fuselage
"""
@enumx EngineLocation Wing Fuselage

"""
    PropSysArch

Propulsion system architecture.

- `TF`: conventional turbofan
- `TE`: turboelectric
- `ConstantTSFC`: simplified constant-TSFC model
- `FuelCellWithDuctedFan`: PEM fuel cell driving a ducted fan
"""
@enumx PropSysArch TF TE ConstantTSFC FuelCellWithDuctedFan

"""
    WingMove

Strategy for positioning the wing during horizontal-tail sizing.

- `Fixed`: wing position fixed
- `FixedCLh`: move wing to achieve `CLh = CLhspec` in cruise
- `MinStaticMargin`: move wing to achieve minimum static margin = `SMmin`
"""
@enumx WingMove Fixed FixedCLh MinStaticMargin

"""
    FuelType

Fuel type carried by the aircraft.

- `JetA`: conventional jet fuel (Jet-A)
- `LH2`: liquid hydrogen
- `CH4`: liquid methane
"""
@enumx FuelType JetA LH2 CH4

# ---------------------------------------------------------------------------
# Canonical string representations (for TOML via save_model)
# These match the strings expected by read_input.jl
# ---------------------------------------------------------------------------

function Base.string(loc::EngineLocation.T)
    loc == EngineLocation.Wing     && return "wing"
    loc == EngineLocation.Fuselage && return "fuselage"
    error("Unknown EngineLocation value: $loc")
end

function Base.string(arch::PropSysArch.T)
    arch == PropSysArch.TF                    && return "tf"
    arch == PropSysArch.TE                    && return "te"
    arch == PropSysArch.ConstantTSFC          && return "constant_tsfc"
    arch == PropSysArch.FuelCellWithDuctedFan && return "fuel_cell_with_ducted_fan"
    error("Unknown PropSysArch value: $arch")
end

function Base.string(move::WingMove.T)
    move == WingMove.Fixed           && return "fixed"
    move == WingMove.FixedCLh        && return "fixed_CLh"
    move == WingMove.MinStaticMargin && return "min_static_margin"
    error("Unknown WingMove value: $move")
end

function Base.string(fuel::FuelType.T)
    fuel == FuelType.JetA && return "JET-A"
    fuel == FuelType.LH2  && return "LH2"
    fuel == FuelType.CH4  && return "CH4"
    error("Unknown FuelType value: $fuel")
end
