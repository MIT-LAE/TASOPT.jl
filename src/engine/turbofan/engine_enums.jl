"""
Enums for turbofan engine calculation options.
"""

using EnumX

"""
    CalcMode

Selects which turbofan calculation to perform.

- `Sizing`: on-design sizing (`tfsize!`)
- `FixedTt4OffDes`: off-design at fixed turbine inlet temperature
- `FixedFeOffDes`: off-design at fixed net thrust
"""
@enumx CalcMode Sizing FixedTt4OffDes FixedFeOffDes

"""
    CoolingOpt

Selects the turbine cooling model.

- `NoCooling`: no cooling air, station 41 == station 4
- `FixedCoolingFlowRatio`: cooling bypass ratios `epsrow` are inputs; `Tmrow` is computed
- `FixedTmetal`: metal temperatures `Tmrow` are inputs; `epsrow` is computed
"""
@enumx CoolingOpt NoCooling FixedCoolingFlowRatio FixedTmetal
