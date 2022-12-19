# Propulsion system

Currently TASOPT uses Numerical Propulsion System Simulation (NPSS) to perform gas turbine calculations.

The following functions are wrappers to access NPSS.

```@docs
startNPSS
endNPSS
NPSS_run
```

The following functions are specifc to the turbo-electric system programmed into NPSS

```@docs
NPSS_TEsys
NPSS_TEsysOD
```