module aerodynamics

export airfun, airtable, trefftz1, cdsum!, cfturb, cditrp, surfcd, surfcd2, surfcm, wingsc, wingpo, wingcl

include("../misc/index.jl")
print(igAR)
include("../../utils/spline.jl")
include("airtable2.jl")
include("airfun2.jl")
include("trefftz.jl")
include("cdsum.jl")
include("surfcd.jl")
include("surfcm.jl")
include("wingpo.jl")
include("wingsc.jl")

end