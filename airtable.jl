using Dierckx

"""
Reads airfoil file and outputs a matrix and spline objects
"""
function airtable(fname)

    # Read file
    f = open(fname)

    nAMa::Int, nAcl::Int, nAτ::Int, nAfun::Int = [parse(Int,ss) for ss in split(readline(f))]  

    AMa = zeros(nAMa)
    Acl = zeros(nAcl)
    Aτ  = zeros(nAτ)
    A   = zeros(nAMa, nAcl, nAτ, nAfun)

    for i = 1:nAMa
       AMa[i] = parse(Float64, readline(f))
    end
    for i=1:nAcl
        Acl[i] = parse(Float64, readline(f))
    end
    for i=1:nAτ
        Aτ[i] = parse(Float64, readline(f))
    end

    #Read in Reynolds number
    ARe = parse(Float64, readline(f))

    #Read in airfoil data
    for k = 1:nAτ
        for j = 1:nAcl
            for i = 1:nAMa
                    A[i,j,k,:] = [parse(Float64,ss) for ss in split(readline(f))]
            end
        end
    end

    close(f)

    # Create spline objects to make it faster for airfun.jl since this does not need to be computed each time

    ∂cdf_∂M  = Array{Spline1D, 2}(undef, nAcl, nAτ)
    ∂cdp_∂M  = Array{Spline1D, 2}(undef, nAcl, nAτ)
     ∂cm_∂M  = Array{Spline1D, 2}(undef, nAcl, nAτ)

    for k = 1:nAτ
        for j = 1:nAcl
            ∂cdf_∂M[j, k] = Spline1D(AMa, A[:,j,k,1], k = 3)
            ∂cdp_∂M[j, k] = Spline1D(AMa, A[:,j,k,2], k = 3)
             ∂cm_∂M[j, k] = Spline1D(AMa, A[:,j,k,3], k = 3)
        end
    end

return nAMa, nAcl, nAτ, nAfun, AMa, Acl, Aτ, A, ∂cdf_∂M, ∂cdp_∂M, ∂cm_∂M
end
