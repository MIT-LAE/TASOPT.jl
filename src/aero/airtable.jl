"""
    airtable(fname)

Reads airfoil file and outputs the matrix and splines as an `airfoil`.

The airfoil data is stored as a function of three variables, typically
Mach number ``\\mathrm{Ma}``, lift coefficient ``c_l``, and thickness to chord ratio ``\\tau``.

    cdf(Ma, cl, τ), cdp(Ma, cl, τ), cm(Ma, cl, τ)

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fname::String`: Path to file.

    **Outputs:**
    - `airf::TASOPT.aerodynamics.airfoil`: `struct` with airfoil performance characteristics.

"""
function airtable(fname)

# Read airfoil data
f = open(fname)

    nAMa::Int, nAcl::Int, nAτ::Int, nAfun::Int = [parse(Int,ss) for ss in split(readline(f))]  

    AMa = zeros(nAMa)     # array of Mach numbers in database
    Acl = zeros(nAcl)     # array of cl           in database     
    Aτ  = zeros(nAτ)      # array of t/c          in database    
    A   = zeros(nAMa, nAcl, nAτ, nAfun)

    @inbounds for  i = 1:nAMa
       AMa[i] = parse(Float64, readline(f))
    end
    @inbounds for  i=1:nAcl
        Acl[i] = parse(Float64, readline(f))
    end
    @inbounds for  i=1:nAτ
        Aτ[i] = parse(Float64, readline(f))
    end

    # Read in Reynolds number
    ARe = parse(Float64, readline(f))

    #Read in airfoil data
    @inbounds for  k = 1:nAτ
        @inbounds for  j = 1:nAcl
            @inbounds for  i = 1:nAMa
                    A[i,j,k,:] = [parse(Float64,ss) for ss in split(readline(f))]
            end
        end
    end

close(f)

# ---------------------
     A_M = zeros(nAMa, nAcl, nAτ, nAfun)
     A_τ = zeros(nAMa, nAcl, nAτ, nAfun)
    A_cl = zeros(nAMa, nAcl, nAτ, nAfun)
   A_M_τ = zeros(nAMa, nAcl, nAτ, nAfun)
  A_M_cl = zeros(nAMa, nAcl, nAτ, nAfun)
  A_cl_τ = zeros(nAMa, nAcl, nAτ, nAfun)
A_M_cl_τ = zeros(nAMa, nAcl, nAτ, nAfun)
# Iterate over the 3 functions - cdf, cdp and cm    
for l = 1:nAfun
    for k = 1:nAτ
        for j = 1:nAcl
            A_M[:, j, k, l] = spline(AMa, A[:, j, k, l]) # This gets the derivatives ∂cdf/∂M, ∂cdp/∂M, ∂cm/∂M
        end
    end

    for k = 1:nAτ
        for i = 1:nAMa
            # for j = 1:nAcl
                f   =   A[i, :, k, l]
                f_M = A_M[i, :, k, l]
                  A_cl[i, :, k, l] = spline(Acl, f)
                A_M_cl[i, :, k, l] = spline(Acl, f_M)
        end
    end

    for j = 1:nAcl
        for i = 1:nAMa
            f   = A[i, j, :, l]
            f_M = A_M[i, j, :, l]
              A_τ[i, j, :, l] = spline(Aτ, f)
            A_M_τ[i, j, :, l] = spline(Aτ, f_M)
        end
    end

    for i = 1:nAMa
        for j = 1:nAcl
            f   = A_cl[i, j, :, l]
            f_M = A_M_cl[i, j, :, l]
              A_cl_τ[i, j, :, l] = spline(Aτ, f)
            A_M_cl_τ[i, j, :, l] = spline(Aτ, f_M)
        end
        for k = 1:nAτ
            f   = A_τ[i, :, k, l]
            f_M = A_M_τ[i, :, k, l]
              A_cl_τ[i, :, k, l] = 0.5*(  A_cl_τ[i, :, k, l] + spline(Acl, f))
            A_M_cl_τ[i, :, k, l] = 0.5*(A_M_cl_τ[i, :, k, l] + spline(Acl, f_M))
        end
    end


end

return  airfoil(AMa, Acl, Aτ, ARe,
        A,
        A_M,
        A_τ,
        A_cl,
        A_M_τ,
        A_M_cl,
        A_cl_τ,
        A_M_cl_τ)
end

