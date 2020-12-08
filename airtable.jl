"""
Reads airfoil file
"""
function airtable(fname)

f = open("./air/C.air")
    
nAMa, nAcl, nAτ, nAfun = [parse(Int,ss) for ss in split(readline(f))]  

AMa = zeros(nAMa)
Acl = zeros(nAcl)
Aτ = zeros(nAτ)
A = zeros(nAMa, nAcl, nAτ, nAfun)

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
return nAMa, nAcl, nAτ, nAfun, AMa, Acl, Aτ, A
end
