using PyPlot

θ1 = 60*pi/180
θ2 = 120*pi/180

ϕlist = LinRange(0.0, π, 360)
k = zeros(length(ϕlist))

for (i,ϕ) in enumerate(ϕlist)
    if 0 ≤ ϕ ≤ θ1
            k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) ) +
                        - (π-θ2)*sin(θ2) +  (π-θ1)*sin(θ1)
    elseif θ1 ≤ ϕ ≤ θ2
            k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) ) +
                        - (π-θ2)*sin(θ2) +  π*sin(ϕ) -θ1*sin(θ1)

    elseif θ2 ≤ ϕ ≤ π
            k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) )  +
                        (θ2*sin(θ2) - θ1*sin(θ1)  )
    end

end
kmax, imax = findmax(abs.(k))
ϕmax = ϕlist[imax]

plt.figure()
plt.plot(ϕlist*180/pi,k)
plt.xlim(0,180)
plt.savefig("k.png")
plt.close()