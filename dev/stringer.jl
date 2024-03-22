using TASOPT

tanktype = "inner"
W = 327.5e3
θ1 = 80*pi/180
Rtank = 1.525
s_a = 129.2e6
ρstiff = 2840

if tanktype == "inner" 
    _, kmax = TASOPT.structures.stiffeners_bendingM(θ1) #Find k = 2πM/(WR)
    Icollapse = 0 #Inner tank cannot collapse as it is pressurized

elseif tanktype == "outer"
    _, kmax = stiffeners_bendingM_outer(θ1, θ2) #Find k = 2πM/(WR)
    pc = 4 * pref #Critical pressure is 4 times atmospheric pressure, Eq. (7.11) in Barron (1985)
    Do = 2 * Rtank #outer diameter

    L = l_cyl/ (Nstiff - 1) #Length of portion between supports
    Icollapse = pc * Do * L / (24 * E) #Second moment of area needed to avoid collapse
end

Mmax = kmax * W * Rtank / (2π) #Maximum bending moment due to loads
Z = Mmax / s_a #required section modulus to withstand bending loads

#Assume sectional properties of a 100 x 100 I-beam
W = 80e-3 #Flange width
t_w = 7.1e-3 #Web thickness
t_f = 8.8e-3 #Flange thickness

#For an I-beam, I > W * H^2 * t_f / 2 + t_f^3 * W / 6
#The required second moment of area is I = Icollapse + Z * (H/2 + t_f/2)
#Find beam height by solving W * H^2 * t_f / 2 + t_f^3 * W / 6 = Icollapse + Z * (H/2 + t_f/2)

a = t_f * W / 2 #Coefficients in quadratic equation
b = -Z/2
c = -1 * Icollapse - Z * t_f/2 + t_f^3 * W / 6

H = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a) #Solve quadratic eq.
S = 2 * W * t_f + (H - t_f) * t_w #Beam cross-sectional area

Wstiff = gee * ρstiff * S * 2π * Rtank #Weight of a single stiffener running along circumference