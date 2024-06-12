"""
    double_bubble_geom(Rfuse::Float64, dRfuse::Float64, wfb::Float64, nfweb, R::Float64 = Rfuse)

Calculates the geometric properties of a double-bubble cross section, such as the aicraft's fuselage. In
addition to the fuselage geometry, this function can calculate properties for a scaled shape with a 
double-bubble radius different tho that of the fuselage.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Rfuse::Float64`: fuselage exterior radius (m)
    - `dRfuse::Float64`: downward shift of double bubble (m)
    - `wfb::Float64`: lateral shift of double bubble (m)
    - `nfweb::Float64`: number of vertical webs in fuselage
    - `R::Float64`: radius of geometrically-similar double bubble (m)

    **Outputs:**
    - `p::Float64`: perimeter
    - `A::Float64`: cross-sectional area (m^2)
    - `lweb::Float64`: length of vertical web (m)
"""
function double_bubble_geom(Rfuse::Float64, dRfuse::Float64, wfb::Float64, nfweb, R::Float64 = Rfuse)
    #cross-section geometric parameters
    wfblim = max( min( wfb , Rfuse) , 0.0 )
    thetafb = asin(wfblim / Rfuse)
    hfb = sqrt(Rfuse^2 - wfb^2)
    sin2t = 2.0*hfb*wfb/Rfuse^2 #sin(2θ) = 2sinθcosθ
    
    #Fuselage area and perimeter
    Afuse = (pi + nfweb*(2.0*thetafb + sin2t))*Rfuse^2 + 2.0*Rfuse*dRfuse #Fuselage cross-sectional area
    p_fuse = (2.0*pi+4.0*nfweb*thetafb)*Rfuse + 2.0*dRfuse #Fuselage perimeter

    #Web length
    lfuse_web = 2.0*hfb + dRfuse

    #Apply scaling: areas scale with R^2 and distances scale with R
    A = Afuse * R^2/Rfuse^2
    p = p_fuse * R/Rfuse
    lweb = lfuse_web * R/Rfuse

    return p, A, lweb
end