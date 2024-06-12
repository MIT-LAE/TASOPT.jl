function double_bubble_geom(Rfuse, dRfuse, wfb, nfweb, R = Rfuse)
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

    A = Afuse * R^2/Rfuse^2
    p = p_fuse * R/Rfuse
    lweb = lfuse_web * R/Rfuse

    return p, A, lweb
end