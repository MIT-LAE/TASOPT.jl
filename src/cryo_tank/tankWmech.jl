"""
Length of the head-morphing window of the inner vessel, as a fraction of the (maximum) inner vessel radius. Over this window of
added cylinder length the heads change continuously from hemispherical (aspect ratio 1) to the design aspect ratio `ARtank`.
"""
const HEAD_MORPH_LEN_FRAC = 1.0

"""Smallest allowed size of the spherical inner vessel, as a fraction of its maximum radius (avoids a degenerate zero-volume tank)."""
const SPHERE_MIN_FRAC = 1e-3

"""
Geometric and weight outputs of the inner cryogenic tank sizing.
"""
struct InnerTankGeometry{T<:Real}
      """Total tank weight (structural + insulation) (N)"""
      Wtank::T
      """Total insulation weight (N)"""
      Winsul_sum::T
      """Volume of fuel (m^3)"""
      Vfuel::T
      """Insulated surface area of each ellipsoidal-cap interface, length N+1 (m^2)"""
      Shead_insul::Vector{T}
      """Outer radius of the inner vessel (m)"""
      Rtank_outer::T
      """Total longitudinal length of the tank (m)"""
      l_tank::T
      """Length of the cylindrical portion (m)"""
      l_cyl::T
      """Effective aspect ratio of the heads (-), between 1 (hemispherical) and the design `ARtank`"""
      AR_head::T
end

#Outer constructor that promotes mixed argument types (e.g., ForwardDiff duals and Float64)
function InnerTankGeometry(Wtank, Winsul_sum, Vfuel, Shead_insul::AbstractVector, Rtank_outer, l_tank, l_cyl, AR_head)
      T = promote_type(typeof(Wtank), typeof(Winsul_sum), typeof(Vfuel), eltype(Shead_insul), 
                       typeof(Rtank_outer), typeof(l_tank), typeof(l_cyl), typeof(AR_head))
      return InnerTankGeometry{T}(T(Wtank), T(Winsul_sum), T(Vfuel), convert(Vector{T}, Shead_insul), 
                                  T(Rtank_outer), T(l_tank), T(l_cyl), T(AR_head))
end

"""
      size_inner_tank(fuse::Fuselage, fuse_tank::fuselage_tank, t_cond::Vector{Float64})

`size_inner_tank` calculates the weight of the cryogenic fuel tank for a LH-fueled aircraft.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuse::Fuselage`: fuselage object.
      - `fuse_tank::fuselage_tank`: fuselage tank object.
      - `t_cond::Vector{Float64}`: Vector with tank insulation layer thickness. Provided separately from fuse_tank as it changes during
      non-linear solve process.

      **Outputs:**
      - `::InnerTankGeometry`: inner tank geometry and weights.

NOTE: Al alloy 2219 has been recommended as tank material (from H2 tank paper in OneNote)

See [here](@ref fueltanks).
"""
function size_inner_tank(fuse::Fuselage, fuse_tank::fuselage_tank, t_cond::AbstractVector{<:Real})
      #---------------------------------
      # Unpack parameters in fuse_tank
      #---------------------------------
      Rfuse = fuse.layout.radius
      fuse_cs = fuse.layout.cross_section #Fuselage cross section

      Wfuel = fuse_tank.Wfuelintank

      ρfuel = fuse_tank.rhofuel
      ρfgas = fuse_tank.rhofuelgas
      ftankadd = fuse_tank.ftankadd
      Δp = fuse_tank.pvent
      sigskin = fuse_tank.inner_material.UTS
      material_insul = fuse_tank.material_insul
      rhoskin = fuse_tank.inner_material.ρ
      clearance_fuse = fuse_tank.clearance_fuse
      AR = fuse_tank.ARtank
      ullage_frac = fuse_tank.ullage_frac
      weld_eff = fuse_tank.ew
      θ = fuse_tank.theta_inner

      #---------------------------------
      # Wall thickness sizing
      #---------------------------------
      thickness_insul = sum(t_cond) #Total insulation thickness
      Rmax_outer = Rfuse - thickness_insul - clearance_fuse #Largest inner vessel outer radius that fits in the fuselage
      Rmax_outer > 0 || throw(DomainError(Rmax_outer, "size_inner_tank: insulation and clearance leave no room for the inner vessel"))

      #Wall thicknesses scale linearly with the vessel radius, so the ratio t/R is a constant for a given pressure
      s_a = sigskin / 4 #Maximum allowable stress is 1/4 Ultimate tensile strength (Barron 1985, p. 359)
      tfrac_cyl = Δp * 2 / (2 * s_a * weld_eff + 0.8 * Δp) #Cylinder wall thickness / outer radius, (7.1) in Barron (1985)
      Rmax_inner = Rmax_outer * (1 - tfrac_cyl) #Largest inner vessel inner radius

      #---------------------------------
      # Required volume
      #---------------------------------
      Wfuel_tot = Wfuel #Wfuel already includes the amount that boils off
      ρfmix = (1 - ullage_frac) * ρfuel + ullage_frac * ρfgas #Density of saturated mixture in tank
      Vfuel = Wfuel_tot / (gee * ρfmix) #Total tank volume taken by saturated mixture
      Vinternal = Vfuel  # required internal volume

      #---------------------------------
      # Continuous tank shape as a function of volume
      #---------------------------------
      #The volume of two end caps of aspect ratio AR on a cross-section of area A and radius R is 2*(2/3)*A*R/AR 
      #(for a standard ellipsoid V = 1/2×(4π/3 ×(abc)) where a,b,c are the semi-axes length. If base is circular,
      #a = R/AR, b=c=R. Since Abase = πR^2, this can also be written as V = 2 * (Abase * R/AR)/3. 
      #
      #The shape grows without discontinuities as the fuel volume increases:
      #  1. V ≤ V_sph: a sphere (two AR = 1 caps, no cylinder), geometrically scaled from zero up to the maximum 
      #     radius that fits in the fuselage (insulation and clearance included).
      #  2. V_sph < V ≤ V_morph: radius is the maximum radius, a cylinder begins to grow, and the aspect ratio of the heads 
      #     moves smoothly (cubic smoothstep) from 1 (hemispherical) to the design value ARtank.
      #  3. V > V_morph: design heads and a cylinder with length l_cyl = (V - V_caps)/A.
      _, Amax = scaled_cross_section(fuse_cs, Rmax_inner)
      AR_design = max(AR, one(AR)) #Heads are never more slender than hemispherical
      V_sph = 4 * Amax * Rmax_inner / 3 #Volume of the largest sphere (two AR = 1 caps)
      V_morph = V_sph + HEAD_MORPH_LEN_FRAC * Rmax_inner * Amax #End of the morphing window

      if Vinternal <= V_sph #Regime 1: sphere
            scale = max(cbrt(Vinternal / V_sph), SPHERE_MIN_FRAC) #Geometric scaling of the largest sphere
            Rtank_outer = Rmax_outer * scale
            AR_head = one(AR_design)
            l_cyl = zero(Rtank_outer)
      else #Regimes 2 and 3
            Rtank_outer = Rmax_outer
            ξ = clamp((Vinternal - V_sph) / (V_morph - V_sph), 0.0, 1.0)
            AR_head = 1 + (AR_design - 1) * ξ^2 * (3 - 2 * ξ) #Smoothstep: zero slope at both ends of the window
            Vcaps = 4 * Amax * Rmax_inner / (3 * AR_head)
            l_cyl = max((Vinternal - Vcaps) / Amax, zero(Vinternal / Amax)) #Required length of cylindrical portion
      end

      #Cylindrical portion wall thickness
      tskin = tfrac_cyl * Rtank_outer 
      Rtank = Rtank_outer - tskin #Inner vessel inner radius

      #Ellipsoid wall thickness
      Lhead = Rtank / AR_head       # eg. for a 2:1 ellipsoid majorax/minorax = 2/1 ⟹ R/Lhead = 2/1 
      K = (1/6) * (AR_head^2 + 2) # Barron pg 359
      t_head = Δp* (2*Rtank_outer) * K/ (2 * s_a * weld_eff + 2 * Δp * (K - 0.1)) #(7.2) in Barron (1985)

      #---------------------------------
      # Vessel geometry
      #---------------------------------
      perim_tank, Atank = scaled_cross_section(fuse_cs, Rtank) #Tank perimeter and cross-sectional area

      # areas
      Scyl = perim_tank*l_cyl # Surface area of cylindrical part

      Shead = 2*Atank * ( 0.333 + 0.667*(Lhead/Rtank)^1.6 )^0.625 # Surface area of ellipsoidal caps
            #This form is better for insul thickness
            # but just as a note to reader this comes from  semi- oblate spheroid surf area is ≈ 2π×R²[1/3 + 2/3×(1/AR)^1.6]^(1/1.6)
      
            # component volumes    
      Vcyl  = Scyl*tskin    # volume of the metal in the cylindrical part

      Vhead = Shead * t_head # volume of head

      Sinternal = Scyl + 2 * Shead

      #---------------------------------
      # Weight calculation
      #---------------------------------
      Whead = rhoskin*gee*Vhead #Head weight
      Wcyl  = rhoskin*gee*Vcyl #Cylinder weight
      Winnertank = Wcyl + 2*Whead + Wfuel #Weight of inner vessel without stiffeners and supports, inc. fuel

      #stiffeners
      Nmain = 2 #Tanks typically have two main support rings
      Wmainstiff = stiffener_weight("inner", Winnertank / Nmain, Rtank, perim_tank, 
                                s_a, rhoskin, θ) #Weight of one main stiffener

      Wstiff = Nmain * Wmainstiff
      Wtank = (Wcyl + 2*Whead + Wstiff) * (1 + ftankadd)

      # Insulation weight
      N = length(t_cond) #Number of insulation layers
      #Initialize storage vectors
      Winsul_sum = zero(eltype(t_cond)) #Element type follows t_cond so that ForwardDiff duals propagate
      Shead_insul = zeros(eltype(t_cond), N + 1) #add one for first (tank wall) surface 
      L = Lhead + tskin #Length of ellipsoid semi-minor axis

      Ro = Ri = Rtank_outer # Start calculating insulation from the outer wall of the metal tank ∴Ri of insul = outer R of tank
      _, Ao = scaled_cross_section(fuse_cs, Ro) #Cross-sectional area of double bubble
      Shead_insul[1] = 2*Ao * ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625 #Surface area of vessel head
      
      for n in 1:N
            #Update radius and semi-minor axis
            Ro = Ro + t_cond[n]
            L  = L  + t_cond[n]

            _, Ao = scaled_cross_section(fuse_cs, Ro)
            _, Ai = scaled_cross_section(fuse_cs, Ri)
            rho_insul = material_insul[n].ρ
            Vcyl_insul  = l_cyl * (Ao - Ai) #Volume of cylindrical layer
            Shead_insul[n+1] = 2*Ao * ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625 #Surface area of ellipsodal cap

            Area_coeff = Shead_insul[n+1] / Ro^2 #coefficient that relates area and radius squared
            Vhead_insul = ((Shead_insul[n] + Shead_insul[n+1])/2 - Area_coeff/(6) * t_cond[n]^2) * t_cond[n] #Closed-form solution
            
            Winsul_sum += (Vcyl_insul + 2*Vhead_insul) * rho_insul * gee #Weight of insulation layer

            Ri = Ro #Update inner radius for next point
      end

      Wtank = (Wtank + Winsul_sum)
      l_tank = l_cyl + 2*Lhead + 2*thickness_insul + 2*t_head #Total longitudinal length of the tank

return InnerTankGeometry(Wtank, Winsul_sum, Vfuel, Shead_insul, Rtank_outer, l_tank, l_cyl, AR_head)
end


"""
Geometric and weight outputs of the outer cryogenic tank (vacuum-jacket vessel) sizing.
"""
struct OuterTankGeometry
      """Total weight of outer vessel (N)"""
      Wtank::Float64
      """Weight of cylindrical portion (N)"""
      Wcyl::Float64
      """Weight of one elliptical head (N)"""
      Whead::Float64
      """Total stiffener weight (N)"""
      Wstiff::Float64
      """Total surface area of outer vessel (m^2)"""
      Souter::Float64
      """Surface area of one outer-vessel head (m^2)"""
      Shead::Float64
      """Surface area of cylindrical portion (m^2)"""
      Scyl::Float64
      """Wall thickness of cylindrical portion (m)"""
      t_cyl::Float64
      """Wall thickness of tank head (m)"""
      t_head::Float64
      """Total length of the outer vessel (m)"""
      l_outer::Float64
end

"""
    size_outer_tank(fuse::Fuselage, fuse_tank::fuselage_tank, Winnertank::Float64, l_cyl::Float64, Ninterm::Float64)
This function sizes the outer vessel and calculates the weights of its components.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse::Fuselage`: fuselage object.
    - `fuse_tank::fuselage_tank`: fuselage tank object.
    - `Winnertank::Float64`: weight of inner vessel and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of outer vessel (m).
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.
    - `Router::Float64`: radius of the outer vessel (m). 
    - `AR_head::Float64`: aspect ratio of the outer vessel heads.

    **Outputs:**
    - `::OuterTankGeometry`: outer tank geometry and weights.
"""
function size_outer_tank(fuse::Fuselage, fuse_tank::fuselage_tank, Winnertank::Float64, l_cyl::Float64, Ninterm::Float64;
      Router::Float64 = fuse.layout.radius - fuse_tank.clearance_fuse, AR_head::Float64 = fuse_tank.ARtank)
      #---------------------------------
      # Unpack parameters in fuse_tank
      #---------------------------------
      Rfuse = fuse.layout.radius
      fuse_cs = fuse.layout.cross_section #Fuselage cross section

      poiss = fuse_tank.inner_material.ν
      Eouter = fuse_tank.inner_material.E
      ρouter = fuse_tank.inner_material.ρ
      UTSouter = fuse_tank.inner_material.UTS
      ftankadd = fuse_tank.ftankadd
      ARtank = AR_head
      θ_outer = fuse_tank.theta_outer

      #---------------------------------
      # Size outer vessel
      #---------------------------------
      θ1 = θ_outer[1]
      θ2 = θ_outer[2]

      Nmain = 2 #Tanks typically have two main support rings
      pc = 4 * pref #4*p_atm; Collapsing pressure, Eq. (7.11) in Barron (1985)
      s_a = UTSouter / 4

      #Calculate outer vessel geometry
      Rtank_outer = Router
      Do = 2 * Rtank_outer #outside diameter

      Nstiff = Nmain + Ninterm #Total number of stiffeners
      L = l_cyl / (Nstiff - 1) #There are two stiffeners at the ends, so effective number of sections in skin is N - 1
      L_Do = L / Do

      #---------------------------------
      # Wall thicknesses
      #---------------------------------

      #Find cylinder wall thickness. This applies to a short cylinder.
      #If the vessel has no cylindrical portion (spherical inner vessel) there is no cylinder to size; this is also the limit of
      #the buckling relation as L/D → 0, so the weight is continuous.
      t_Do = l_cyl > 0 ? outer_cylinder_thickness_ratio(Eouter, poiss, L_Do, pc) : 0.0 #Bracketed root of the buckling equation
      t_cyl = t_Do * Do

      #Find head wall thickness
      K1 = find_K1_head(ARtank) #Find coefficient K1
      t_head = K1 * Do * sqrt(pc * sqrt(3*(1 - poiss^2))/ (0.5*Eouter))

      #---------------------------------
      # Vessel geometry
      #---------------------------------
      perim_vessel, Avessel = scaled_cross_section(fuse_cs, Rtank_outer) #Tank perimeter and cross-sectional area

      Shead = 2*Avessel*(0.333 + 0.667*(1/ARtank)^1.6 )^0.625 #ellipsoid surface area
      Scyl  = perim_vessel*l_cyl  # surface area

      Souter = Scyl + 2*Shead

      ## Volume and Weight
      Vcyl  = Scyl*t_cyl
      Vhead = Shead*t_head

      Wcyl  = Vcyl*ρouter*gee
      Whead =  Vhead*ρouter*gee

      Wtank_no_stiff = Wcyl + 2 * Whead

      #---------------------------------
      # Size stiffeners
      #---------------------------------
      tanktype = "outer"

      Wmainstiff = stiffener_weight(tanktype, Winnertank / Nmain, Rtank_outer, perim_vessel, #Weight of one main stiffener, each one 
                                    s_a, ρouter, θ1, θ2, Nstiff, l_cyl, Eouter)  #carries half the inner vessel load
                                                                              
      Wintermstiff = stiffener_weight(tanktype, 0.0, Rtank_outer, perim_vessel, 
                                    s_a, ρouter, θ1, θ2, Nstiff, l_cyl, Eouter) #Weight of one intermediate stiffener, which carries no load

      #---------------------------------
      # Weight estimation
      #---------------------------------
      Wstiff = Nmain * Wmainstiff + Ninterm * Wintermstiff #Total stiffener weight

      Wtank = (Wtank_no_stiff + Wstiff) * (1 + ftankadd) #Find total tank weight, including additional mass factor

      l_outer = l_cyl + Do / ARtank + 2*t_head #Total length of outer vessel

      return OuterTankGeometry(Wtank, Wcyl, Whead, Wstiff, Souter, Shead, Scyl, t_cyl, t_head, l_outer)
end

"""
      outer_cylinder_thickness_ratio(E, ν, L_Do, pc)

Finds the thickness-to-diameter ratio `t/D` of a short cylindrical outer vessel that resists a collapsing pressure `pc`, 
using the buckling relation of Barron (1985):

``p_c = \\frac{2.42 E (t/D)^{5/2}}{(1-ν^2)^{3/4} (L/D - 0.45 \\sqrt{t/D})}``

The right-hand side is `-pc` at `t/D = 0`, increases monotonically, and diverges at `t/D = (L/D / 0.45)^2`, so a single root
always exists in between. It is found with a bracketed bisection, which cannot step outside the valid interval 
(unlike an unbracketed secant iteration).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `E::Float64`: Young's modulus of the outer vessel (Pa).
    - `ν::Float64`: Poisson ratio of the outer vessel.
    - `L_Do::Float64`: ratio of unsupported length to outer diameter. Must be positive and finite.
    - `pc::Float64`: collapsing pressure (Pa).

    **Outputs:**
    - `t_Do::Float64`: wall thickness to outer diameter ratio.
"""
function outer_cylinder_thickness_ratio(E::Real, ν::Real, L_Do::Real, pc::Real)
      (isfinite(L_Do) && L_Do > 0) || throw(DomainError(L_Do, "outer_cylinder_thickness_ratio: L/D must be positive and finite"))

      t_pole = (L_Do / 0.45)^2 #Location of the pole of the buckling relation
      pressure_res(t_D) = 2.42*E*t_D^(5/2) / ( (1 - ν^2)^(3/4) * (L_Do - 0.45*sqrt(t_D)) ) - pc

      #Stay slightly away from the pole, where the residual is infinite. The residual is -pc at 0 and positive near the pole.
      t_hi = t_pole * (1 - 1e-9)
      return find_zero(pressure_res, (zero(t_hi), t_hi), Bisection())
end

"""
      stiffener_weight(tanktype::String, W::Float64, Rtank::Float64, perim::Float64, s_a::Float64, 
      ρstiff::Float64, θ1::Float64, θ2::Float64 = 0.0, Nstiff::Float64 = 2.0, l_cyl::F
This function calculates the weight of a single stiffener in an inner or outer vessel for a given inner vessel weight.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `tanktype::String`: type of tank, options are "inner" or "outer".
    - `W::Float64`: load carried by a stiffener ring (N).
    - `Rtank::Float64`: tank radius (m).
    - `perim::Float64`: tank perimeter (m).
    - `s_a::Float64`: maximum allowable stress in stiffener material (Pa).
    - `ρstiff::Float64`: stiffener density (kg/m^3).
    - `θ1::Float64`: angular position of bottom tank supports, measured from the bottom of the tank (rad).
    - `θ2::Float64`: angular position of top tank supports, measured from the bottom of the tank (rad). Only used with "outer" tank.
    - `Nstiff::Float64`: total number of stiffeners on outer vessel. Only used with "outer" tank.
    - `l_cyl::Float64`: length of cylindrical portion of tank (m). Only used with "outer" tank.
    - `E::Float64`: Young's modulus of stiffener material (Pa). Only used with "outer" tank.

    **Outputs:**
    - `Wstiff::Float64`: weight of a stiffener ring (N).
"""
function stiffener_weight(tanktype::String, W::Real, Rtank::Real, perim::Real, s_a::Float64, 
      ρstiff::Float64, θ1::Float64, θ2::Float64 = 0.0, Nstiff::Float64 = 2.0, l_cyl::Float64 = 0.0, E::Float64 = 0.0)
    
      if tanktype == "inner" 
            _, kmax = stiffeners_bendingM(θ1) #Find k = 2πM/(WR)
            Icollapse = 0.0 #inner vessel cannot collapse as it is pressurized

      elseif tanktype == "outer"
            _, kmax = stiffeners_bendingM_outer(θ1, θ2) #Find k = 2πM/(WR)
            pc = 4 * pref #Critical pressure is 4 times atmospheric pressure, Eq. (7.11) in Barron (1985)
            Do = 2 * Rtank #outer diameter

            L = l_cyl/ (Nstiff - 1) #Length of portion between supports
            Icollapse = pc * Do^3 * L / (24 * E) #Second moment of area needed to avoid collapse
      end

      Mmax = kmax * W * Rtank / (2π) #Maximum bending moment due to loads
      Z = Mmax / s_a #required section modulus to withstand bending loads

      #Assume sectional properties of a 100 x 100 I-beam
      W = 100e-3 #Flange width
      t_w = 7.1e-3 #Web thickness
      t_f = 8.8e-3 #Flange thickness

      #For an I-beam, I > W * H^2 * t_f / 2 + t_f^3 * W / 6
      #The required second moment of area is I = Icollapse + Z * (H/2 + t_f/2)
      #Find beam height by solving W * H^2 * t_f / 2 + t_f^3 * W / 6 = Icollapse + Z * (H/2 + t_f/2)

      a = t_f * W / 2 #Coefficients in quadratic equation
      b = -Z/2
      c = -1 * Icollapse - Z * t_f/2 + t_f^3 * W / 6

      #For a lightly loaded ring the flanges alone exceed the requirement and the discriminant can be negative.
      #In that case use the minimum physical beam height (flanges touching).
      disc = max(b^2 - 4 * a * c, zero(b^2 - 4 * a * c))
      H = max((-b + sqrt(disc)) / (2 * a), 2 * t_f) #Solve quadratic eq.
      S = 2 * W * t_f + (H - t_f) * t_w #Beam cross-sectional area

      Wstiff = gee * ρstiff * S * perim #Weight of a single stiffener running along perimeter
      return Wstiff
end

"""
    stiffeners_bendingM(θ)
This function can be used to calculate the bending moment distribution in a stiffener ring for an inner cryogenic tank.
It applies Eqs. (7.4) and (7.5) in Barron (1985) to find the bending moment distribution. The function returns the
maximum value of ``k = 2πM/(WR)`` on the ring's circumference.

Note: An attempt was made to use the NLopt.jl package to find the maximum value of k, 
but this was more expensive and less robust than the brute-force approach.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `θ::Float64`: angular position of tank supports, measured from the bottom of the tank (rad).

    **Outputs:**
    - `ϕmax::Float64`: angular position of maximum bending moment on ring circumference (rad).
    - `kmax::Float64`: Maximum value of the ratio ``k = 2πM/(WR)`` on ring circumference.
"""
function stiffeners_bendingM(θ::Float64)
      ϕlist = [θ; LinRange(1.1,1.2,21);pi] #Take advantage of known form of maximum
      k = zeros(length(ϕlist))

      for (i,ϕ) in enumerate(ϕlist)
            if 0 ≤ ϕ < θ
                  k[i] = 0.5*cos(ϕ) + ϕ*sin(ϕ) - (π - θ)*sin(θ) + cos(θ) + cos(ϕ)*(sin(θ)^2)
            elseif θ ≤ ϕ ≤ π
                  k[i] = 0.5*cos(ϕ) - (π - ϕ)*sin(ϕ) + θ  + cos(θ) + cos(ϕ)*(sin(θ)^2)
                  #TODO this equation is Eq. (7.5) in Barron; however, this equation does not match the curves in Fig. 7.3
                  #Suspect error in Barron.
            end
      end
      kmax, imax = findmax(abs.(k))
      ϕmax = ϕlist[imax]
      return ϕmax, kmax
end

"""
    stiffeners_bendingM_outer(θ1,θ2)
This function can be used to calculate the bending moment distribution in a stiffener ring for an outer cryogenic tank.
It applies Eqs. (7.13)-(7.15) in Barron (1985) to find the bending moment distribution. The function returns the
maximum value of ``k = 2πM/(WR)`` on the ring's circumference.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `θ1::Float64`: angular position of bottom tank supports, measured from the bottom of the tank (rad).
    - `θ2::Float64`: angular position of top tank supports, measured from the bottom of the tank (rad).

    **Outputs:**
    - `ϕmax::Float64`: angular position of maximum bending moment on ring circumference (rad).
    - `kmax::Float64`: Maximum value of the ratio ``k = 2πM/(WR)`` on ring circumference.
"""
function stiffeners_bendingM_outer(θ1::Float64, θ2::Float64)
      ϕlist = LinRange(0.0, π, 361)
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
      return ϕmax, kmax
end


"""
    optimize_outer_tank(fuse::Fuselage, fuse_tank::fuselage_tank, Winnertank, l_cyl, θ1, θ2)
This function optimizes the number of intermediate stiffener rings to minimize the weight of an outer vessel.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse::Fuselage`: fuselage object.
    - `fuse_tank::fuselage_tank`: fuselage tank object.
    - `Winnertank::Float64`: weight of inner vessel and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of outer vessel (m).

    **Outputs:**
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.
"""
function optimize_outer_tank(fuse::Fuselage, fuse_tank::fuselage_tank, Winnertank::Float64, l_cyl::Float64;
      Router::Float64 = fuse.layout.radius - fuse_tank.clearance_fuse, AR_head::Float64 = fuse_tank.ARtank)

      obj(x, grad) = size_outer_tank(fuse, fuse_tank, Winnertank, l_cyl, x[1]; Router, AR_head).Wtank #Minimize Wtank

      #Set bounds
      lower = 0.0
      upper = 50.0

      #Start strictly inside the bounds: a Nelder-Mead simplex started on the lower bound (e.g. from a previous 
      #solution of Ninterm = 0) can collapse there and never explore larger values.
      initial_x = [clamp(fuse_tank.Ninterm, lower .+ 0.5, upper .- 0.5)]
      
      #Use NLopt.jl to minimize function 
      opt = Opt(:LN_NELDERMEAD, length(initial_x))
      opt.lower_bounds = lower
      opt.upper_bounds = upper
      opt.ftol_rel = 1e-9
      opt.maxeval = 100  # Set the maximum number of function evaluations

      opt.min_objective = obj

      (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
      Ninterm = xopt[1]
      return Ninterm
end




"""
    optimize_outer_tank(fuse::Fuselage, fuse_tank::fuselage_tank, Winnertank, l_cyl; xtol = 0.05, halfwidth = 3.0, maxexpand = 6)
This function optimizes the number of intermediate stiffener rings to minimize the weight of an outer vessel.
It uses a bounded golden-section search in a bracket centered on the previous solution (`fuse_tank.Ninterm`),
which is shifted outward if the minimum lands on a bracket edge. Typically requires ~10-15 evaluations of 
`size_outer_tank` instead of up to 100 with Nelder-Mead.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse::Fuselage`: fuselage object.
    - `fuse_tank::fuselage_tank`: fuselage tank object.
    - `Winnertank::Float64`: weight of inner vessel and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of outer vessel (m).
    - `xtol::Float64`: tolerance on `Ninterm` for the golden-section search.
    - `halfwidth::Float64`: half-width of the initial bracket around the warm-start guess.
    - `maxexpand::Int`: maximum number of times the bracket is shifted outward.

    **Outputs:**
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.
"""
function optimize_outer_tank_new(fuse::Fuselage, fuse_tank::fuselage_tank, Winnertank::Float64, l_cyl::Float64;
      xtol::Float64 = 0.05, halfwidth::Float64 = 3.0, maxexpand::Int = 6,
      Router::Float64 = fuse.layout.radius - fuse_tank.clearance_fuse, AR_head::Float64 = fuse_tank.ARtank)

      ### This is not nessarily faster than previous iteration - both bug out at the same time within optimize iteration.
      obj(N) = size_outer_tank(fuse, fuse_tank, Winnertank, l_cyl, N; Router, AR_head).Wtank #Minimize Wtank

      #Global bounds on the number of intermediate stiffeners
      Nmin, Nmax = 0.0, 50.0

      #Warm-started bracket centered on the previous solution
      x0 = clamp(fuse_tank.Ninterm, Nmin, Nmax)
      lo = max(Nmin, x0 - halfwidth)
      hi = min(Nmax, x0 + halfwidth)

      invϕ = (sqrt(5.0) - 1) / 2 #1/golden ratio
      xopt = x0
      for _ in 0:maxexpand
            #Golden-section search on [lo, hi]
            a, b = lo, hi
            c = b - invϕ * (b - a)
            d = a + invϕ * (b - a)
            fc, fd = obj(c), obj(d)
            while (b - a) > xtol
                  if fc < fd
                        b, d, fd = d, c, fc
                        c = b - invϕ * (b - a)
                        fc = obj(c)
                  else
                        a, c, fc = c, d, fd
                        d = a + invϕ * (b - a)
                        fd = obj(d)
                  end
            end
            xopt = (a + b) / 2

            #If the minimum lies at a bracket edge that is not a global bound, the true optimum may be
            #outside the bracket: shift/expand the bracket in that direction and search again
            at_lo = (xopt - lo) < 2xtol && lo > Nmin
            at_hi = (hi - xopt) < 2xtol && hi < Nmax
            if at_lo
                  hi = lo
                  lo = max(Nmin, lo - 2halfwidth)
            elseif at_hi
                  lo = hi
                  hi = min(Nmax, hi + 2halfwidth)
            else
                  break
            end
      end
      return xopt
end
"""
      find_K1_head(AR)

This function calculates the equivalent radius for hemi-ellipsoidal heads under pressure.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `AR::Float64`: head aspect ratio.

      **Outputs:**
      - `K1::Float64`: equivalent radius parameter, such that `R = K1 * D``, where `R` is 
      the equivalent radius and `D` is the ellipse major diameter.
"""
function find_K1_head(AR)
      # See table 7.6 for data in Barron (1985) p. 367
      ARs = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
      K1s = [0.50,0.57,0.65,0.73,0.81,0.9,0.99,1.08,1.18,1.27,1.36]
      if (AR > maximum(ARs)) || (AR < minimum(ARs)) #If case is outside of range
            println("ARtank of heads not supported, see size_outer_tank()")
            K1 = 1.0
      else
            i = 1
            while AR > ARs[i] #Find index of entry that is not smaller than AR
                  i += 1
            end
            if AR ≈ ARs[i] #If it is exactly that entry
                  K1 = K1s[i]
            else #Otherwise interpolate linearly
                  K1 = K1s[i-1] + (K1s[i] - K1s[i-1]) / (ARs[i] - ARs[i-1]) * (AR - ARs[i-1])
            end
      end
      return K1
end