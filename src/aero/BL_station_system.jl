using StaticArrays

"""
    _BL_station_system(is_selfsimilar, is_laminar, is_wake, solves_direct, Mach, uinv, hksep,
          x, b, rn, th, ds, ue,
          h , h_th, h_ds,
          hk, hk_th, hk_ds, hk_ue,
          hc, hc_th, hc_ds, hc_ue,
          hs, hs_th, hs_ds, hs_ue,
          cf, cf_th, cf_ds, cf_ue,
          di, di_th, di_ds, di_ue,
          xm,bm,rnm,thm,dsm,uem, 
          hm , hm_thm, hm_dsm,
          hkm, hkm_thm, hkm_dsm, hkm_uem,
          hcm, hcm_thm, hcm_dsm, hcm_uem,
          hsm, hsm_thm, hsm_dsm, hsm_uem,
          cfm, cfm_thm, cfm_dsm, cfm_uem,
          dim, dim_thm, dim_dsm, dim_uem)

Computes Jacobian matrices for BL solution at an axial station. Called repeatedly by [`_axisymm_BL`](@ref). Formerly, `blsys!()`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `is_selfsimilar::Bool`: Self-similar BL profile flag.
      - `is_laminar::Bool`: Laminar flow flag.
      - `is_wake::Bool`: In wake? Flag.
      - `solves_direct::Bool`: Direct solution flag, with prescribed inviscid velocity ``u_e = u_\\mathrm{inv}``
      - `Mach::Float64`: Mach number for compressibility.
      - `uinv::Float64`: Inviscid velocity.
      - `x::Float64`: Arc length.
      - `b::Float64`: Lateral width of BL.
      - `rn::Float64`: ``dr/dn``, ``= 0`` for 2D.
      - `th::Float64`: Momentum thickness.
      - `ds::Float64`: Displacement thickness.
      - `ue::Float64`: Edge velocity.
      - `h::Float64`: Shape parameter.
      - `hk::Float64`: Kinematic shape parameter.
      - `hc::Float64`: density shape parameter (Whitfield).
      - `hs::Float64`: kinetic energy shape parameter.
      - `cf::Float64`: Skin friction factor.
      - `di::Float64`: Dissipation factor.

      `m` denotes the previous point (minus one) in the upstream.
      `_z` denotes partial derivative with respect to `z` (`z` = `th`, `ds`, `ue`).

      **Outputs:**
      - `aa::Array{Float64, 3, 3}`: Jacobian matrix (wrt current point vars).
      - `bb::Array{Float64, 3, 3}`: Jacobian matrix (wrt previous point vars).
      - `rr::Array{Float64, 3}`: Residual.

See Section 4 of [Simplified Viscous/Inviscid Analysis for Nearly-Axisymmetric Bodies](../assets/drela_TASOPT_2p16/axibl.pdf).
"""
function _BL_station_system(is_selfsimilar, is_laminar, is_wake, solves_direct, Mach, uinv,hksep,
                      x,b,rn,th,ds,ue,
                      h , h_th, h_ds,
                      hk, hk_th, hk_ds, hk_ue,
                      hc, hc_th, hc_ds, hc_ue,
                      hs, hs_th, hs_ds, hs_ue,
                      cf, cf_th, cf_ds, cf_ue,
                      di, di_th, di_ds, di_ue,
                      xm,bm,rnm,thm,dsm,uem, 
                      hm , hm_thm, hm_dsm,
                      hkm, hkm_thm, hkm_dsm, hkm_uem,
                      hcm, hcm_thm, hcm_dsm, hcm_uem,
                      hsm, hsm_thm, hsm_dsm, hsm_uem,
                      cfm, cfm_thm, cfm_dsm, cfm_uem,
                      dim, dim_thm, dim_dsm, dim_uem)
      
      
      aa = @MMatrix zeros(3,3)
      bb = @MMatrix zeros(3,3)
      rr = @MVector zeros(3)
     
      ɣ = 1.4
      gmi = ɣ - 1.0

      trat = 1.0 + 0.5*gmi*Mach^2 * (1.0-ue^2)
      trat_ue =       -gmi*Mach^2 *      ue

      amsq  = 1.0 + 0.5*gmi*Mach^2 * (1.0-uem^2)
      amsq_uem =       -gmi*Mach^2 *      uem

      rh = trat^(1.0/gmi)
      rh_ue = rh/(gmi*trat) * trat_ue

      rhm = amsq^(1.0/gmi)
      rhm_uem = rhm/(gmi*amsq) * amsq_uem

      if (is_selfsimilar)
       xl = 1.0
       bl = 1.0
       rl = 0.0
       ul = 1.0
       tl = 0.0
       hl = 0.0
       bl_ds = 0.
       rl_ue = 0.
       ul_ue = 0.
       tl_th = 0.
       
       hl_hs = 0.
       bl_dsm = 0.
       rl_uem = 0.
       ul_uem = 0.
       tl_thm = 0.
       hl_hsm = 0.

       hl_th = hl_hs*hs_th
       hl_ds = hl_hs*hs_ds
       hl_ue = hl_hs*hs_ue
       hl_thm = 0.
       hl_dsm = 0.
       hl_uem = 0.

       cfxa    = 0.5*cf   *x/th
       cfxa_th = 0.5*cf_th*x/th - cfxa/th
       cfxa_ds = 0.5*cf_ds*x/th
       cfxa_ue = 0.5*cf_ue*x/th
       cfxa_thm = 0.
       cfxa_dsm = 0.
       cfxa_uem = 0.

       dcxa    = di   *x/th - 0.5*cf   *x/th           
       dcxa_th = di_th*x/th - 0.5*cf_th*x/th - dcxa/th
       dcxa_ds = di_ds*x/th - 0.5*cf_ds*x/th
       dcxa_ue = di_ue*x/th - 0.5*cf_ue*x/th
       dcxa_thm = 0.
       dcxa_dsm = 0.
       dcxa_uem = 0.

       ha     = h
       ha_th  = h_th
       ha_ds  = h_ds
       ha_thm = 0.
       ha_dsm = 0.

       hsa     = hs
       hsa_th  = hs_th
       hsa_ds  = hs_ds
       hsa_ue  = hs_ue
       hsa_thm = 0.
       hsa_dsm = 0.
       hsa_uem = 0.

       hca     = hc
       hca_th  = hc_th
       hca_ds  = hc_ds
       hca_ue  = hc_ue
       hca_thm = 0.
       hca_dsm = 0.
       hca_uem = 0.

      else
       bd  = b  + 2.0*π*ds *rn
       bdm = bm + 2.0*π*dsm*rnm

       xl = log(x/xm)
       bl = log(bd/bdm)
       rl = log(rh/rhm)
       ul = log(ue/uem)
       tl = log(th/thm)
       hl = log(hs/hsm)
       bl_ds  =  1.0/bd * 2.0*π*rn
       rl_ue  =  1.0/rh * rh_ue
       ul_ue  =  1.0/ue
       tl_th  =  1.0/th
       hl_hs  =  1.0/hs
       bl_dsm = -1.0/bdm * 2.0*π*rnm
       rl_uem = -1.0/rhm * rhm_uem
       ul_uem = -1.0/uem
       tl_thm = -1.0/thm
       hl_hsm = -1.0/hsm

       hl_th = hl_hs*hs_th
       hl_ds = hl_hs*hs_ds
       hl_ue = hl_hs*hs_ue
       hl_thm = hl_hsm*hsm_thm
       hl_dsm = hl_hsm*hsm_dsm
       hl_uem = hl_hsm*hsm_uem

       cfx      = 0.5*cf     *x /th
       cfx_th   = 0.5*cf_th  *x /th  - cfx/th
       cfx_ds   = 0.5*cf_ds  *x /th
       cfx_ue   = 0.5*cf_ue  *x /th
       cfxm     = 0.5*cfm    *xm/thm
       cfxm_thm = 0.5*cfm_thm*xm/thm - cfxm/thm
       cfxm_dsm = 0.5*cfm_dsm*xm/thm
       cfxm_uem = 0.5*cfm_uem*xm/thm

       cfxa     = 0.5*(cfx + cfxm)
       cfxa_th  = 0.5* cfx_th
       cfxa_ds  = 0.5* cfx_ds
       cfxa_ue  = 0.5* cfx_ue
       cfxa_thm = 0.5*       cfxm_thm
       cfxa_dsm = 0.5*       cfxm_dsm
       cfxa_uem = 0.5*       cfxm_uem


       dcx      = di     *x /th  - 0.5*cf     *x /th           
       dcx_th   = di_th  *x /th  - 0.5*cf_th  *x /th - dcx/th
       dcx_ds   = di_ds  *x /th  - 0.5*cf_ds  *x /th
       dcx_ue   = di_ue  *x /th  - 0.5*cf_ue  *x /th
       dcxm     = dim    *xm/thm - 0.5*cfm    *xm/thm           
       dcxm_thm = dim_thm*xm/thm - 0.5*cfm_thm*xm/thm - dcxm/thm
       dcxm_dsm = dim_dsm*xm/thm - 0.5*cfm_dsm*xm/thm
       dcxm_uem = dim_uem*xm/thm - 0.5*cfm_uem*xm/thm
 
       dcxa     = 0.5*(dcx + dcxm)
       dcxa_th  = 0.5* dcx_th
       dcxa_ds  = 0.5* dcx_ds
       dcxa_ue  = 0.5* dcx_ue
       dcxa_thm = 0.5*       dcxm_thm
       dcxa_dsm = 0.5*       dcxm_dsm
       dcxa_uem = 0.5*       dcxm_uem

 
       ha     = 0.5*(h + hm)
       ha_th  = 0.5* h_th
       ha_ds  = 0.5* h_ds
       ha_thm = 0.5*     hm_thm
       ha_dsm = 0.5*     hm_dsm

       hsa     = 0.5*(hs + hsm)
       hsa_th  = 0.5* hs_th
       hsa_ds  = 0.5* hs_ds
       hsa_ue  = 0.5* hs_ue
       hsa_thm = 0.5*      hsm_thm
       hsa_dsm = 0.5*      hsm_dsm
       hsa_uem = 0.5*      hsm_uem

       hca     = 0.5*(hc + hcm)
       hca_th  = 0.5* hc_th
       hca_ds  = 0.5* hc_ds
       hca_ue  = 0.5* hc_ue
       hca_thm = 0.5*      hcm_thm
       hca_dsm = 0.5*      hcm_dsm
       hca_uem = 0.5*      hcm_uem

#       hca = 0.
#       hca_th  = 0.
#       hca_ds  = 0.
#       hca_ue  = 0.
#       hca_thm = 0.
#       hca_dsm = 0.
#       hca_uem = 0.

      end


      rr[1]   = tl     - cfxa    *xl + (ha + 2.0)*ul  +  bl + rl
      aa[1,1] = tl_th  - cfxa_th *xl +  ha_th    *ul
      aa[1,2] =        - cfxa_ds *xl +  ha_ds    *ul  +  bl_ds
      aa[1,3] =        - cfxa_ue *xl + (ha + 2.0)*ul_ue     + rl_ue
      bb[1,1] = tl_thm - cfxa_thm*xl +  ha_thm   *ul
      bb[1,2] =        - cfxa_dsm*xl +  ha_dsm   *ul  +  bl_dsm
      bb[1,3] =        - cfxa_uem*xl + (ha + 2.0)*ul_uem    + rl_uem


      btmp     =  2.0*hca/hsa + 1.0 - ha #This is the (2H**/H* + 1 - H) factor in the K.E. eqn
      btmp_hca =  2.0    /hsa
      btmp_hsa = -2.0*hca/hsa^2
      btmp_ha  =                    - 1.0

      btmp_th = btmp_hca*hca_th + btmp_hsa*hsa_th + btmp_ha*ha_th
      btmp_ds = btmp_hca*hca_ds + btmp_hsa*hsa_ds + btmp_ha*ha_ds
      btmp_ue = btmp_hca*hca_ue + btmp_hsa*hsa_ue

      btmp_thm = btmp_hca*hca_thm + btmp_hsa*hsa_thm + btmp_ha*ha_thm
      btmp_dsm = btmp_hca*hca_dsm + btmp_hsa*hsa_dsm + btmp_ha*ha_dsm
      btmp_uem = btmp_hca*hca_uem + btmp_hsa*hsa_uem

      rr[2]   = hl     - dcxa    *xl + btmp    *ul
      aa[2,1] = hl_th  - dcxa_th *xl + btmp_th *ul
      aa[2,2] = hl_ds  - dcxa_ds *xl + btmp_ds *ul
      aa[2,3] = hl_ue  - dcxa_ue *xl + btmp_ue *ul + btmp * ul_ue
      bb[2,1] = hl_thm - dcxa_thm*xl + btmp_thm*ul
      bb[2,2] = hl_dsm - dcxa_dsm*xl + btmp_dsm*ul
      bb[2,3] = hl_uem - dcxa_uem*xl + btmp_uem*ul + btmp * ul_uem

      if(solves_direct)
       rr[3]   = ue - uinv
       aa[3,1] = 0.
       aa[3,2] = 0.
       aa[3,3] = 1.0
      else
       rr[3]   = hk - hksep
       aa[3,1] = hk_th
       aa[3,2] = hk_ds
       aa[3,3] = hk_ue
      end

      bb[3,1] = 0.
      bb[3,2] = 0.
      bb[3,3] = 0.

      return aa, bb, rr

end # _BL_station_system



"""
    _BL_station_vars(is_selfsimilar, is_laminar, is_wake, Reyn, Mach, fexcr, x, θ, δs, ue)

Returns the boundary layer variables needed for solution. Formerly, `blvar!()`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `is_selfsimilar::Bool`: Self-similar BL profile flag.
      - `is_laminar::Bool`: Laminar flow flag.
      - `is_wake::Bool`: In wake flag.
      - `Reyn::Float64`: Reynolds number.
      - `Mach::Float64`: Mach number for compressibility.
      - `fexcr::Float64`: Excrescence factor.

      **Outputs:**
      - `h::Float64` : Shape parameter.
      - `hk::Float64`: Kinematic shape parameter.
      - `hc::Float64`: Density shape parameter (Whitfield).
      - `hs::Float64`: Kinetic energy shape parameter.
      - `cf::Float64`: Skin friction factor.
      - `cd::Float64`: Dissipation factor and their derivatives.

"""
function _BL_station_vars(is_selfsimilar,is_laminar,is_wake, Reyn,Mach, fexcr,
                      x, θ ,δs ,ue )

      acon = 6.0
      bcon = 0.72

      ɣ = 1.4
      gmi = ɣ - 1.0

      trat = 1.0 + 0.5*gmi*Mach^2 * (1.0-ue^2)
      trat_ue =       -gmi*Mach^2 *      ue    #deriv of trat wrt to ue d(trat)/d(ue)

      msq = (ue*Mach)^2 / trat
      msq_ue = 2.0*ue*Mach^2 - (msq/trat)*trat_ue

      h    =  δs/θ #H = delta star/ theta (2D shape parameter)
      ∂h_∂θ =  -h/θ
      ∂h_∂δs = 1.0/θ

      (hk, hk_h, hk_msq) = hkin( h , msq)
      ∂hk_∂θ = hk_h*∂h_∂θ
      ∂hk_∂δs = hk_h*∂h_∂δs
      hk_ue = hk_msq*msq_ue

      hk = max( hk , 1.005 )

      rh = trat^(1.0/gmi)
      rh_ue = rh/(gmi*trat) * trat_ue

#---- mu ~ T assumed here.  Could be replaced with Sutherland's Law
      mu = trat/Reyn
      mu_ue = trat_ue/Reyn

      rt    = rh*ue*θ/mu
      rt_ue = rh   *θ/mu + rh_ue*ue*θ/mu  - (rt/mu)*mu_ue
      rt_θ = rh*ue   /mu

      if(is_laminar)
       (hs, hs_hk, hs_rt, hs_msq) = hsl( hk, rt, msq)
      else
       (hs, hs_hk, hs_rt, hs_msq) = hst( hk, rt, msq)
      end
      ∂hs_∂θ = hs_hk*∂hk_∂θ + hs_rt*rt_θ
      ∂hs_∂δs = hs_hk*∂hk_∂δs
      hs_ue = hs_hk*hk_ue + hs_rt*rt_ue + hs_msq*msq_ue

#      (hc, hc_hk, hc_msq) = hct( hk, msq )
#      ∂hc_∂θ = hc_hk*∂hk_∂θ
#      ∂hc_∂δs = hc_hk*∂hk_∂δs
#      hc_ue = hc_hk*hk_ue + hc_msq*msq_ue

      hc    = 0.5*gmi*msq   *h
      hc_ue = 0.5*gmi*msq_ue*h
      ∂hc_∂θ = 0.5*gmi*msq   *∂h_∂θ
      ∂hc_∂δs = 0.5*gmi*msq   *∂h_∂δs

      if(is_wake) 
       cf    = 0.
       ∂cf_∂θ = 0.
       ∂cf_∂δs = 0.
       cf_ue = 0.
      else
       if(is_laminar)
        (cf, cf_hk, cf_rt, cf_msq) = cfl( hk, rt, msq )
       else
        (cf, cf_hk, cf_rt, cf_msq) = cft( hk, rt, msq )
       end
       cf     = fexcr * cf
       cf_hk  = fexcr * cf_hk
       cf_rt  = fexcr * cf_rt
       cf_msq = fexcr * cf_msq
       
       ∂cf_∂θ = cf_hk*∂hk_∂θ + cf_rt*rt_θ
       ∂cf_∂δs = cf_hk*∂hk_∂δs
       cf_ue = cf_hk*hk_ue + cf_rt*rt_ue + cf_msq*msq_ue
      end

      if(is_laminar) 
       (𝒟ᵢ, 𝒟ᵢ_hk, 𝒟ᵢ_rt) = 𝒟ᵢl( hk, rt )
       ∂𝒟ᵢ_∂θ = 𝒟ᵢ_hk*∂hk_∂θ + 𝒟ᵢ_rt*rt_θ
       ∂𝒟ᵢ_∂δs = 𝒟ᵢ_hk*∂hk_∂δs
       𝒟ᵢ_ue = 𝒟ᵢ_hk*hk_ue + 𝒟ᵢ_rt*rt_ue
      else
       hrat = (hk-1.0)/(acon*hk)
       hrat_hk = 1.0/(acon*hk^2)

       fc = sqrt(1.0 + 0.5*gmi*msq)
       fc_ue = (0.25*gmi/fc)*msq_ue

       uq      = (0.5*cf - hrat^2/fc  ) / (bcon*hk)
       uq_cf   =  0.5                    / (bcon*hk)
       uq_hrat =      -2.0*hrat   /fc    / (bcon*hk)
       uq_fc   =           hrat^2/fc^2 / (bcon*hk)

       uq_hk = uq_hrat*hrat_hk - uq/hk

       uq_θ = uq_cf*∂cf_∂θ + uq_hk*∂hk_∂θ
       uq_δs = uq_cf*∂cf_∂δs + uq_hk*∂hk_∂δs
       uq_ue = uq_cf*cf_ue + uq_hk*hk_ue + uq_fc*fc_ue

       𝒟ᵢ = 0.5*cf - (hk-1.0)*uq
       𝒟ᵢ_cf = 0.5
       𝒟ᵢ_hk = -uq
       𝒟ᵢ_uq = -(hk-1.0)

       ∂𝒟ᵢ_∂θ = 𝒟ᵢ_cf*∂cf_∂θ + 𝒟ᵢ_hk*∂hk_∂θ + 𝒟ᵢ_uq*uq_θ
       ∂𝒟ᵢ_∂δs = 𝒟ᵢ_cf*∂cf_∂δs + 𝒟ᵢ_hk*∂hk_∂δs + 𝒟ᵢ_uq*uq_δs
       𝒟ᵢ_ue = 𝒟ᵢ_cf*cf_ue + 𝒟ᵢ_hk*hk_ue + 𝒟ᵢ_uq*uq_ue
      end

      if(is_wake) 
       wfac = 2.0
       𝒟ᵢ    = wfac*𝒟ᵢ
       ∂𝒟ᵢ_∂θ = wfac*∂𝒟ᵢ_∂θ
       ∂𝒟ᵢ_∂δs = wfac*∂𝒟ᵢ_∂δs
       𝒟ᵢ_ue = wfac*𝒟ᵢ_ue
      end
 
      return h , ∂h_∂θ, ∂h_∂δs,
             hk, ∂hk_∂θ, ∂hk_∂δs, hk_ue,
             hc, ∂hc_∂θ, ∂hc_∂δs, hc_ue,
             hs, ∂hs_∂θ, ∂hs_∂δs, hs_ue,
             cf, ∂cf_∂θ, ∂cf_∂δs, cf_ue,
             𝒟ᵢ, ∂𝒟ᵢ_∂θ, ∂𝒟ᵢ_∂δs, 𝒟ᵢ_ue 
    
      end # _BL_station_vars



 
 
      function hkin( H, MSQ )
#
#---- calculate kinematic shape parameter (assuming air)
#     (from Whitfield )

      HK     =    (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
      HK_H   =     1.0          /(1.0 + 0.113*MSQ)
      HK_MSQ = (-.29 - 0.113*HK)/(1.0 + 0.113*MSQ)
#
      return HK, HK_H, HK_MSQ
      end
 


      function dil( HK, RT )
#
#---- Laminar dissipation function  ( 2 CD/H* )     (from Falkner-Skan)
      if(HK<4.0) 
       DI    = ( 0.00205  *  (4.0-HK)^5.5 + 0.207 ) / RT
       DI_HK = ( -.00205*5.5*(4.0-HK)^4.5         ) / RT
      else
       HKB = HK - 4.0
       DEN = 1.0 + 0.02*HKB^2
       DI    = ( -.0016  *  HKB^2  /DEN   + 0.207             ) / RT
       DI_HK = ( -.0016*2.0*HKB*(1.0/DEN - 0.02*HKB^2/DEN^2) ) / RT
      end
      DI_RT = -DI/RT
#
      return DI, DI_HK, DI_RT
      end


      function DILW( HK, RT )
#
      MSQ = 0.
      HS, HS_HK, HS_RT, HS_MSQ = hsl( HK, RT, MSQ)
#
#---- Laminar wake dissipation function  ( 2 CD/H* )
      RCD    =  1.10 * (1.0 - 1.0/HK)^2  / HK
      RCD_HK = -1.10 * (1.0 - 1.0/HK)*2.0 / HK^3-
	 RCD/HK
#
      DI    = 2.0*RCD   /(HS*RT)
      DI_HK = 2.0*RCD_HK/(HS*RT) - (DI/HS)*HS_HK
      DI_RT = -DI/RT             - (DI/HS)*HS_RT
#
      return DI, DI_HK, DI_RT
      end


      function hsl( HK, RT, MSQ)
#
#---- Laminar HS correlation
      if(HK<4.35) 
       TMP = HK - 4.35
       HS    = 0.0111*TMP^2/(HK+1.0) -
	 0.0278*TMP^3/(HK+1.0)  + 1.528 -
	 0.0002*(TMP*HK)^2
       HS_HK = 0.0111*(2.0*TMP    - TMP^2/(HK+1.0))/(HK+1.0) -
	 0.0278*(3.0*TMP^2 - TMP^3/(HK+1.0))/(HK+1.0) -
	 0.0002*2.0*TMP*HK * (TMP + HK)
      else
       HS    = 0.015*    (HK-4.35)^2/HK + 1.528
       HS_HK = 0.015*2.0*(HK-4.35)   /HK -
	 0.015*    (HK-4.35)^2/HK^2
      end
#
      HS_RT  = 0.
      HS_MSQ = 0.
#
      return HS, HS_HK, HS_RT, HS_MSQ
      end


      function cfl( HK, RT, MSQ )
#
#---- Laminar skin friction function  ( Cf )    ( from Falkner-Skan )
      if(HK < 5.5)
       TMP = (5.5-HK)^3 / (HK+1.0)
       CF    = ( 0.0727*TMP                      - 0.07       )/RT
       CF_HK = ( -.0727*TMP*3.0/(5.5-HK) - 0.0727*TMP/(HK+1.0))/RT
      else
       TMP = 1.0 - 1.0/(HK-4.5)
       CF    = ( 0.015*TMP^2      - 0.07  ) / RT
       CF_HK = ( 0.015*TMP*2.0/(HK-4.5)^2 ) / RT
      end
      CF_RT = -CF/RT
      CF_MSQ = 0.0
#
      return CF, CF_HK, CF_RT, CF_MSQ
      end



      function DIT( HS, US, CF, ST) 
#
#---- Turbulent dissipation function  ( 2 CD/H* )
      DI    =  ( 0.5*CF*US + ST*ST*(1.0-US) ) * 2.0/HS
      DI_HS = -( 0.5*CF*US + ST*ST*(1.0-US) ) * 2.0/HS^2
      DI_US =  ( 0.5*CF    - ST*ST          ) * 2.0/HS
      DI_CF =  ( 0.5   *US                  ) * 2.0/HS
      DI_ST =  (            2.0*ST*(1.0-US) ) * 2.0/HS
#
      return DI, DI_HS, DI_US, DI_CF, DI_ST
      end


      function hst( HK, RT, MSQ) 
#
#---- Turbulent HS correlation
#
      HSMIN = 1.50
      DHSINF= 0.015      
#
#---- ###  12/4/94
#---- limited Rtheta dependence for Rtheta < 200
#
#
      if(RT>400.0) 
       HO    = 3.0 + 400.0/RT
       HO_RT =     - 400.0/RT^2
      else
       HO    = 4.0
       HO_RT = 0.
      end
#
      if(RT>200.0) 
       RTZ    = RT
       RTZ_RT = 1.
      else
       RTZ    = 200.0
       RTZ_RT = 0.
      end
#
      if(HK<HO) 
#----- attached branch
# =======================================================
#----- old correlation
#-     (from Swafford profiles)
#       SRT = sqrt(RT)
#       HEX = (HO-HK)^1.6
#       RTMP = 0.165 - 1.6/SRT
#       HS    = HSMIN + 4.0/RT + RTMP*HEX/HK
#       HS_HK = RTMP*HEX/HK*(-1.6/(HO-HK) - 1.0/HK)
#       HS_RT = -4.0/RT^2 + HEX/HK*0.8/SRT/RT
#     &             + RTMP*HEX/HK*1.6/(HO-HK)*HO_RT
# =======================================================
#----- new correlation  29 Nov 91
#-     (from  arctan(y+) + Schlichting  profiles)
       HR    = ( HO - HK)/(HO-1.0)
       HR_HK =      - 1.0/(HO-1.0)
       HR_RT = (1.0 - HR)/(HO-1.0) * HO_RT
       HS    = (2.0-HSMIN-4.0/RTZ)*HR^2  * 1.5/(HK+0.5) + HSMIN + 4.0/RTZ
       HS_HK =-(2.0-HSMIN-4.0/RTZ)*HR^2  * 1.5/(HK+0.5)^2 +
	 (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_HK
       HS_RT = (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_RT +
	 (HR^2 * 1.5/(HK+0.5) - 1.0)*4.0/RTZ^2 * RTZ_RT
#
      else
#
#----- separated branch
       GRT = log(RTZ)
       HDif = HK - HO 
       RTMP = HK - HO + 4.0/GRT
       HTMP    = 0.007*GRT/RTMP^2 + DHSINF/HK
       HTMP_HK = -.014*GRT/RTMP^3 - DHSINF/HK^2
       HTMP_RT = -.014*GRT/RTMP^3 * (-HO_RT - 4.0/GRT^2/RTZ * RTZ_RT) + 0.007RTMP^2 / RTZ * RTZ_RT
       HS    = HDif^2 * HTMP + HSMIN + 4.0/RTZ
       HS_HK = HDif*2.0* HTMP + HDif^2 * HTMP_HK
       HS_RT = HDif^2 * HTMP_RT - 4.0/RTZ^2 * RTZ_RT + HDif*2.0* HTMP * (-HO_RT)
#
      end
#
#---- fudge HS slightly to make sure   HS -> 2   as   HK -> 1
#-    (unnecessary with new correlation)
#      HTF    = 0.485/9.0 * (HK-4.0)^2/HK  +  1.515
#      HTF_HK = 0.485/9.0 * (1.0-16.0/HK^2)
#      ARG = max( 10.0*(1.0 - HK) , -15.0 )
#      HXX = exp(ARG)
#      HXX_HK = -10.0*HXX
#C
#      HS_HK  = (1.0-HXX)*HS_HK  +  HXX*HTF_HK
#     &       + (        -HS     +      HTF    )*HXX_HK
#      HS_RT  = (1.0-HXX)*HS_RT
#      HS     = (1.0-HXX)*HS     +  HXX*HTF
#
#---- Whitfield's minor additional compressibility correction
      FM = 1.0 + 0.014*MSQ
      HS     = ( HS + 0.028*MSQ ) / FM
      HS_HK  = ( HS_HK          ) / FM
      HS_RT  = ( HS_RT          ) / FM
      HS_MSQ = 0.028/FM  -  0.014*HS/FM
#
      return HS, HS_HK, HS_RT, HS_MSQ
      end
 
 
"""
CFT returns the turbulent skin friction factor using the turbulent compressible version
of the Coles formula 


MSQ: (Mach number)^2
HK:  Kinematic shape parameter
RT:  Momentum thickness (theta) Reynolds number 

!!! details "🔃 Inputs and Outputs"
      **Inputs:**

      **Outputs:**
      - CF:   Skin friction factor
      - CF\\_HK: Derivative wrt to HK -> ``\\frac{\\partial C_f}{\\partial H_k}``
      - CF\\_RT: Derivative wrt to RT -> dCf
      - CF\\_MSQ: Derivative wrt to ``M^2``

"""
      function cft( HK, RT, MSQ )
      
      ɣ = 1.4
      gmi = ɣ - 1.0
      CFFAC = 1.0

#---- Turbulent skin friction function  ( Cf )    (Coles)
      
      FC = sqrt(1.0 + 0.5*gmi*MSQ)
      GRT = log(RT/FC)
      GRT = max(GRT,3.0)
 
      GEX = -1.74 - 0.31*HK
 
      ARG = -1.33*HK
      ARG = max(-20.0, ARG )

      THK = tanh(4.0 - HK/0.875)

#       CFO =  CFFAC * 0.3*exp(ARG) * (GRT/log(10))^GEX
      CFO =  CFFAC * 0.3*exp(ARG) * (GRT/2.3026)^GEX
      CF     = ( CFO  +  1.1E-4*(THK-1.0) ) / FC
      CF_HK  = (-1.33*CFO - 0.31*log(GRT/2.3026)*CFO -
                 1.1E-4*(1.0-THK^2) / 0.875    ) / FC
    
      CF_RT  = GEX*CFO/(FC*GRT) / RT
      CF_MSQ = GEX*CFO/(FC*GRT) * (-0.25*gmi/FC^2) - 0.25*gmi*CF/FC^2
#
      return CF, CF_HK, CF_RT, CF_MSQ
      end # CFT



 
      function hct( HK, MSQ) 
#
#---- density shape parameter    (from Whitfield)
      HC     = MSQ * (0.064/(HK-0.8) + 0.251)
      HC_HK  = MSQ * (-.064/(HK-0.8)^2     )
      HC_MSQ =        0.064/(HK-0.8) + 0.251
#
      return HC, HC_HK, HC_MSQ
      end

