# Fuselage layout
@kwdef mutable struct FuselageLayout
    fuse_radius::Float64 = 1.9558 # = ac.parg[igRfuse] #Rfuse
    bubble_lower_downward_shift::Float64 = 0.381#dRfuse 
    bubble_center_y_offset::Float64 = 0 #wfb
    n_webs::Float64 = 1 #nfwebs
    h_floor::Float64 = 0.127# = ac.parg[ighfloor] #hfloor
    x_nose::Float64  = 0# = ac.parg[igxnose] #xnose
    x_pressure_shell_fwd::Float64 = 5.1816# = ac.parg[igxshell1] #xshell1
    x_pressure_shell_aft::Float64 = 31.0896# = ac.parg[igxshell2] #xshell2
    x_start_cylinder::Float64 = 6.096# = ac.parg[igxblend1] #xblend1
    x_end_cylinder::Float64 = 29.5656# = ac.parg[igxblend2] #xblend2
    x_cone_end::Float64 = 35.6616# = ac.parg[igxend] #xconeend
    tailcone_taper_ratio::Float64 = 0.3# lambdac
    floor_height::Float64 = 0.127
    nose_radius::Float64 = 1.65
    tail_radius::Float64 = 2.0
    taper_fuse::Int64 = 1 # 0 = point ; 1 = edge
end

# function FuselageLayout(;default = true)
#     #TODO add read input
#     # if default

#     # else

#     # end
#     return FuselageLayout(fuse_radius = 3, )
# end