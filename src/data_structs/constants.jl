
const ft_to_m = 0.3048
const in_to_m = 0.0254
const nmi_to_m = 1852.0
const deg_to_rad = π / 180.0
const lbf_to_N = 4.448222
const kts_to_mps = 0.51444
const hp_to_W = 745.7
const lb_N = 0.22480902473348890

const gee = 9.81
const gamSL = 1.4
const cpSL = 1004.0
const μAir = 1.78e-5
const pref = 101320.0
const Tref = 288.2
const σ_SB = 5.670374419e-8
const μ₀ = 1.25663706127e-6 #N⋅A⁻² https://physics.nist.gov/cgi-bin/cuu/Value?mu0 
const p_atm = 101325.0 #Pa in one atm
const t_h = 3600.0 #s in one hour

const seat_layouts = Dict{Int64, Vector{Int64}}(
    1 => [1, 0],
    2 => [1, 1],
    3 => [2, 1],
    4 => [2, 2],
    5 => [3, 2],
    6 => [3, 3],
    7 => [2, 3, 2],
    8 => [2, 4, 2],
    9 => [3, 3, 3],
    10 => [3, 4, 3],
    11 => [3, 5, 3],
    12 => [3, 3, 3, 3],
    13 => [3, 4, 3, 3],
    14 => [3, 4, 4, 3],
    15 => [3, 5, 4, 3],
    16 => [3, 5, 5, 3]
)

const aerodrome_codes = Dict{Int64, Vector{String}}(
    #max_wingspan => [ICAO code, FAA code]
    15 => ["A", "I"],
    24 => ["B", "II"],
    36 => ["C", "III"],
    52 => ["D", "IV"],
    65 => ["E", "V"],
    80 => ["F", "VI"],
)