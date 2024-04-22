"""
    RK4(dy_dx, x, y0, u, p) 

This function uses a 4th-order Runge-Kutta method to integrate a vector ODE in space or time.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `dy_dx::Function`: function of the form dy_dx(x, y, u, p) for derivative
    - `x::Vector{Float64}`: vector with the spatial coordinates for integration
    - `y0::Vector{Float64}`: vector with the initial conditions
    - `u`: object with inputs 
    - `p`: object with parameters 
    
    **Outputs:**
    - `y::Matrix{Float64}`: matrix with values of y for every point in x 
    - `yend::Vector{Float64}`: vector with final conditions at x[end]
"""
function RK4(dy_dx, x, y0, u, p)
    y = zeros(Float64, length(y0), length(x))
    y[:, 1] .= y0
    for i = 1:(length(x) - 1)
        Δx = x[i+1] - x[i]
        k_1 = dy_dx(x[i], y[:, i], u, p)
        k_2 = dy_dx(x[i] + Δx/2, y[:, i] + Δx/2 * k_1, u, p)
        k_3 = dy_dx(x[i] + Δx/2, y[:, i] + Δx/2 * k_2, u, p)
        k_4 = dy_dx(x[i] + Δx, y[:, i] + Δx * k_3, u, p)

        y[:, i+1] = y[:, i] + 1/6 * Δx * (k_1 + 2*k_2 + 2*k_3 + k_4)
    end
    yend = y[:,end]
    return y, yend
end