using LinearAlgebra, Plots, Printf

gr()  
default(fontfamily="times")

# Initial conditions
x_0 = 0 
v_0 = 2*π 

# System parameters
c = 0.2   
k = 1     
m = 2     

# Derived parameters
ω_n = sqrt(k/m)
c_cr = 2 * m * ω_n
ξ = c / c_cr  
ω_d = ω_n * sqrt(1 - ξ^2)

# Constants for the solution
c1 = ((ξ * ω_n * x_0) + v_0) / ω_d
c2 = x_0
Φ = atan(c1/c2)
C = sqrt(c1^2 + c2^2)

# Time values
dt = 0.1  # Time step
t_values = 0:dt:100

# Displacement and velocity
x = [C * exp(-ξ * ω_n * t) * cos(ω_d * t - Φ) for t in t_values]
v = [C * (-ξ * ω_n) * exp(-ξ * ω_n * t) * cos(ω_d * t - Φ) - C * exp(-ξ * ω_n * t) * ω_d * sin(ω_d * t - Φ) for t in t_values]

# Energy calculations
PE = 0.5 * k .* x .^ 2
KE = 0.5 * m .* v .^ 2
TE = PE .+ KE  # Total Energy

# Proper Dissipated Energy using integral of cv²
DE = cumsum(c .* v.^2) .* (t_values[2] - t_values[1])  # Numerical integration using summation

# Plot displacement
p1 = plot(t_values, x, label="x(t)", xlabel="Time (s)", ylabel="Displacement", linewidth=2, title="Displacement vs Time", color=:black)

# Plot potential energy
p2 = plot(t_values, PE, label="Potential Energy", xlabel="Time (s)", ylabel="Energy", linewidth=2, title="Potential Energy vs Time", color=:blue)

# Plot kinetic energy
p3 = plot(t_values, KE, label="Kinetic Energy", xlabel="Time (s)", ylabel="Energy", linewidth=2, title="Kinetic Energy vs Time", color=:red)

# Plot total energy
p4 = plot(t_values, TE, label="Total Energy", xlabel="Time (s)", ylabel="Energy", linewidth=2, title="Total Energy vs Time", linestyle=:dash, color=:green)

# Combined energy plot
p5 = plot(t_values, PE, label="Potential Energy", linewidth=2)
plot!(t_values, KE, label="Kinetic Energy", linewidth=2)
plot!(t_values, TE, label="Total Energy", linestyle=:dash, linewidth=2, title="Energy Components Over Time")

# Plot for dissipated energy (computed properly using integral of c*v²)
p6 = plot(t_values, DE, label="Dissipated Energy", xlabel="Time (s)", ylabel="Energy Dissipated", linewidth=2, title="Dissipated Energy Over Time", linestyle=:dashdot, color=:red)
@printf "Natural Frequency (ω_n): %.2f\n" ω_n
@printf "Damping Ratio (ξ): %.2f\n" ξ
@printf "Damped Frequency (ω_d): %.3f\n" ω_d


# Display all plots
display(p1)
display(p2)
display(p3)
display(p4)
display(p5)
display(p6)


savefig(p1, "2a_Displacement_vs_Time.png")
savefig(p2, "2a_Potential_Energy_vs_Time.png")
savefig(p3, "2a_Kinetic_Energy_vs_Time.png")
savefig(p4, "2a_Total_Energy_vs_Time.png")
savefig(p5, "2a_Energy_Components_Over_Time.png")
savefig(p6, "2a_Dissipated_Energy_vs_Time.png")

