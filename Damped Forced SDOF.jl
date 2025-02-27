using LinearAlgebra, Plots, Printf

gr()  
default(fontfamily="times")

# Initial conditions
x_0 = 0 
v_0 = 2*π 

# System parameters
c = 0.2    # Damping coefficient
k = 1      # Stiffness
m = 2      # Mass

# Forcing function parameters
F = 0.5    # Force amplitude
ω_f = 0.2  # Forcing frequency

# Derived parameters
ω_n = sqrt(k/m)  # Natural frequency
r = ω_f / ω_n    # Frequency ratio

c_cr = 2 * m * ω_n
ξ = c / c_cr      # Damping ratio
ω_d = ω_n * sqrt(1 - ξ^2)

# Constants for the homogeneous solution
c1 = ((ξ * ω_n * x_0) + v_0) / ω_d
c2 = x_0
Φ = atan(c1/c2)
C = sqrt(c1^2 + c2^2)

# Particular solution constants
denom = ((1 - r^2)^2 + (2 * ξ * r)^2)
H = 1 / denom
Xp = (F/m) * sqrt(H)
φ_p = atan((2 * ξ * r) / (1 - r^2))

@printf "Natural Frequency (ω_n): %.2f\n" ω_n
@printf "Damping Ratio (ξ): %.2f\n" ξ
@printf "Damped Frequency (ω_d): %.3f\n" ω_d
@printf("r = %.4f\n", r)
@printf("H = %.6f\n", H)

# Time values
dt = 0.05 
t_values = 0:dt:300

# Force function
force = [F * cos(ω_f * t) for t in t_values]

# Displacement and velocity
x_transient = [C * exp(-ξ * ω_n * t) * cos(ω_d * t - Φ) for t in t_values]  # Transient response
x_steady = [Xp * sin(ω_f * t - φ_p) for t in t_values]  # Steady-state response
x_total = x_transient .+ x_steady  # Total displacement

# Velocity calculation
v = [C * (-ξ * ω_n) * exp(-ξ * ω_n * t) * cos(ω_d * t - Φ) - 
     C * exp(-ξ * ω_n * t) * ω_d * sin(ω_d * t - Φ) + 
     Xp * ω_f * cos(ω_f * t - φ_p) for t in t_values]

# Energy calculations
PE = 0.5 * k .* x_total .^ 2
KE = 0.5 * m .* v .^ 2
TE = PE .+ KE  # Total Energy
DE = cumsum(c .* v.^2) .* dt  # Numerical integration using summation

# Compute external input energy
E_input = cumsum(force .* v) .* dt  # Work done by the force over time

# Plot force
p_force = plot(t_values, force, label="Applied Force", color=:purple, linewidth=2)
title!("Applied Force vs Time")
xlabel!("Time (s)")
ylabel!("Force")
savefig("2b_force_vs_time.png")

# Plot displacement (Transient and Steady-State response separately)
p1 = plot(t_values, x_transient, label="Transient Response", color=:blue, linewidth=2, linestyle=:dash)
plot!(t_values, x_steady, label="Steady-State Response", color=:green, linewidth=2, linestyle=:dash)
plot!(t_values, x_total, label="Total Displacement", color=:red, linewidth=1)
title!("Displacement vs Time")
xlabel!("Time (s)")
ylabel!("Displacement")
savefig("2b_displacement_vs_time.png")

# Plot potential energy
p2 = plot(t_values, PE, label="Potential Energy", xlabel="Time (s)", ylabel="Energy", linewidth=2, title="Potential Energy vs Time", color=:blue)
savefig("2b_potential_energy.png")

# Plot kinetic energy
p3 = plot(t_values, KE, label="Kinetic Energy", xlabel="Time (s)", ylabel="Energy", linewidth=2, title="Kinetic Energy vs Time", color=:red)
savefig("2b_kinetic_energy.png")

# Plot total energy
p4 = plot(t_values, TE, label="Total Energy", xlabel="Time (s)", ylabel="Energy", linewidth=2, title="Total Energy vs Time", linestyle=:dash, color=:green)
savefig("2b_total_energy.png")

# Combined energy plot
p5 = plot(t_values, PE, label="Potential Energy", linewidth=2)
plot!(t_values, KE, label="Kinetic Energy", linewidth=2)
plot!(t_values, TE, label="Total Energy", linestyle=:dash, linewidth=2, title="Energy Components Over Time")
savefig("2b_energy_components.png")

# Plot for dissipated energy
p6 = plot(t_values, DE, label="Dissipated Energy", xlabel="Time (s)", ylabel="Energy Dissipated", linewidth=2, title="Dissipated Energy Over Time", linestyle=:dashdot, color=:red)
savefig("2b_dissipated_energy.png")

# Plot for external input energy
p7 = plot(t_values, E_input, label="External Input Energy", xlabel="Time (s)", ylabel="Energy Input", linewidth=2, title="External Input Energy Over Time", color=:black)
savefig("2b_external_input_energy.png")

# Display all plots
display(p_force)
display(p1)
display(p2)
display(p3)
display(p4)
display(p5)
display(p6)
display(p7)
