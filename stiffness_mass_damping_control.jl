using LinearAlgebra, Plots, Printf

gr()  # Activate GR backend

#Times New Roman
default(fontfamily="times")

# Define parameter variations
k_values = [0.2, 1, 5]  # Varying stiffness
m_values = [0.2, 1, 5]  # Varying mass
c_values = [0.2, 1, 5]  # Varying damping

# Fixed values for other parameters
k_fixed = 1
m_fixed = 1
c_fixed = 1

colors = [:red, :blue, :green]  # Colors for three variations
titles = ["Stiffness Controlled", "Mass Controlled", "Damping Controlled"]  # Titles
filenames = ["Stiffness_Controlled.png", "Mass_Controlled.png", "Damping_Controlled.png"]

yy_lim = 3
ω = 0:0.01:5  # Frequency range

# Function to compute G(iω)
function compute_G(ω, k, m, c)
    ω_n = sqrt(k / m)
    c_cr = 2 * sqrt(k * m)
    ξ = c / c_cr

    r_values = ω ./ ω_n
    G_values = [sqrt(1 + (2 * ξ * r)^2) / sqrt((1 - r^2)^2 + (2 * ξ * r)^2) for r in r_values]

    return G_values, ω_n, ξ
end

# Function to find the maximum value of G(iω)
function find_max_G(G_values)
    return maximum(G_values)
end

# Function to calculate percentage change
function percentage_change(new_val, ref_val)
    return ((new_val - ref_val) / ref_val) * 100
end

# Function to highlight influential region
function plot_influential_region!(plt, min_ω_n, max_ω_n, mode)
    if mode == "k"
        x1, x2 = 0, max_ω_n * 0.9
    elseif mode == "m"
        x1, x2 = min_ω_n * 1.1, ω[end]
    elseif mode == "c"
        x1, x2 = 0.8 * max_ω_n, 1.2 * max_ω_n
    end
    
    plot!([x1, x2, x2, x1], [0, 0, yy_lim, yy_lim], 
          fillcolor=:yellow, fillalpha=0.25, seriestype=:shape, linecolor=:grey, lw=1, label="Influential region")
    vline!([max_ω_n], linestyle=:dash, color=:black, label="ω_n")
end

# Reference values for k=1, m=1, c=1
G_ref, ω_n_ref, ξ_ref = compute_G(ω, k_fixed, m_fixed, c_fixed)
G_max_ref = find_max_G(G_ref)

# Vary k while keeping m, c fixed
ω_n_values = [sqrt(k / m_fixed) for k in k_values]
min_ω_n, max_ω_n = minimum(ω_n_values), maximum(ω_n_values)

@printf("Stiffness Controlled\n \n")

plt1 = plot(xlabel="Frequency (ω)", ylabel="G(iω)", title=titles[1], legend=:topright, linewidth=2, ylims=(0, yy_lim), xlims=(0, ω[end]), framestyle=:box)
for (i, k) in enumerate(k_values)
    G_values, ω_n, ξ = compute_G(ω, k, m_fixed, c_fixed)
    G_max = find_max_G(G_values)

    ΔG_max = percentage_change(G_max, G_max_ref)
    Δω_n = percentage_change(ω_n, ω_n_ref)

    @printf("k = %.1f, m = %.1f, c = %.1f, ω_n = %.2f, ξ = %.4f, G_max = %.4f\n", k, m_fixed, c_fixed, ω_n, ξ, G_max)
    @printf("Percentage Change: Δω_n = %.2f%%, ΔG_max = %.2f%%\n\n", Δω_n, ΔG_max)

    plot!(ω, G_values, label="k=$k", color=colors[i])
    vline!([ω_n], linestyle=:dash, color=:black, label=false)
end
plot_influential_region!(plt1, min_ω_n, max_ω_n, "k")
savefig(plt1, filenames[1])
display(plt1)

# Vary m while keeping k, c fixed
ω_n_values = [sqrt(k_fixed / m) for m in m_values]
min_ω_n, max_ω_n = minimum(ω_n_values), maximum(ω_n_values)
@printf("Mass Controlled\n \n")
plt2 = plot(xlabel="Frequency (ω)", ylabel="G(iω)", title=titles[2], legend=:topright, linewidth=2, ylims=(0, yy_lim), xlims=(0, ω[end]), framestyle=:box)
for (i, m) in enumerate(m_values)
    G_values, ω_n, ξ = compute_G(ω, k_fixed, m, c_fixed)
    G_max = find_max_G(G_values)

    ΔG_max = percentage_change(G_max, G_max_ref)
    Δω_n = percentage_change(ω_n, ω_n_ref)

    @printf("k = %.1f, m = %.1f, c = %.1f, ω_n = %.2f, ξ = %.4f, G_max = %.4f\n", k_fixed, m, c_fixed, ω_n, ξ, G_max)
    @printf("Percentage Change: Δω_n = %.2f%%, ΔG_max = %.2f%%\n\n", Δω_n, ΔG_max)

    plot!(ω, G_values, label="m=$m", color=colors[i])
    vline!([ω_n], linestyle=:dash, color=:black, label=false)
end
plot_influential_region!(plt2, min_ω_n, max_ω_n, "m")
savefig(plt2, filenames[2])
display(plt2)

# Vary c while keeping k, m fixed
ω_n_fixed = sqrt(k_fixed / m_fixed)
@printf("Damping Controlled\n \n")
plt3 = plot(xlabel="Frequency (ω)", ylabel="G(iω)", title=titles[3], legend=:topright, linewidth=2, ylims=(0, yy_lim), xlims=(0, ω[end]), framestyle=:box)
for (i, c) in enumerate(c_values)
    G_values, ω_n, ξ = compute_G(ω, k_fixed, m_fixed, c)
    G_max = find_max_G(G_values)

    ΔG_max = percentage_change(G_max, G_max_ref)
    Δω_n = percentage_change(ω_n, ω_n_ref)

    @printf("k = %.1f, m = %.1f, c = %.1f, ω_n = %.2f, ξ = %.4f, G_max = %.4f\n", k_fixed, m_fixed, c, ω_n, ξ, G_max)
    @printf("Percentage Change: Δω_n = %.2f%%, ΔG_max = %.2f%%\n\n", Δω_n, ΔG_max)

    plot!(ω, G_values, label="c=$c", color=colors[i])
    vline!([ω_n], linestyle=:dash, color=:black, label=false)
end
plot_influential_region!(plt3, ω_n_fixed, ω_n_fixed, "c")
savefig(plt3, filenames[3])
display(plt3)
