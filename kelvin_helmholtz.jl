using Oceananigans

grid = RectilinearGrid(size=(64, 64), x=(-5, 5), z=(-5, 5),
                       topology=(Periodic, Flat, Bounded))

shear_flow(x, y, z, t) = tanh(z)

stratification(x, y, z, t, p) = p.h * p.Ri * tanh(z / p.h)

U = BackgroundField(shear_flow)

B = BackgroundField(stratification, parameters=(Ri=0.1, h=1/4))

# Our basic state thus has a thin layer of stratification in the center of
# the channel, embedded within a thicker shear layer surrounded by unstratified fluid.

using Plots

zF = znodes(Face, grid)
zC = znodes(Center, grid)

Ri, h = B.parameters


model = NonhydrostaticModel(timestepper = :RungeKutta3,
                              advection = UpwindBiasedFifthOrder(),
                                   grid = grid,
                               coriolis = nothing,
                      background_fields = (u=U, b=B),
                                closure = IsotropicDiffusivity(ν=2e-4, κ=2e-4),
                               buoyancy = BuoyancyTracer(),
                                tracers = :b)

amplitude = 1e-2
noise(x, y, z) = amplitude * randn()
set!(model, u=noise, w=noise)

simulation = Simulation(model, Δt=0.1, stop_time=200)

u, v, w = model.velocities
b = model.tracers.b

total_vorticity = ComputedField(∂z(u) + ∂z(model.background_fields.velocities.u) - ∂x(w))

total_b = ComputedField(b + model.background_fields.tracers.b)

simulation.output_writers[:vorticity] =
    JLD2OutputWriter(model, (Ω=total_vorticity, b=b, B=total_b),
                     schedule = TimeInterval(1),
                     prefix = "kelvin_helmholtz_instability",
                     force = true)

# And now we...

@info "*** Running a simulation of Kelvin-Helmholtz instability..."
run!(simulation)

# ## Pretty things
#
# Load it; plot it. First the nonlinear equilibration of the perturbation fields together with the evolution of the
# kinetic energy.

using JLD2

file = jldopen(simulation.output_writers[:vorticity].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

@info "Making a neat movie of stratified shear flow..."

time = []

xF, yF, zF = nodes(total_vorticity)

xC, yC, zC = nodes(b)


using Printf
function eigenplot(ω, b, σ, t; ω_lim=maximum(abs, ω)+1e-16, b_lim=maximum(abs, b)+1e-16)

    kwargs = (xlabel="x", ylabel="z", linewidth=0, label=nothing, color = :balance, aspectratio = 1,)

    ω_title(t) = t == nothing ? @sprintf("vorticity") : @sprintf("vorticity at t = %.2f", t)

    plot_ω = contourf(xF, zF, clamp.(ω, -ω_lim, ω_lim)';
                      levels = range(-ω_lim, stop=ω_lim, length=20),
                       xlims = (xF[1], xF[grid.Nx]),
                       ylims = (zF[1], zF[grid.Nz]),
                       clims = (-ω_lim, ω_lim),
                       title = ω_title(t), kwargs...)

    b_title(t) = t == nothing ? @sprintf("buoyancy") : @sprintf("buoyancy at t = %.2f", t)

    plot_b = contourf(xC, zC, clamp.(b, -b_lim, b_lim)';
                    levels = range(-b_lim, stop=b_lim, length=20),
                     xlims = (xC[1], xC[grid.Nx]),
                     ylims = (zC[1], zC[grid.Nz]),
                     clims = (-b_lim, b_lim),
                     title = b_title(t), kwargs...)

    return plot(plot_ω, plot_b, layout=(1, 2), size=(800, 380))
end

anim_total = @animate for (i, iteration) in enumerate(iterations)

    @info "Plotting frame $i from iteration $iteration..."

    t = file["timeseries/t/$iteration"]
    ω_snapshot = file["timeseries/Ω/$iteration"][:, 1, :]
    b_snapshot = file["timeseries/B/$iteration"][:, 1, :]

    push!(time, t)

    eigenmode_plot = eigenplot(ω_snapshot, b_snapshot, nothing, t; ω_lim=1, b_lim=0.05)

    plot(eigenmode_plot, size=(800, 600))
end

mp4(anim_total, "kelvin_helmholtz_instability_total.mp4", fps = 8) # hide
