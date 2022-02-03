using Oceananigans
using Oceananigans.Units
using Printf

L = 10
grid = RectilinearGrid(size=(64, 64, 128), x=(-L/2, L/2), y=(-L/2, L/2), z=(-L/2, L/2),
                       topology=(Periodic, Periodic, Bounded))


# Our basic state thus has a thin layer of stratification in the center of
# the channel, embedded within a thicker shear layer surrounded by unstratified fluid.

Ri, h = 0.1, 1/4

IsotropicDiffusivity(ν=2e-4, κ=2e-4)
closure = SmagorinskyLilly(C=16)
model = NonhydrostaticModel(timestepper = :RungeKutta3,
                              advection = UpwindBiasedFifthOrder(),
                                   grid = grid,
                               coriolis = nothing,
                                closure = closure,
                               buoyancy = BuoyancyTracer(),
                                tracers = :b)

amplitude = 1e-4
noise(x, y, z) = amplitude * randn()

shear_flow(x, y, z) = tanh(z) + noise(x, y, z)
stratification(x, y, z) = h * Ri * tanh(z / h) + noise(x, y, z)
set!(model, u=shear_flow, v=noise, w=noise, b=stratification)

stop_time = 200
Δt_diffusive = grid.Δzᵃᵃᶜ^2 / maximum(model.diffusivity_fields.νₑ)
Δt_advective = grid.Δzᵃᵃᶜ / maximum(model.velocities.u)
simulation = Simulation(model, Δt=0.5 * min(Δt_advective, Δt_diffusive),
                        stop_time=stop_time)


# Simple logging
using Oceanostics: SingleLineProgressMessenger
progress = SingleLineProgressMessenger(LES=true)
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

# Adaptive time-stepping
wizard = TimeStepWizard(cfl=0.8, diffusive_cfl=0.5, max_change=1.02, min_change=0.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

u, v, w = model.velocities
b = model.tracers.b

total_vorticity = Field(∂z(u) - ∂x(w))

simulation.output_writers[:snapshots] =
    JLD2OutputWriter(model, (Ω=total_vorticity, b=b, νₑ=model.diffusivity_fields.νₑ),
                     schedule = TimeInterval(1),
                     field_slicer = FieldSlicer(j=1),
                     prefix = "kelvin_helmholtz_instability",
                     force = true)

@info "*** Running a simulation of Kelvin-Helmholtz instability..."
run!(simulation)



# Now we plot stuff
using Plots
using Printf
using JLD2

file = jldopen(simulation.output_writers[:snapshots].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

@info "Making a neat movie of stratified shear flow..."

xF, yF, zF = nodes(total_vorticity)
xC, yC, zC = nodes(b)

function eigenplot(ω, b, σ, t; ω_lim=maximum(abs, ω)+1e-16, b_lim=maximum(abs, b)+1e-16)

    kwargs = (xlabel="x", ylabel="z", linewidth=0, label=nothing, color = :balance, aspectratio = 1,)

    ω_title(t) = t == nothing ? @sprintf("vorticity") : @sprintf("vorticity at t = %.2f", t)

    plot_ω = heatmap(xF, zF, clamp.(ω, -ω_lim, ω_lim)';
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
    b_snapshot = file["timeseries/b/$iteration"][:, 1, :]

    eigenmode_plot = eigenplot(ω_snapshot, b_snapshot, nothing, t; ω_lim=1, b_lim=0.05)

    plot(eigenmode_plot, size=(800, 600))
end

mp4(anim_total, "kelvin_helmholtz_instability_total.mp4", fps = 8) # hide
