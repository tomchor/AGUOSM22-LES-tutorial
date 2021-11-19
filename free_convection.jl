using Oceananigans
using Oceananigans.Units
using Oceananigans: fields
using Printf
using GLMakie

function run_free_convection(; N=32, H=64, Qᵇ=1e-8, N²=1e-6, advection=UpwindBiasedFifthOrder(), model_kwargs...)

    grid = RectilinearGrid(size=(N, N, N), x=(0, 2H), y=(0, 2H), z=(-H, 0), halo=(3, 3, 3))

    buoyancy_boundary_conditions =
        FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵇ), bottom=GradientBoundaryCondition(N²))                               
    model = NonhydrostaticModel(; grid = grid,
                                  tracers = :b,
                                  advection = advection,
                                  timestepper = :RungeKutta3,
                                  buoyancy = BuoyancyTracer(),
                                  boundary_conditions = (; b=buoyancy_boundary_conditions),
                                  closure = AnisotropicMinimumDissipation(),
                                  model_kwargs...)

    # Initial condition: linear stratification + surface-concentrated random noise
    bᵢ(x, y, z) = N² * z + 1e-3 * N² * H * exp(16z / H) * randn()
    set!(model, b=bᵢ)

    stop_time = (H/2)^2 * N² / 3Qᵇ # run till boundary layer depth ≈ H/2
    simulation = Simulation(model; Δt=10.0, stop_time)

    # Adaptive time-stepping
    simulation.callbacks[:wizard] = Callback(TimeStepWizard(cfl=0.5, max_Δt=1minute), IterationInterval(10))

    # Simple logging
    progress(sim) = @info @sprintf("[%.2f %%] iter %d, t = %s, Δt = %s, max(|w|): %.2e",
                                   100time(sim) / stop_time,
                                   iteration(sim),
                                   prettytime(sim),
                                   prettytime(sim.Δt),
                                   maximum(abs, sim.model.velocities.w))

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    simulation.output_writers[:fields] = JLD2OutputWriter(model, fields(model),
                                                          prefix = "free_convection_$N",
                                                          schedule = TimeInterval(stop_time/10),
                                                          field_slicer = nothing,
                                                          force = true)
        
    κₑ = model.diffusivity_fields.κₑ.b
    wb = AveragedField(w * b, dims=(1, 2))
    qᵇ = AveragedField(κₑ * b, dims=(1, 2))
    simulation.

    run!(simulation)

    return simulation.output_writers[:fields].filepath
end

path = run_free_convection()

# Visualize
w = FieldTimeSeries(path, "w")
n = Node(1)
wⁿ = @lift Array(interior(w[$n]))[1, :, :]
title = @lift "Vertical velocity at t = " * prettytime(w.times[$n])

fig = Figure(resolution=(800, 600))
ax = Axis(fig[1, 1]; title)
heatmap!(ax, wⁿ)

record(fig, basename(path)[1:end-5] * ".mp4", 1:length(w.times), framerate=8) do nn
    @info "Animating frame $nn of $(length(w.times))..."
    n[] = nn
end
