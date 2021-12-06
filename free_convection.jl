using Oceananigans
using Oceananigans.Units
using Oceananigans: fields
using Printf
using GLMakie
using JLD2

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

    prefix = "free_convection_$N"

    simulation.output_writers[:fields] = JLD2OutputWriter(model, fields(model),
                                                          prefix = prefix * "_fields",
                                                          schedule = TimeInterval(stop_time/10),
                                                          field_slicer = nothing,
                                                          force = true)
        
    u, v, w, b = fields(model)
    κₑ = model.diffusivity_fields.κₑ.b
    wb = AveragedField(w * b, dims=(1, 2))
    qᵇ = AveragedField(- ∂z(b) * κₑ, dims=(1, 2))
    simulation.output_writers[:averages] = JLD2OutputWriter(model, (; wb, qᵇ),
                                                            prefix = prefix * "_averages",
                                                            schedule = TimeInterval(stop_time/10),
                                                            init = file -> (file["surface_flux"] = Qᵇ),
                                                            force = true)

    run!(simulation)

    return prefix
end

# Run the free convection experiment
# prefix = run_free_convection()

# Generate file paths from the file prefix
fields_path = prefix * "_fields.jld2"
averages_path = prefix * "_averages.jld2"

# Load vertical velocity data and coordinates
w = FieldTimeSeries(fields_path, "w")
x, y, z = nodes(w)
Nt = size(w, 4)
Nz = size(w, 3) - 1
wmax = maximum(abs, w[Nt])
n = Node(1) # `n` represents the "snapshot index" (n varies from 1 to 11 here)

# Build the figure
fig = Figure(resolution=(1600, 800))

# A heatmap of vertical velocity
# Note: @lift interprets the node `n` so wⁿ can be updated dynamically to produce an animation
w_title = @lift "Vertical velocity at t = " * prettytime(w.times[$n])
wⁿ = @lift Array(interior(w[$n]))[1, :, :]

ax = Axis(fig[1, 1:3]; title=w_title, xlabel="x (m)", ylabel="z (m)")
hm = heatmap!(ax, y, z, wⁿ, limits=(-wmax, wmax), colormap=:balance) 
cb = Colorbar(fig[1, 4], limits=(-wmax, wmax), colormap=:balance)

# Line plots of the vertical fluxes
fluxes_label = @lift "xy-averaged fluxes at t = " * prettytime(w.times[$n])
averages_file = jldopen(averages_path)
Qᵇ = averages_file["surface_flux"]
iterations = parse.(Int, keys(averages_file["timeseries/t"]))
wb = @lift averages_file["timeseries/wb/" * string(iterations[$n])][:]
qᵇ = @lift vcat(averages_file["timeseries/qᵇ/" * string(iterations[$n])][1:Nz], [Qᵇ])

ax = Axis(fig[1, 5]; xlabel=fluxes_label, ylabel="z (m)")
lines!(ax, wb, z, label="Resolved, ⟨wb⟩")
lines!(ax, qᵇ, z, label="Unresolved, -⟨κₑ ∂z(b)⟩")
axislegend(ax, position=:rb)
xlims!(ax, -1e-8, 1e-8)

# Update figure data to produce an animation
record(fig, prefix * ".mp4", 1:Nt, framerate=8) do nn
    @info "Animating frame $nn of $Nt..."
    n[] = nn
end

# Close the file with horizontal averages
close(averages_file)
