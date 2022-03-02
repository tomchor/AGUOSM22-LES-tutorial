using Oceananigans
using Oceananigans.Units
using Plots
using Printf
using JLD2

L = 10
grid = RectilinearGrid(size=(32, 32, 64), x=(-L/2, L/2), y=(-L/2, L/2), z=(-L/2, L/2),
                       halo=(3,3,3,),
                       topology=(Periodic, Periodic, Bounded))

#
# Our basic state thus has a thin layer of stratification in the center of
# the channel, embedded within a thicker shear layer surrounded by unstratified fluid.
function run_kelvin_helmholtz(closure; grid=grid, 
                              Ri=0.1, h=1/4, amplitude=1e-3,
                              stop_time=200,
                              advection=UpwindBiasedFifthOrder(), 
                              model_kwargs...,
                              )

    model = NonhydrostaticModel(timestepper = :RungeKutta3,
                                  advection = advection,
                                       grid = grid,
                                   coriolis = nothing,
                                    closure = closure,
                                   buoyancy = BuoyancyTracer(),
                                    tracers = :b)
    @show model

    noise(x, y, z) = amplitude * randn()

    shear_flow(x, y, z) = tanh(z) + noise(x, y, z)
    stratification(x, y, z) = h * Ri * tanh(z / h) + noise(x, y, z)
    set!(model, u=shear_flow, v=noise, w=noise, b=stratification)

    Δt_diffusive = grid.Δzᵃᵃᶜ^2 / maximum(model.diffusivity_fields.νₑ)
    Δt_advective = grid.Δzᵃᵃᶜ / maximum(model.velocities.u)
    simulation = Simulation(model, Δt=0.5 * min(Δt_advective, Δt_diffusive),
                            stop_time=stop_time)


    # Simple logging
    start_time = 1e-9*time_ns()
    # Simple logging
    progress(sim) = @info @sprintf("[%.2f %%] iter %d, t = %s, Δt = %s, max(|w|): %.2e",
                                   100time(sim) / stop_time,
                                   iteration(sim),
                                   prettytime(sim),
                                   prettytime(sim.Δt),
                                   maximum(abs, sim.model.velocities.w))
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    # Adaptive time-stepping
    wizard = TimeStepWizard(cfl=0.8, diffusive_cfl=0.5, max_change=1.05, min_change=0.1)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

    u, v, w = model.velocities
    b = model.tracers.b

    total_vorticity = Field(∂z(u) - ∂x(w))

    closure_name = string(typeof(closure).name.wrapper) # Get closure's name
    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, (ω=total_vorticity, b=b, νₑ=model.diffusivity_fields.νₑ),
                                                             schedule = TimeInterval(1),
                                                             field_slicer = FieldSlicer(j=1),
                                                             prefix = "kelvin_helmholtz_instability_$closure_name",
                                                             force = true)

    @info "*** Running a simulation of Kelvin-Helmholtz instability with $closure_name closure"
    run!(simulation)

    return simulation

end 


function plot_video(simulation; fps=14, size=(800, 600))
    jld2writer = simulation.output_writers[:snapshots]
    file = jldopen(jld2writer.filepath)

    iterations = parse.(Int, keys(file["timeseries/t"]))

    @info "Making a neat movie from $(jld2writer.filepath)"

    xF, yF, zF = nodes(jld2writer.outputs[:ω])
    xC, yC, zC = nodes(jld2writer.outputs[:b])

    function eigenplot(ω, b, σ, t; ω_lim=maximum(abs, ω)+1e-16, b_lim=maximum(abs, b)+1e-16)

        kwargs = (xlabel="x", ylabel="z", linewidth=0, label=nothing, color = :balance,)

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

        return plot(plot_ω, plot_b, layout=(1, 2),)
    end

    anim_total = @animate for (i, iteration) in enumerate(iterations)

        @info "Plotting frame $i from iteration $iteration..."

        t = file["timeseries/t/$iteration"]
        ω_snapshot = file["timeseries/ω/$iteration"][:, 1, :]
        b_snapshot = file["timeseries/b/$iteration"][:, 1, :]
        ν_snapshot = file["timeseries/νₑ/$iteration"][:, 1, :]

        eigenmode_plot = eigenplot(ω_snapshot, b_snapshot, nothing, t; ω_lim=1, b_lim=0.05)

        plot(eigenmode_plot, size=size)
    end

    videofile_name = replace(jld2writer.filepath, "jld2"=>"mp4")
    mp4(anim_total, videofile_name, fps = fps) # hide
end


closures = [AnisotropicMinimumDissipation(), SmagorinskyLilly()]
for closure in closures
    @info "\nStarting simulation with closure" closure
    global simulation = run_kelvin_helmholtz(closure, grid=grid, stop_time=150)
    plot_video(simulation, fps=14, size=(800, 300))
end
