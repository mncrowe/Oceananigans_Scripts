using Oceananigans, NCDatasets, CUDA

# Set parameters:

ν = 0.01
κS, κT = 0.04, 0.01
ΔT, ΔS = 2, 1
N = 1e5			# number of iterations
Ns, Np = 200, 100	# number of saves, number of progress updates
filename = "big_S_diff"

Nx, Nz = 512, 101
Lx, Lz = 100, 1
Δt = 0.5 * (Lz / Nz)^2 / max(ν, κS, κT)    # diffusive CFL

Lf = Lx / 20

grid = RectilinearGrid(GPU(),
                       size = (Nx, Nz),
                       x = (-Lx/2, Lx/2),
                       z = (-Lz/2, Lz/2),
                       topology = (Periodic, Flat, Bounded))

# NonhydrostaticModel

parameters = (ΔT = ΔT, ΔS = ΔS, Lx = Lx)
Tb(x, z, t, p) = 2 * p.ΔT * x / p.Lx
Sb(x, z, t, p) = 2 * p.ΔS * x / p.Lx

T_background = BackgroundField(Tb, parameters = parameters)
S_background = BackgroundField(Sb, parameters = parameters)

model = NonhydrostaticModel(grid;
                            advection = Centered(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S),
                            buoyancy = SeawaterBuoyancy(gravitational_acceleration = 1,
                                                        equation_of_state = LinearEquationOfState(thermal_expansion = 1, 
                                                                                                  haline_contraction = 1)),
                            background_fields = (T = T_background, S = S_background),
                            closure = VerticalScalarDiffusivity(; ν = ν, κ = (S = κS, T = κT)))

#T₀(x, z) = ΔT * (1 + tanh((x - Lx/4) / Lf) - tanh((x + Lx/4) / Lf))
#S₀(x, z) = ΔS * (1 + tanh((x - Lx/4) / Lf) - tanh((x + Lx/4) / Lf))

T₀(x, z) = ΔT * tanh(x / Lf) - 2 * ΔT * x / Lx
S₀(x, z) = ΔS * tanh(x / Lf) - 2 * ΔS * x / Lx

set!(model, T = T₀, S = S₀)

simulation = Simulation(model, Δt = Δt, stop_iteration = N)

print_progress(sim) = @info "Iter $(iteration(sim)): $(prettytime(sim))"
add_callback!(simulation, print_progress, IterationInterval(Int(N / Np)), name=:progress)

simulation.output_writers[:velocities] = NetCDFWriter(model, model.velocities; filename = filename * "_velocity.nc",
                                                    schedule = IterationInterval(Int(N / Ns)),
                                                    overwrite_existing = true)

simulation.output_writers[:tracers] = NetCDFWriter(model, model.tracers; filename = filename * "_tracers.nc",
                                                    schedule = IterationInterval(Int(N / Ns)),
                                                    overwrite_existing = true)

run!(simulation)
