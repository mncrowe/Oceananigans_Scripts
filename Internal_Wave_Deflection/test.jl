using Oceananigans, NCDatasets

Nx, Ny, Nz = 256, 256, 256
Δt = 2π / Nx

grid = RectilinearGrid(GPU(), size=(Nx, Ny, Nz), x=(-π, π), y=(-π, π), z=(-π, π),
                            topology=(Periodic, Periodic, Periodic))

coriolis = FPlane(f=0.5)

N = 1
B_func(x, y, z, t, N) = N^2 * z
B = BackgroundField(B_func, parameters=N)

model = NonhydrostaticModel(grid; coriolis,
                            advection = Centered(order=4),
                            closure = ScalarDiffusivity(ν=1e-4, κ=1e-4),
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            background_fields = (; b=B))

m = 10
f = coriolis.f
ω = f

u₀(x, y, z) = (1 - 2*y^2) * exp(-(x^2 + y^2))
v₀(x, y, z) = - 2 * x * y * exp(-(x^2 + y^2))
w₀(x, y, z) = 0.1 * ω   * cos(m * z)
b₀(x, y, z) = 0.1 * N^2 * sin(m * z)

set!(model, u=u₀, v=v₀, w=w₀, b=b₀)

simulation = Simulation(model, Δt = Δt, stop_iteration = 100)

filename = "internal_wave.nc"
simulation.output_writers[:velocities] = NetCDFWriter(model, model.velocities; filename,
                                                    schedule = IterationInterval(10),
                                                    overwrite_existing = true,
                                                    indices=(:, :, :))    # indices=(:, :, grid.Nz)

run!(simulation)
