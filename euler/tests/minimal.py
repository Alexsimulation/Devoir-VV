# A minimal example of the euler module
import euler

# Read mesh
mesh = euler.gmshToFvm("naca0012-m08/naca0012-coarse.msh")

# Define problem constants
constants = euler.constants(gamma=1.4)

# Generate initial conditions
q_inf = euler.var(rho=1.4, u=0.8, v=0., p=1.)

q0 = euler.vectorVar()
for i in range(mesh.size()):
    q0.append(q_inf)

# Define solver with initial condition
solver = euler.steadyRk5Solver(mesh=mesh, init=q0, constants=constants)

# Set solver properties
solver.set_cfl(1.)
solver.set_print_interval(1000)
solver.set_tolerance(1e-16)
solver.set_max_steps(30000)
solver.set_order(2)
solver.set_smoother("implicit", 0.6, 2)
solver.set_limiter("venkatakrishnan")

# Define boundary conditions
solver.set_bc("airfoil", "wall")
solver.set_bc("farfield", "farfield", q_inf)

# Run simulation
solver.simulate()

# Write all results to file
solver.writeVtk("results-o2-coarse.vtk")

# Write residuals to file
solver.writeResiduals("residuals-o2-coarse.plot")

# Write variables on airfoil to file
solver.writeField("airfoil", "airfoil-o2-coarse.plot", q_inf)
