import numpy as np
import euler

# Load mesh from gmsh
mesh = euler.mesh("circle.su2")

# Define problem constants
constants = euler.constants(
    gamma=1.4
)

# Generate initial conditions
mach = 2.
aoa = 0.0 * np.pi/180.
q_inf = euler.var(
    rho=constants.gamma, 
    u=mach*np.cos(aoa), 
    v=mach*np.sin(aoa), 
    p=1.
)
q_init = euler.var(
    rho=constants.gamma, 
    u=0.5*np.cos(aoa), 
    v=0.5*np.sin(aoa), 
    p=1.
)

q0 = euler.vectorVar()
for i in range(mesh.size()):
    q0.append(q_init)

# Define solver with initial condition
solver = euler.steadyRk5Solver(
    mesh=mesh,
    init=q0, 
    constants=constants
)

# Set solver properties
solver.set_cfl(5.)
solver.set_print_interval(100)
solver.set_tolerance(1e-5)
solver.set_max_steps(10000)
solver.set_order(1)
solver.set_smoother("implicit", 0.6, 2)

# Define boundary conditions
solver.set_bc("circle", "wall")
solver.set_bc("farfield", "farfield", q_inf)

# Run simulation
solver.simulate()

# Write all results to file
solver.writeVtk("circle.vtk")

