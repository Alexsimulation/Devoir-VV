import numpy as np
import matplotlib.pyplot as plt
import euler

# Load mesh from gmsh
mesh = euler.mesh("line.su2", True)

# Define problem constants
constants = euler.constants(
    gamma=1.4
)

# Generate initial conditions
q_left = euler.var(
    rho=1., 
    u=0., 
    v=0., 
    p=1.
)
q_right = euler.var(
    rho=0.125, 
    u=0., 
    v=0., 
    p=0.1
)

q0 = euler.vectorVar()
for i in range(len(mesh)):
    x, y = mesh.get_xy(i)
    if (x < 0.5):
        q0.append(q_left)
    else:
        q0.append(q_right)

# Define solver with initial condition
solver_o1 = euler.transientEulerSolver(
    mesh=mesh,
    init=q0, 
    constants=constants
)

# Set solver properties
solver_o1.set_cfl(0.1)
solver_o1.set_print_interval(1000)
solver_o1.set_max_time(0.2)
solver_o1.set_max_steps(100000)
solver_o1.set_order(1)


# Define boundary conditions
solver_o1.set_bc("empty", "wall")
solver_o1.set_bc("left", "farfield", q_left)
solver_o1.set_bc("right", "farfield", q_right)


# Run simulation
solver_o1.simulate(True)

with open("solver.log", "wt") as f:
    f.write(solver_o1.log())


solver_o1.writeVtk("test.vtk")

# Get results into vector
sol_o1 = solver_o1.get_solution()

rho1 = np.zeros((len(sol_o1), 1))
x = np.zeros((len(sol_o1), 1))
for i in range(len(rho1)):
    rho1[i] = sol_o1[i][0]
    x[i], _ = mesh.get_xy(i)



fig, ax = plt.subplots(nrows=1, ncols=1)


ax.plot(x, rho1, label="Order 1")
ax.legend()


ax.set_xlabel("X position")
ax.set_ylabel("Density")



plt.show()
