import numpy as np
import euler


# Load mesh from gmsh
mesh = euler.mesh("naca0012.su2")

# Define problem constants
constants = euler.constants(
    gamma=1.4
)

# Generate initial conditions
mach = 0.8
aoa = 1.25 * np.pi/180.
q_inf = euler.var(
    rho=constants.gamma, 
    u=mach*np.cos(aoa), 
    v=mach*np.sin(aoa), 
    p=1.
)

q0 = euler.vectorVar()
for i in range(mesh.size()):
    q0.append(q_inf)

# Define solver with initial condition
solver = euler.steadyRk5Solver(
    mesh=mesh,
    init=q0, 
    constants=constants
)

# Set solver properties
solver.set_cfl(2.5)
solver.set_print_interval(500)
solver.set_tolerance(1e-16)
solver.set_max_steps(16000)
solver.set_order(2)
solver.set_smoother("implicit", 0.6, 2)
solver.set_limiter("venkatakrishnan")
solver.set_coefficient_save_field("airfoil", q_inf)

# Define boundary conditions
solver.set_bc("airfoil", "wall")
solver.set_bc("farfield", "farfield", q_inf)

# Run simulation
solver.simulate()

log = solver.log()


# Write all results to file
solver.writeVtk("results.vtk")

# Write residuals to file
solver.writeResiduals("residuals.data")

# Write variables on airfoil to file
solver.writeField("airfoil", "airfoil.data", q_inf)

# Write coefficients
solver.writeSavedCoefficients("coeffs.data")


# Compute and print coefficients
coeffs = solver.computeCoefficients("airfoil", q_inf)

print("(CD, CL, CM) =", coeffs)

# Write solver log to file
log += "\n"
log += "Aerodynamic coefficients\n"
log += " - CD = " + str(coeffs[0]) + "\n"
log += " - CL = " + str(coeffs[1]) + "\n"
log += " - CM = " + str(coeffs[2]) + "\n"

with open("solver.log", "wt") as f:
    f.write(log)
