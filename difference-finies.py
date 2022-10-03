import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# Définition de l'équation différentielle
r, t, c, cr, crr, ct = sp.symbols("r, t, c, cr, crr, ct")
k, D = sp.symbols("k, D")

L = -ct + D*(1/r*cr + crr) - k*c

# Définition du schéma de différentiation
c_tp1_rp1, c_tp1, c_tp1_rm1, c_t = sp.symbols("c_tp1_rp1, c_tp1, c_tp1_rm1, c_t")
dr, dt = sp.symbols("dr, dt")

L = L.subs(cr, (c_tp1_rp1 - c_tp1)/dr)
L = L.subs(crr, (c_tp1_rp1 - 2*c_tp1 + c_tp1_rm1)/dr**2)
L = L.subs(ct, (c_tp1 - c_t)/dt)
L = L.subs(c, c_tp1)
L = L.expand()


# Génération des fonctions pour chaque coefficients
a_ii   = sp.lambdify([r, dr, dt, k, D], L.coeff(c_tp1))
a_iim1 = sp.lambdify([r, dr, dt, k, D], L.coeff(c_tp1_rm1))
a_iip1 = sp.lambdify([r, dr, dt, k, D], L.coeff(c_tp1_rp1))

b_i = sp.lambdify([r, dr, dt, k, D, c_t], -sp.simplify(L-L.coeff(c_tp1)*c_tp1-L.coeff(c_tp1_rm1)*c_tp1_rm1-L.coeff(c_tp1_rp1)*c_tp1_rp1))

# Discrétisation du domaine et définition des coefficients
N = 200
R = 0.5
k = 4e-9
D = 1e-10
C_e = 10

r = np.linspace(0, R, N)
dt = 1e8
dr = r[1] - r[0]

c = np.zeros(N)

# Initialisation de la figure
fig, ax = plt.subplots()

# Boucle temporelle
residual = 1.
tolerance = 1e-3
while residual > tolerance:
    # Remplissage de la matrice
    A = np.zeros((N, N))
    b = np.zeros((N, 1))
    for i in range(1, N-1):
        A[i, i] = a_ii(r[i], dr, dt, k, D)
        A[i, i-1] = a_iim1(r[i], dr, dt, k, D)
        A[i, i+1] = a_iip1(r[i], dr, dt, k, D)

        b[i] = b_i(r, dr, dt, k, D, c[i])

    # Conditions limites
    A[N-1, N-1] = 1
    b[N-1] = C_e

    A[0, 0] = -3
    A[0, 1] = 4
    A[0, 2] = -1
    b[0] = 0

    # Résoudre
    dc = np.linalg.solve(A, b).flatten() - c
    c = c + dc

    # Calcul du résidu
    residual = np.linalg.norm(dc)
    print(residual)

    ax.cla()
    ax.plot(r, c)
    plt.pause(1e-4)

plt.show()
