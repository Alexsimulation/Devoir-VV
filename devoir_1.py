import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Fonction qui résous le problème pour un nombre d'élément N et un ordre donné
def solve_N(N, plot=False, ordre=1, constant_source=True, extra_source=lambda r, t : 0, tf=1e200, init=["none"]):
    # Définition de l'équation différentielle
    r, t, c, cr, crr, ct = sp.symbols("r, t, c, cr, crr, ct")
    k, S, D = sp.symbols("k, S, D")

    L = -ct + D*(1/r*cr + crr) - k*c - S

    # Définition du schéma de différentiation
    c_tp1_rp1, c_tp1, c_tp1_rm1, c_t = sp.symbols("c_tp1_rp1, c_tp1, c_tp1_rm1, c_t")
    dr, dt = sp.symbols("dr, dt")

    if ordre == 1:
        L = L.subs(cr, (c_tp1_rp1 - c_tp1)/dr)
    else:
        L = L.subs(cr, (c_tp1_rp1 - c_tp1_rm1)/(2*dr))
    L = L.subs(crr, (c_tp1_rp1 - 2*c_tp1 + c_tp1_rm1)/dr**2)
    L = L.subs(ct, (c_tp1 - c_t)/dt)
    L = L.subs(c, c_tp1)
    L = L.expand()


    # Génération des fonctions pour chaque coefficients
    a_ii   = sp.lambdify([r, dr, dt, k, S, D], L.coeff(c_tp1))
    a_iim1 = sp.lambdify([r, dr, dt, k, S, D], L.coeff(c_tp1_rm1))
    a_iip1 = sp.lambdify([r, dr, dt, k, S, D], L.coeff(c_tp1_rp1))

    b_i = sp.lambdify([r, dr, dt, k, S, D, c_t], -sp.simplify(L-L.coeff(c_tp1)*c_tp1-L.coeff(c_tp1_rm1)*c_tp1_rm1-L.coeff(c_tp1_rp1)*c_tp1_rp1))

    # Discrétisation du domaine et définition des coefficients
    R = 0.5
    if constant_source:
        k = 0 # 4e-9
        S = 1e-8 # 1e-8    
    else:
        k = 4e-9 # 4e-9
        S = 0 # 1e-8
    D = 1e-10
    C_e = 10

    r = np.linspace(0, R, N)
    dt = 1e4
    dr = r[1] - r[0]

    c = np.zeros(N)
    if init[0] != "none":
        c = init[1](r)


    # Initialisation de la figure
    if plot:
        fig, ax = plt.subplots()
        ax.set_title("Solution")

    # Boucle temporelle
    residual = 1.
    tolerance = 1e-14
    time = 0
    while (residual > tolerance)&(time < tf):
        # Remplissage de la matrice
        A = np.zeros((N, N))
        b = np.zeros((N, 1))
        for i in range(1, N-1):
            A[i, i] = a_ii(r[i], dr, dt, k, S, D)
            A[i, i-1] = a_iim1(r[i], dr, dt, k, S, D)
            A[i, i+1] = a_iip1(r[i], dr, dt, k, S, D)

            b[i] = b_i(r[i], dr, dt, k, S, D, c[i]) + extra_source(r[i], time)

        # Conditions limites
        A[N-1, N-1] = 1
        b[N-1] = C_e

        if ordre == 1:
            A[0, 0] = -1
            A[0, 1] = 1
            b[0] = 0
        else:
            A[0, 0] = -3
            A[0, 1] = 4
            A[0, 2] = -1
            b[0] = 0

        # Résoudre
        dc = np.linalg.solve(A, b).flatten() - c
        c = c + dc

        # Calcul du résidu
        residual = np.linalg.norm(dc)

        time += dt

        if plot:
            ax.cla()
            ax.set_title("Solution")
            ax.plot(r, c)
            plt.pause(1e-4)



    # Calcul de l'erreur selon
    c_theo = 0.25*S/D*R**2*(r**2/R**2-1) + C_e
    erreur = c - c_theo

    if plot:
        ax.plot(r, c_theo)

        fig, ax = plt.subplots()
        ax.plot(r, erreur)
        ax.set_title("Erreur absolue")

    L1 = np.trapz(np.abs(c - c_theo), r)/R
    L2 = np.sqrt(np.trapz((c - c_theo)**2, r)/R)
    Linf = np.max(np.abs(c - c_theo))
    return {
        "L1":L1,
        "L2":L2,
        "Linf":Linf,
        "r":r,
        "c":c,
        "t":time
    }


# Fonction qui résous pour N=5 et ordre=1
def solve_and_show(n=5, ordre=1):
    solve_N(n, plot=True, ordre=ordre)
    plt.show()


# Fonction qui génère un graphique de l'erreur selon le nombre d'éléments
def plot_error(ordre=1):
    N = [round(5*10**(t/4)) for t in range(1, 8)]
    print("Nombre d'éléments:",N)

    es = {"L1":[], "L2":[], "Linf":[]}
    for n in N:
        e = solve_N(n, plot=False, ordre=ordre)
        print(n, ":", e)
        es["L1"].append(e["L1"])
        es["L2"].append(e["L2"])
        es["Linf"].append(e["Linf"])

    ordre = (np.log(es["L2"][0]) - np.log(es["L2"][len(N)-1]))/(np.log(N[len(N)-1]) - np.log(N[0]))
    print("Ordre = ",ordre)
    fig, ax = plt.subplots()
    ls = []
    for key, value in es.items():
        ax.plot(N, value)
        ls.append(key)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(ls)
    ax.set_xlabel("Nombre d'éléments")
    ax.set_ylabel("Norme de l'erreur")
    plt.show()

if __name__=="__main__":
    solve_and_show(n=50, ordre=2)
