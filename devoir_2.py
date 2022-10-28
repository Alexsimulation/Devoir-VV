from tkinter import E
from devoir_1 import solve_N
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def comsol(N, ordre=1, plot=False):
    def read_comsol(filename):
        d = {"r":[], "c":[]}
        with open(filename, "rt") as f:
            lines = f.readlines()
            for l in lines:
                if l[0] != "%":
                    ls = " ".join(l.split())
                    ls = ls.split(" ")
                    d["r"].append(float(ls[0].replace("\n", "")))
                    d["c"].append(float(ls[1].replace("\n", "")))
        
        d["r"] = np.array(d["r"])
        d["c"] = np.array(d["c"])
        return d


    d = solve_N(
        N=N, 
        plot=False, 
        ordre=ordre, 
        constant_source=False
    )


    comsol = read_comsol("data.txt")

    if plot:
        fig, ax = plt.subplots()

        ax.plot(comsol["r"], comsol["c"], label="Comsol")
        ax.plot(d["r"], d["c"], label="Notre code")
        ax.legend()

        plt.show()

    c_theo = []
    for ri in d["r"]:
        c_theo.append(np.interp(ri, comsol["r"], comsol["c"]))
    c_theo = np.array(c_theo)

    pi = sp.pi
    R = 0.5
    D = 1e-10
    K = 4e-9
    CE = 10
    r = d["r"]
    c = d["c"]
    
    L2 = np.sqrt(np.trapz((c - c_theo)**2, r)/R)
    return L2


def mms(N, ordre=1, plot=False, dt=1e6):
    pi = sp.pi
    R = 0.5
    D = 1e-10
    K = 4e-9
    CE = 10


    r, t = sp.symbols("r, t")

    c = CE - sp.cos(pi/2*r/R) * sp.exp(t/1e10)
    cr = sp.diff(c, r)
    crr = sp.diff(cr, r)
    ct = sp.diff(c, t)

    L = -ct + D*(1/r*cr + crr) - K*c
    print(L)
    L = sp.lambdify([r, t], L)
    MS = sp.lambdify([r, t], c)

    d = solve_N(
        N=N, 
        plot=False, 
        ordre=ordre,
        constant_source=False,
        extra_source=L,
        tf=1e8,
        init=["func", lambda r : MS(r, 0)],
        c_theo=["func",MS],
        dt=dt
    )

    d["sc"] = []
    for ri in d["r"]:
        d["sc"].append(MS(ri, d["t"]))
    d["sc"] = np.array(d["sc"])

    if plot:
        fig, ax = plt.subplots()

        ax.plot(d["r"], d["c"], label="Notre code")
        ax.plot(d["r"], d["sc"], label="MS")
        
        ax.legend()

        plt.show()
    
    return d["L2"]


if __name__=="__main__":

    # Comsol

    N = np.array([5, 10, 15, 20])
    L2_o1s = []
    L2_o2s = []
    for n in N:
        print(n)
        L2_o1s.append( comsol(n, 1) )
        L2_o2s.append( comsol(n, 2) )
    
    L = len(N) - 1
    ordre_o1 = (np.log(L2_o1s[0]) - np.log(L2_o1s[L]))/(np.log(N[L]) - np.log(N[0]))
    ordre_o2 = (np.log(L2_o2s[0]) - np.log(L2_o2s[L]))/(np.log(N[L]) - np.log(N[0]))
    print(ordre_o1)
    print(ordre_o2)

    fig, ax = plt.subplots()
    ax.plot(1/N, L2_o1s, label="Order 1")
    ax.plot(1/N, L2_o2s, label="Order 2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("Delta x")
    ax.set_ylabel("Norme de l'erreur")
    ax.set_title("Ordre de la solution vs Comsol avec variation de dx")

    # MMS dx variation

    N = np.array([5, 10, 15, 20])
    L2_o1s = []
    L2_o2s = []
    for n in N:
        L2_o1s.append( mms(n, 1, dt=1e5) )
        L2_o2s.append( mms(n, 2, dt=1e5) )
    
    L = len(N) - 1
    ordre_o1 = (np.log(L2_o1s[0]) - np.log(L2_o1s[L]))/(np.log(N[L]) - np.log(N[0]))
    ordre_o2 = (np.log(L2_o2s[0]) - np.log(L2_o2s[L]))/(np.log(N[L]) - np.log(N[0]))
    print(ordre_o1)
    print(ordre_o2)

    fig, ax = plt.subplots()
    ax.plot(1/N, L2_o1s, label="Order 1")
    ax.plot(1/N, L2_o2s, label="Order 2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("Delta x")
    ax.set_ylabel("Norme de l'erreur")
    ax.set_title("Ordre de la solution MMS avec variation de dx")

    # MMS dt variation
    N = np.array([1e5, 2e5, 5e5, 1e6])
    L2_o1s = []
    L2_o2s = []
    for n in N:
        L2_o1s.append( mms(100, 1, dt=n) )
        L2_o2s.append( mms(100, 2, dt=n) )
    
    L = len(N) - 1
    ordre_o1 = (np.log(L2_o1s[0]) - np.log(L2_o1s[L]))/(np.log(N[0]) - np.log(N[L]))
    ordre_o2 = (np.log(L2_o2s[0]) - np.log(L2_o2s[L]))/(np.log(N[0]) - np.log(N[L]))
    print(ordre_o1)
    print(ordre_o2)

    fig, ax = plt.subplots()
    ax.plot(N, L2_o1s, label="Order 1")
    ax.plot(N, L2_o2s, label="Order 2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("Pas de temps")
    ax.set_ylabel("Norme de l'erreur")
    ax.set_title("Ordre de la solution MMS avec variation de dt")
    plt.show()