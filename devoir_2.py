from tkinter import E
from devoir_1 import solve_N
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def comsol():
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
        N=100, 
        plot=False, 
        ordre=2, 
        constant_source=False
    )


    comsol = read_comsol("data.txt")

    fig, ax = plt.subplots()

    ax.plot(comsol["r"], comsol["c"], label="Comsol")
    ax.plot(d["r"], d["c"], label="Notre code")
    ax.legend()

    plt.show()


def mms():
    pi = sp.pi
    R = 0.5
    D = 1e-10
    K = 4e-9
    CE = 10


    r, t = sp.symbols("r, t")

    c = CE - sp.cos(pi/2*r/R) * sp.exp(t/1e9)
    cr = sp.diff(c, r)
    crr = sp.diff(cr, r)
    ct = sp.diff(c, t)

    L = -ct + D*(1/r*cr + crr) - K*c
    print(L)
    L = sp.lambdify([r, t], L)
    MS = sp.lambdify([r, t], c)

    d = solve_N(
        N=10, 
        plot=False, 
        ordre=2,
        constant_source=False,
        extra_source=L,
        tf=1e7,
        init=["func", lambda r : MS(r, 0)]
    )

    d["sc"] = []
    for ri in d["r"]:
        d["sc"].append(MS(ri, d["t"]))
    d["sc"] = np.array(d["sc"])

    fig, ax = plt.subplots()

    ax.plot(d["r"], d["c"], label="Notre code")
    ax.plot(d["r"], d["sc"], label="MS")
    
    ax.legend()

    plt.show()


if __name__=="__main__":
    mms()