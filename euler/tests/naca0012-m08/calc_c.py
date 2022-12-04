import numpy as np
import matplotlib.pyplot as plt

def read_plot_file(filename):

    d = {}
    
    with open(filename, "rt") as f:
        lines = f.readlines()
        for i, l in enumerate(lines):
            if i == 0:
                if "\t" in l:
                    ls = l.split("\t")
                else:
                    ls = l.split(",")
                for key in ls:
                    key = key.replace("\n", "").replace(" ", "")
                    d[key] = []
            else:
                if "\t" in l:
                    ls = l.split("\t")
                else:
                    ls = l.split(",")
                for j, key in enumerate(d.keys()):
                    d[key].append(float(ls[j]))
    
    for key in d.keys():
        d[key] = np.array(d[key])
    return d

def naca0012(x, u):
    return u*0.594689181*(0.298222773*np.sqrt(x) - 0.127125232*x - 0.357907906*x**2 + 0.291984971*x**3 - 0.105174606*x**4)

def add_y(d):
    d['y'] = []
    
    for i in range(len(d['x'])):
        x = d['x'][i]
        if i < 5:
            y = naca0012(x, 1)
        elif d['x'][i-1] < x:
            y = naca0012(x, 1)
        else:
            y = naca0012(x, -1)
        d['y'].append(y)
    return d

def rescale_x(d):
    x_min = np.min(d['x'])
    x_max = np.max(d['x'])
    d['x'] = (d['x'] - x_min)/(x_max - x_min)
    return d

def integrate(x, y, v, s=1):
    V = np.array([0., 0.])
    M = 0.
    center = np.array([0.25, 0.])
    for i in range(len(x)-1):
        r = np.array([x[i] + x[i+1], y[i] + y[i+1]])*0.5
        dxi = x[i+1] - x[i]
        dyi = y[i+1] - y[i]
        vm = (v[i] + v[i+1])*0.5
        # Normal vector
        nv = np.array([-1*dyi, dxi])
        V += s*nv * vm
        M -= s*np.cross(r - center, nv*vm)
    return V[0], V[1], M


files = [
    "Data_1.csv",
    "Data_3.csv",
    "Data_4.csv",
    "Data_5.csv"
]

print("CD, CL, CM")
for f in files:
    d = read_plot_file("Data_1.csv")
    d = rescale_x(d)
    d = add_y(d)

    aoa = 1.25 * np.pi/180

    D, L, CM = integrate(d['x'], d['y'], d['Cp'])
    CD = D*np.cos(aoa) + L*np.sin(aoa)
    CL = -D*np.sin(aoa) + L*np.cos(aoa)

    print(CD, CL, CM, sep=", ")

