import csv
from functions import coefs, total_ins, row, derivs, condsiniciais
from os import mkdir, remove
from math import comb, floor
from random import uniform
from tqdm import tqdm
import numpy as np
import scipy as sp
from numba import jit
import numba_scipy

# -------------- Declaration of Global Variables --------------- #

V = 4
nviz = V + 1  # Number of cells in the region
zero = 0.
dt = np.sqrt(2)*1E-3   # Time step interval
tmax = 15  # Maximum execution time of the program

# ------------------- Setting Rules and Integration Parameters ---------------------- #

f = np.empty((2, nviz), dtype=float)  # Fraction of type-s cells with k living neighbors
soma = np.empty(2, dtype=float)  # Total number of cells of each type
t0 = 0  # Time
tf = tmax
step = dt
max_num = tf / step

rules = ["01100"]
pcs = ["0.363606306911"]

for i in range(0, len(rules)):

    p_range_menos = []
    p_range_mais = []

    regra1 = rules[i]
    regra0 = regra1[::-1]
    regras = [regra0, regra1]

    pc = float(pcs[i])

    p_vs_pc = [10**-6, 10**-5, 5*10**-5, 10**-4, 5*10**-4, 10**-3, 5*10**-3, 10**-2, 5*10**-2]

    for dp in p_vs_pc:
        p_range_mais.append(pc + dp)
        p_range_menos.append(pc - dp)

    lst = [p_range_menos, p_range_mais]

    try:
        mkdir(f"./Output/desvios_sp/{regra1}")
    except FileExistsError:
        pass

    # ------------------------- Generating Coefficients ----------------------------#
    N = coefs(nviz, regras)
    u = N[0]
    w = N[1]
    x = N[2]
    y = N[3]
    z = N[4]

    # ------------------- Setting Execution Parameters ---------------------- #
    for p_range in lst:

        if lst.index(p_range) == 0:
            str_desvio = "menos"
        else:
            str_desvio = "mais"

        try:
            mkdir(f"./Output/desvios_sp/{regra1}/{str_desvio}")
        except FileExistsError:
            pass

        with open(
            f"./Output/desvios_sp/{regra1}/{str_desvio}/ppcx.csv",
            'a') as ppcx:
            ppcx_w = csv.writer(ppcx, delimiter=',')
            ppcx_w.writerow(["|p - pc|", "T", "pu_final_1", "pu_final_0"])

            for r1 in p_range:

                try:
                    mkdir(f"./Output/desvios_sp/{regra1}/{str_desvio}/{p_vs_pc[p_range.index(r1)]}")
                except FileExistsError:
                    pass

                print(f"\n regra: {regra1}", f"| r1: {r1:.12f}\n")
                condsiniciais(nviz, r1, f)
                dy = f.flatten().copy()
                soma = f.sum(axis=1)

                # ----------------------- Configuring Writers (CSV) ------------------------ #
                total_ins_s = [total_ins(nviz, regras, dy)[0], total_ins(nviz, regras, dy)[1]]

                with open(f"./Output/desvios_sp/{regra1}/{str_desvio}/{p_vs_pc[p_range.index(r1)]}/s0.csv", 'w') as csvfile0,                      open(f"./Output/desvios_sp/{regra1}/{str_desvio}/{p_vs_pc[p_range.index(r1)]}/s1.csv", 'w') as csvfile1:

                    writer0 = csv.writer(csvfile0, delimiter=',')
                    writer0.writerow(row(nviz, total_ins_s, t0, f, soma, 0, True))

                    writer1 = csv.writer(csvfile1, delimiter=',')
                    writer1.writerow(row(nviz, total_ins_s, t0, f, soma, 1, True))

                    writer0.writerow(row(nviz, total_ins_s, t0, f, soma, 0))
                    writer1.writerow(row(nviz, total_ins_s, t0, f, soma, 1))

                    # ----------------------- Main Code  ------------------------ #
                    def deriva(t, y, r1):
                        return derivs(nviz, regras, y, N, 1. - r1)

                    def stop(t, y, r1):
                        total_ins_s = [total_ins(nviz, regras, y)[0], total_ins(nviz, regras, y)[1]]
                        if (total_ins_s[0] <= 0. or total_ins_s[1] <= 0.) or np.any(y < 0):
                            return 0.
                        else:
                            return 1

                    stop.terminal = True
                    yout = sp.integrate.solve_ivp(deriva, [t0, tf], dy, args=(r1,), max_step=dt,
                                                  method="DOP853",
                                                  events=stop)

                    for i in range(0, len(yout.t)-1):
                        f = np.array([yout.y.T[i][0:5], yout.y.T[i][5:10]])
                        soma = f.sum(axis=1)
                        total_ins_s = [total_ins(nviz, regras, yout.y.T[i])[0], total_ins(nviz, regras, yout.y.T[i])[1]]
                        writer0.writerow(row(nviz, total_ins_s, yout.t[i], f, soma, 0))
                        writer1.writerow(row(nviz, total_ins_s, yout.t[i], f, soma, 1))

                    print(yout.t[-1])
                    ppcx_w.writerow([f"{p_vs_pc[p_range.index(r1)]:.10f}", yout.t[-1], total_ins_s[1], total_ins_s[0]])
