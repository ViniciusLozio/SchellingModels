import csv
from functions import coefs, total_ins, row, derivs, condsiniciais
import numpy as np
import scipy as sp

# -------------- Declaration of Global Variables --------------- #

V = 4
nviz = V + 1  # Number of cells in the region
zero = 0.
dt = 1E-2  # Integration time step dt
tmax = 10  # Maximum execution time of the program

# ------------------- Setting Rules ---------------------- #

regra1 = "10000"
regra0 = regra1[::-1]
regras = [regra0, regra1]

r1 = 0.5
r0 = 1 - r1

print(f"\n regra: {regra1}", f"| r1: {r1:.12f}\n")

# ------------------------- Generating Coefficients ----------------------------#
N = coefs(nviz, regras)
u = N[0]
w = N[1]
x = N[2]
y = N[3]
z = N[4]

# ------------------- Setting Execution Parameters ---------------------- #

f = np.empty((2, nviz), dtype=float)  # Fraction of type-s cells with k living neighbors
soma = np.empty(2, dtype=float)  # Total number of cells of each type
t0 = 0  # Time
condsiniciais(nviz, r1, f)

dy = f.flatten().copy()
soma = f.sum(axis=1)

# ----------------------- Configuring Writers (CSV) ------------------------ #

total_ins_s = [total_ins(nviz, regras, dy)[0], total_ins(nviz, regras, dy)[1]]

with open("./Output/Integration/s0.csv", 'w') as csvfile0, open("./Output/Integration/s1.csv", 'w') as csvfile1, open(
        "./Output/Derivation/s0_deriv.csv", 'w') as csvfile2, open("./Output/Derivation/s1_deriv.csv", 'w') as csvfile3:
    writer0 = csv.writer(csvfile0, delimiter=',')
    writer0.writerow(row(nviz, total_ins_s, t0, f, soma, 0, True))

    writer1 = csv.writer(csvfile1, delimiter=',')
    writer1.writerow(row(nviz, total_ins_s, t0, f, soma, 1, True))

    writer2 = csv.writer(csvfile2, delimiter=',')
    writer3 = csv.writer(csvfile3, delimiter=',')

    writer0.writerow(row(nviz, total_ins_s, t0, f, soma, 0))
    writer1.writerow(row(nviz, total_ins_s, t0, f, soma, 1))


    # ----------------------- Main Code  ------------------------ #

    def deriva(t, y):
        return derivs(nviz, regras, y, N, 1. - r1)


    def stop(t, y):
        total_ins_s = [total_ins(nviz, regras, y)[0], total_ins(nviz, regras, y)[1]]
        if (total_ins_s[0] <= 0. or total_ins_s[1] <= 0.) or np.any(y < 0):
            return 0.
        else:
            return 1.


    stop.terminal = True

    yout = sp.integrate.solve_ivp(deriva, [t0, tmax], dy, max_step=dt, method="DOP853", events=stop, dense_output=True)

    for i in range(0, len(yout.t)):
        f = np.array([yout.y.T[i][0:5], yout.y.T[i][5:10]])
        soma = f.sum(axis=1)
        total_ins_s = [total_ins(nviz, regras, yout.y.T[i])[0], total_ins(nviz, regras, yout.y.T[i])[1]]
        writer0.writerow(row(nviz, total_ins_s, yout.t[i], f, soma, 0))
        writer1.writerow(row(nviz, total_ins_s, yout.t[i], f, soma, 1))

    print("")
    print(f"\nT: {yout.t[-1]:.10f}", f"| regra: {regra1}",
          f"| r1: {r1:.10f}",
          f"| Ins 0: {total_ins_s[0]:.16f}",
          f"| Ins 1: {total_ins_s[1]:.16f}")
    print(
        "\n "
        "-------------------------------------------------------------------------------------------------------------------------------------------------------------- ")
