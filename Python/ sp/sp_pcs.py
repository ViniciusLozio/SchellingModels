from functions_sp import coefs, total_ins, condsiniciais, execute_w, row, pc_finder
import csv
from os import mkdir, remove
import numpy as np

# This script runs the program for a chosen initial rho range
# and searches for critical points if phase transitions occur in the range.

# -------------- Declaration of Global Variables --------------- #

V = 4
nviz = V + 1  # Number of cells in the region
zero = 1E-15

# -------------------------------------------------------------- #
regras = ["01100"]

for rule in regras:

    # ------------------- Setting Rules ---------------------- #
    regra1 = rule
    regra0 = regra1[::-1]
    regras = [regra0, regra1]

    try:
        mkdir(f"./Output/pcs_sp/{regra1}")
    except FileExistsError:
        continue

    # ------------------------- Generating Coefficients ----------------------------#
    N = coefs(nviz, regras)
    u = N[0]
    w = N[1]
    x = N[2]
    y = N[3]
    z = N[4]

    # ------------------------- Setting Loop for Initial rhos ----------------------------#

    p_range = [0.1, 0.2, 0.3, 0.4, 0.5]
    ativa_lista = []

    for r in p_range:
        # Define r1 and r0 and the number of agents of each type:
        r1 = r
        r0 = 1 - r1

        # Create folders organized by r1:
        try:
            mkdir(f"./Output/pcs_sp/{regra1}/{r1:.12f}")
        except FileExistsError:
            continue

        # File detailing the stop condition for each rule and r1:
        with open(f"./Output/pcs_sp/{regra1}/{r1:.12f}/end.txt", "w") as end:

            # ------------------- Setting Execution Parameters ---------------------- #
            f = np.empty((2, nviz), dtype=float)  # Fraction of type-s cells with k living neighbors
            soma = np.empty(2, dtype=float)       # Total number of cells of each type

            t = 0.   # Time
            dt = np.sqrt(2)*1E-3  # Time step dt
            tmax = 20  # Maximum execution time

            print(f"\n regra: {regra1}", f"| r1: {r1:.12f}\n")

            condsiniciais(nviz, r1, f)
            dy = f.flatten().copy()

            # ----------------------- Configuring Writers (CSV) ------------------------ #
            total_ins_s = [total_ins(nviz, regras, dy)[0], total_ins(nviz, regras, dy)[1]]

            with open(f"./Output/pcs_sp/{regra1}/{r1:.12f}/s0.csv", 'w') as csvfile0,                  open(f"./Output/pcs_sp/{regra1}/{r1:.12f}/s1.csv", 'w') as csvfile1:
                writer0 = csv.writer(csvfile0, delimiter=',')
                writer0.writerow(row(nviz, total_ins_s, t, f, soma, 0, True))

                writer1 = csv.writer(csvfile1, delimiter=',')
                writer1.writerow(row(nviz, total_ins_s, t, f, soma, 1, True))

                writer0.writerow(row(nviz, total_ins_s, t, f, soma, 0))
                writer1.writerow(row(nviz, total_ins_s, t, f, soma, 1))

                # ----------------------------- Main Code  ---------------------------------- #
                T = execute_w(nviz, regras, N, dt, tmax, r1, f, [writer0, writer1], cp_find=False)

                dy = f.flatten().copy()

                if T < tmax:
                    ativa = 0
                    parada = "pu = 0"
                else:
                    ativa = 1
                    parada = f"T = tmax = {tmax}"

                total_ins_s = [total_ins(nviz, regras, dy)[0], total_ins(nviz, regras, dy)[1]]

                end.write(
                    f"Dinamica Interrompida: {parada} | t = {T} | ρu0 = {total_ins_s[0]:.14f} | "
                    f"ρu1 = {total_ins_s[1]:.14f} | ρu: {total_ins_s[0] * r0 + total_ins_s[1] * r1:.14f}\n")
                end.write(f"Tipo 0: {f[0]}\n")
                end.write(f"Tipo 1: {f[1]}\n")

        print("")
        print(f[0], soma[0])
        print(f[1], soma[1])
        print(f"\nT: {T:.10f}", f"| Parada: {parada}, "f"| regra: {regra1}",
              f"| r1: {r1:.10f}",
              f"| Ins 0: {total_ins_s[0]:.16f}",
              f"| Ins 1: {total_ins_s[1]:.16f} | Total Ins: "
              f"{(total_ins_s[0] * r0 + total_ins_s[1] * r1):.14f}")
        print(
            "\n "
            "-------------------------------------------------------------------------------------------------------------------------------------------------------------- ")

        # --------------------------------- Finding Critical Points  --------------------------------------- #
        if len(ativa_lista) > 0:
            if ativa != ativa_lista[-1]:
                try:
                    mkdir(f"./Output/pcs_sp/{regra1}/pcs")
                except FileExistsError:
                    pass

                if ativa == 1:
                    r_ativa = r
                    r_fr = r - 0.1
                else:
                    r_ativa = r - 0.1
                    r_fr = r

                print(f"\nStarting critical point search in ]{r - 0.1}, {r}[")
                pc = pc_finder(nviz, regras, N, dt, tmax, r_ativa, r_fr, [writer0, writer1])
                r_b, emax = pc

                with open(f"./Output/pcs_sp/{regra1}/pcs/pcs.txt", "a") as txt:
                    txt.write(f"{r_b:.12f} ± {emax:.12f}\n")

        ativa_lista.append(ativa)
