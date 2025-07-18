from math import comb, floor
from random import uniform
import scipy as sp
import numpy as np
from tqdm import tqdm
from os import mkdir

"""---------------------------------- Functions for coefficient calculation ---------------------------------"""


def estados(nviz: int, n: int):
    """
    Returns a list representing neighborhood state n.
    Ex: [0,1,1,1,0] (The first value indicates the state of the central cell)

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param n: int (State index)

    :return e0: list
    """

    e0 = [0]
    aux = n
    for i in range(1, nviz):
        e0.append(floor((aux % 2)))
        aux = aux / 2

    return e0


def vizinhos(nviz: int):
    """
    Returns a list of how many neighbors in common each neighboring cell 
    shares with the central cell, and a list of their positions (indices).

    :param nviz: int (Number of cells in the region – Neighborhood + Central)

    :return L: list
    :return v: list
    """
    V = nviz - 1
    L = []
    v = []

    # Bethe lattice:
    # L = [1, 1, 1, 1, ...]
    # v = [[0], [0], [0], [0], ...]

    for i in range(0, V):
        L.append(1)
        v.append([0])

    return L, v


def zera(nviz: int):
    """
    Creates a three-dimensional matrix (s, k, l) initialized with zeros.

    :param nviz: int (Number of cells in the region – Neighborhood + Central)

    :return aux: list
    """

    aux = []

    for s in range(0, 2):
        aux.append([])
        for k in range(0, nviz):
            aux[s].append([])
            for l in range(0, nviz):
                aux[s][k].append([])
                for h in range(0, nviz):
                    aux[s][k][l].append(0)

    return aux


def sat(rule: str, nx: int):
    """
    Returns True or False depending on whether the specified cell
    is satisfied or not.

    :param nx: int (Cell position)
    :param rule: str (Satisfaction rule)

    :return: bool
    """
    return rule[nx] == '1'


def minmaxvizvivos(nviz: int, e: list, L: list, v: list):
    """
    Returns two lists containing the minimum and maximum number of same‑type neighbors
    that each neighbor of the central cell can have.

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param e: list (State)
    :param L: list (Number of common neighbors with the central cell)
    :param v: list (Positions (indices) of the cells)

    :return mink: list
    :return maxk: list
    """

    V = nviz - 1
    mink = []
    maxk = []

    for i in range(0, V):
        n0 = 0
        n1 = 0
        for j in range(0, L[i]):
            if e[v[i][j]] == 0:
                n0 += 1
            elif e[v[i][j]] == 1:
                n1 += 1
        mink.append(n1)
        maxk.append(V - n0)

    return mink, maxk


def coefs(nviz: int, rules: list):
    """
    Returns an array containing the coefficients for computing dN(s,k).

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param rules: list (Satisfaction rules)

    :return N: np.array
    """

    Nu = zera(nviz)[:]
    Nw = zera(nviz)[:]
    Nx = zera(nviz)[:]
    Ny = zera(nviz)[:]
    Nz = zera(nviz)[:]
    V = nviz - 1

    # Configurations of dissatisfied central dead and living cells:

    for n in range(0, 2 ** V):
        e0 = estados(nviz, n)[:]
        n1 = [0, 0]
        for i in range(1, nviz):
            n1[0] += e0[i]  # number of living neighbors around a dead central cell
        if sat(rules[0], n1[0]):  # skip if the dead central cell is satisfied
            continue
        mink0, maxk0 = minmaxvizvivos(nviz, e0, L=vizinhos(nviz)[0], v=vizinhos(nviz)[1])

        for m in range(0, 2 ** V):
            e1 = estados(nviz, m)[:]
            e1[0] = 1
            n1[1] = sum(e1[1:])
            if sat(rules[1], n1[1]):  # skip if the living central cell is satisfied
                continue
            mink1, maxk1 = minmaxvizvivos(nviz, e1, L=vizinhos(nviz)[0], v=vizinhos(nviz)[1])

            # u-coefficient updates
            Nu[0][n1[0]][n1[0]][n1[1]] -= 1
            Nu[1][n1[0]][n1[0]][n1[1]] += 1
            Nu[0][n1[1]][n1[0]][n1[1]] += 1
            Nu[1][n1[1]][n1[0]][n1[1]] -= 1

            # w and x coefficient updates
            for i in range(0, V):
                for k in range(mink0[i], maxk0[i] + 1):
                    Nw[e0[i + 1]][k + 1][n1[0]][n1[1]] += 1
                    Nx[e0[i + 1]][k][n1[0]][n1[1]] -= 1

            # y and z coefficient updates
            for i in range(0, V):
                for k in range(mink1[i], maxk1[i] + 1):
                    Ny[e1[i + 1]][k][n1[0]][n1[1]] -= 1
                    Nz[e1[i + 1]][k - 1][n1[0]][n1[1]] += 1

    # Normalize by binomial factors
    for s in range(0, 2):
        for k in range(0, nviz):
            for h in range(0, nviz):
                for l in range(0, nviz):
                    factor = float(comb(V, h) * comb(V, l))
                    Nu[s][k][h][l] /= factor
                    Nw[s][k][h][l] /= factor
                    Nx[s][k][h][l] /= factor
                    Ny[s][k][h][l] /= factor
                    Nz[s][k][h][l] /= factor

    return np.array([Nu, Nw, Nx, Ny, Nz], dtype=float)


"""---------------------------------- Functions for dynamic calculation  ---------------------------------"""


def condsiniciais(nviz: int, r1: float, f: np.array):
    """
    Sets initial conditions for f(s,k) using the binomial distribution.

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param r1: float (Fraction of type-1 cells in the system)
    :param f: np.array (Fractions f(s, k) – fractions of type-s cells with k living neighbors)

    :return f: np.array
    """
    r0 = 1 - r1
    for s in range(0, 2):
        for k in range(0, nviz):
            f[s, k] = comb(nviz - 1, k) * (r0 ** (nviz - 1 - k)) * (r1 ** k)

    return f


def total_ins(nviz: int, rules: list, dy: np.array):
    """
    Returns a list representing the total number of dissatisfied cells of each type.

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param rules: list (Satisfaction rules)
    :param dy: np.array (Fractions f(s, k))

    :return total: list
    """

    f = np.array([dy[0:5], dy[5:10]])
    total = [0, 1]
    for n in range(0, nviz):
        if sat(rules[1], n):
            total[1] -= f[1][n]
        if not sat(rules[0], n):
            total[0] += f[0][n]

    return total


def derivs(nviz: int, rules: list, dy: np.array, N: np.array, r0: float):
    """
    Returns dy/dt – the derivatives of all f(s, k).

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param rules: list (Satisfaction rules)
    :param dy: np.array (Current flattened f(s, k))
    :param N: np.array (Array of coefficients)
    :param r0: float (Fraction of type-0 cells)

    :return dy: np.array (Flattened derivatives)
    """
    V = nviz - 1
    u, w, x, y, z = N
    r1 = 1 - r0
    f = np.array([dy[0:​nviz], dy[nviz:2*​nviz]])

    # Coefficients for 0-neighbor transitions
    c0 = np.zeros(nviz, dtype=float)
    for j in range(1, nviz):
        c0[j] = comb(V - 1, j - 1) / comb(V, j)

    # Coefficients for max-neighbor transitions
    c8 = np.zeros(nviz, dtype=float)
    for k in range(0, V):
        c8[k] = comb(V - 1, k) / comb(V, k)

    normaliza0 = np.array([np.dot(c0, f[s]) for s in range(2)])
    normaliza8 = np.array([np.dot(c8, f[s]) for s in range(2)])

    f0 = np.vstack([np.append(c0 * f[s] / normaliza0[s], 0) for s in range(2)])
    f8 = np.vstack([np.append(c8 * f[s] / normaliza8[s], 0) for s in range(2)])

    total_ins_1 = total_ins(nviz, rules, dy)[1]
    total_ins_0 = total_ins(nviz, rules, dy)[0]
    p = f[1] / total_ins_1
    q = f[0] / total_ins_0

    c = np.array([
        r0 / (r0 * total_ins_0 + r1 * total_ins_1),
        r1 / (r0 * total_ins_0 + r1 * total_ins_1)
    ])

    df = np.zeros((2, nviz), dtype=float)

    for s in range(0, 2):
        for k in range(0, nviz):
            for h in range(0, nviz):
                for l in range(0, nviz):
                    df[s][k] += (
                        u[s][k][h][l]
                        + w[s][k][h][l] * f8[s][k - 1]
                        + y[s][k][h][l] * f0[s][k]
                        + x[s][k][h][l] * f8[s][k]
                        + z[s][k][h][l] * f0[s][k + 1]
                    ) * p[l] * q[h]

    dydx = np.zeros((2, nviz), dtype=float)
    for s in range(0, 2):
        dydx[s] = df[s] / c[s]

    return dydx.flatten().copy()


def execute_w(nviz, regras, N, dt, tmax, r1, f, writers, cp_find=False):
    """
    Executes the integration using SciPy’s solve_ivp, writing results to CSV.
    """
    writer0 = writers[0]
    writer1 = writers[1]

    dy = f.flatten().copy()

    def deriva(t, y):
        return derivs(nviz, regras, y, N, 1. - r1)

    def stop(t, y):
        total_ins_s = [total_ins(nviz, regras, y)[0], total_ins(nviz, regras, y)[1]]
        # Stop if any type becomes fully satisfied or any fraction drops negative
        return 0. if (total_ins_s[0] <= 0. or total_ins_s[1] <= 0. or np.any(y < 0)) else 1.

    stop.terminal = True

    yout = sp.integrate.solve_ivp(
        deriva,
        [0, tmax],
        dy,
        max_step=dt,
        method="DOP853",
        events=stop,
        dense_output=True
    )

    for i in range(0, len(yout.t)):
        f = np.array([yout.y.T[i][0:​nviz], yout.y.T[i][nviz:2*​nviz]])
        soma = f.sum(axis=1)
        total_ins_s = [total_ins(nviz, regras, yout.y.T[i])[0], total_ins(nviz, regras, yout.y.T[i])[1]]

        if not cp_find:
            writer0.writerow(row(nviz, total_ins_s, yout.t[i], f, soma, 0))
            writer1.writerow(row(nviz, total_ins_s, yout.t[i], f, soma, 1))

    return yout.t[-1]


def pc_finder(nviz, regras, N, dt, tmax, r_ativa, r_fr, writers):
    """
    Finds the critical fraction by bisection where the system’s behavior changes.
    """
    emax = abs(r_ativa - r_fr)
    r_b = (r_fr + r_ativa) / 2

    while emax >= 1E-8:
        f = np.empty((2, nviz), dtype=float)  # fractions f(s, k)
        soma = np.empty(2, dtype=float)        # total cells per type

        r1 = r_b
        r0 = 1 - r_b

        print(f"\n rule: {regras[1]} | r1: {r_b:.12f}\n")

        condsiniciais(nviz, r1, f)

        T = execute_w(nviz, regras, N, dt, tmax, r1, f, writers, cp_find=True)
        dy = f.flatten().copy()

        total_ins_s = [total_ins(nviz, regras, dy)[0], total_ins(nviz, regras, dy)[1]]

        atv = T >= tmax
        parada = "T = tmax" if atv else "pu = 0"

        if not atv:
            r_fr = r_b
        else:
            r_ativa = r_b

        print(
            f"\nT: {T:.10f} | Stop reason: {parada}" +
            f" | rule: {regras[1]} | r1: {r_b:.10f}" +
            f" | Ins 0: {total_ins_s[0]:.14f} | Ins 1: {total_ins_s[1]:.14f}\n"
        )

        emax = abs(r_ativa - r_fr)
        r_b = (r_fr + r_ativa) / 2

    print(f"\nr = {r_b:.12f} with max error = {emax:.12f}\n")

    return [r_b, emax]


"""---------------------------------- Auxiliary functions  ---------------------------------"""


def show(nviz: int, N: np.array):
    """
    Prints the calculated coefficient matrices to the console.
    """
    Nstr = ["u", "w", "x", "y", "z"]
    for h in range(0, nviz):
        for l in range(0, nviz):
            print(f"\nh = {h}, l = {l}\n")
            for i, c in enumerate(Nstr):
                print(c, *N[i][:, :, h, l])


def row(nviz: int, total_ins_s: list, t: int, f: np.array, soma: np.array, s: int, header=False):
    """
    Creates a CSV row based on the current state.

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param total_ins_s: list ([total_ins_0, total_ins_1])
    :param t: float (Time step)
    :param f: np.array (Fractions f(s, k))
    :param soma: np.array (Sum of fractions [sum_0, sum_1])
    :param s: int (Agent type 0 or 1)
    :param header: bool (True = header row; False = data row)

    :return: row list
    """
    if header:
        linha = ["t", f"total_ins_{s}"] + [f"f({s}, {m})" for m in range(nviz)] + [f"Sum({s})"]
    else:
        linha = [f"{t:.10f}", f"{total_ins_s[s]:.10f}"]
        linha += [f"{f[s][m]:.10f}" for m in range(nviz)]
        linha.append(f"{soma[s]:.10f}")
    return linha


def row2(nviz: int, N: np.array, coef: int, s: int, k: int, h: int, header=False):
    """
    Creates a CSV row for the coefficients.

    :param nviz: int (Number of cells in the region – Neighborhood + Central)
    :param N: np.array (Coefficient arrays)
    :param coef: int (Index of the desired coefficient)
    :param s: int (Agent type)
    :param k: int
    :param h: int
    :param header: bool (True = header row; False = data row)

    :return: row list
    """
    if header:
        return [f"l={m}" for m in range(nviz)]
    else:
        return [f"{N[coef][s][k][h][m]:.10f}" for m in range(nviz)]


def rules_config():
    """
    Selects the rules to run and creates output directories.
    """
    tipos = ["binomiais", "blinkers", "frozen", "mixeds", "segregadas"]

    all_rules1 = ['00000', '00001', '00010', '00011', '00100', '00101', '00110', '00111', '01000', '01001', '01010',
                  '01011', '01100', '01101', '01110', '01111', '10000', '10001', '10010', '10011', '10100', '10101',
                  '10110', '10111', '11000', '11001', '11010', '11011', '11100', '11101', '11110']

    rules_segreg1 = ["00011", "01011", "10011", "11000", "11010", "11001"]
    rules_blinkers1 = ["00110", "01100", "01110"]
    rules_mixeds1 = ["01101", "10110"]
    rules_frozen1 = ["00111", "01111", "10111", "11100", "11110", "11101"]
    rules_binomiais1 = [e for e in all_rules1 if
                        all(e not in lst for lst in [rules_frozen1, rules_mixeds1, rules_segreg1, rules_blinkers1])]

    all_rules0 = [e[::-1] for e in all_rules1]
    rules_segreg0 = [e[::-1] for e in rules_segreg1]
    rules_blinkers0 = [e[::-1] for e in rules_blinkers1]
    rules_mixeds0 = [e[::-1] for e in rules_mixeds1]
    rules_frozen0 = [e[::-1] for e in rules_frozen1]
    rules_binomiais0 = [e[::-1] for e in rules_binomiais1]

    tipos_lista1 = [rules_binomiais1, rules_blinkers1, rules_frozen1, rules_mixeds1, rules_segreg1]

    for e in tipos:
        try:
            mkdir(f"./Output/ScriptOutput_TransitionRules/{e}")
            mkdir(f"./Output/ScriptOutput_TransitionRules/{e}/flags")
        except FileExistsError:
            pass

    rules0 = []
    rules1 = []
    j = 0

    # Selection of rules to execute:
    while True:
        r = input("Enter rule: ")

        if r == "all":
            rules0 = all_rules0[:]
            rules1 = all_rules1[:]
            break
        elif r == "segregadas":
            rules0 += rules_segreg0
            rules1 += rules_segreg1
        elif r == "blinkers":
            rules0 += rules_blinkers0
            rules1 += rules_blinkers1
        elif r == "mixeds":
            rules0 += rules_mixeds0
            rules1 += rules_mixeds1
        elif r == "frozen":
            rules0 += rules_frozen0
            rules1 += rules_frozen1
        elif r == "binomiais":
            rules0 += rules_binomiais0
            rules1 += rules_binomiais1
        elif r == "":
            break
        else:
            rules1.append(r)
            rules0.append(all_rules0[all_rules1.index(r)])

    regras = [rules0, rules1]

    return [regras, tipos_lista1]
