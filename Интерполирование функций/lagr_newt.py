from math import asin, pi, cos
import matplotlib.pyplot as plt
from prettytable import PrettyTable

a = -0.8
b = 1.2
num = [2, 9, 29]


def f(args):
    ans = []
    for i in range(len(args)):
        y = args[i] ** 2 - asin(args[i] - 0.2)
        ans.append(y)
    return ans


def ravn_points(a, b, n):
    ans = []
    for i in range(n + 1):
        x_i = a + i * (b - a) / n
        ans.append(x_i)
    return ans



def optim_points(a, b, n):
    ans = []
    for i in range(n + 1):
        x_i = ((b - a) * cos(pi * (2 * i + 1) / (2 * (n + 1))) + (b + a)) / 2
        ans.append(x_i)
    ans.reverse()
    return ans


def l_k(x, arg, k):
    l_k = 1
    n = len(arg)
    for i in range(n):
        if i != k:
            l_k *= (x - arg[i]) / (arg[k] - arg[i])
    return l_k


def LIP(x, arg, val):  # x - индекс
    n = len(arg)
    ans = 0
    for i in range(n):
        ans += l_k(x, arg, i) * val[i]

    return ans


def n_k(arg, n):
    k = 0
    val = f(arg)
    for i in range(n):
        temp = val[i]
        for j in range(n):
            if j != i:
                temp /= arg[i] - arg[j]
        k += temp

    return k


def n_ks(arg):  # возвращает массив коэффициентов
    val = f(arg)

    ks = [val[0]]
    for i in range(2, len(arg) + 1):
        ks.append(n_k(arg[:i], i))
    return ks


def NIP(x, arg, ks):
    n = len(arg)
    ans = 0

    for i in range(n):
        temp = ks[i]
        for j in range(i):
            temp *= x - arg[j]
        ans += temp
    return ans


def otklon(val_, arg):
    maxx = 0
    y = f(arg)

    for i in range(len(val_)):
        temp = abs(val_[i] - y[i])
        if temp > maxx:
            maxx = temp

    return maxx


table_L = PrettyTable()
table_L.field_names = ["Количество узлов (n)", "Количество проверочных точек (m)", "Максимальное отклонение (RLn)",
                       "Максимальное отклонение (RLoptn)", "Максимальное отклонение (RNn)",
                       "Максимальное отклонение (RNoptn)"]

plt.figure(figsize=(12, 6))
for i in range(len(num)):
    n = num[i]

    arg = ravn_points(a, b, n)
    val = f(arg)

    arg_opt = optim_points(a, b, n)
    val_opt = f(arg_opt)

    for m in [n * 3, n * 6, n * 9]:  # [n * 5, n * 10, n * 15]
        arg_m = ravn_points(a, b, m)

        L = []  # хранятся новые значения функции(y) для m точек
        N = []

        L_opt = []
        N_opt = []
        for j in range(m):
            L.append(LIP(arg_m[j], arg, val))
            ks = n_ks(arg)  # был rg_m
            N.append(NIP(arg_m[j], arg, ks))

            L_opt.append(LIP(arg_m[j], arg_opt, val_opt))
            ks = n_ks(arg_opt)
            N_opt.append((NIP(arg_m[j], arg_opt, ks)))

        RLn = otklon(L, arg_m)
        RNn = otklon(N, arg_m)

        RLopt = otklon(L_opt, arg_m)
        RNopt = otklon(N_opt, arg_m)
        # print("L_opt = ", L_opt)

        table_L.add_row([n + 1, m, RLn, RLopt, RNn, RNopt])

    # строим графики
    k = 200
    x = ravn_points(a, b, k)
    y = f(x)
    L = []
    N = []

    L_opt = []
    N_opt = []
    for j in range(k + 1):
        L.append(LIP(x[j], arg, val))
        ks = n_ks(arg)  # был rg_m
        N.append(NIP(x[j], arg, ks))

        L_opt.append(LIP(x[j], arg_opt, val_opt))
        ks = n_ks(arg_opt)
        N_opt.append((NIP(x[j], arg_opt, ks)))
    # print(len(x))
    # print(len(L))

    plt.subplot(231 + i)
    plt.grid(True)
    plt.plot(x, L, label="L")
    plt.plot(x, L_opt, label="L_opt")
    plt.plot(x, y, label="f(x)")
    plt.text(-0.8, -0.2, f'Узлов: {n + 1}')
    plt.legend()

    plt.subplot(234 + i)
    plt.grid(True)
    plt.plot(x, N, label="N")
    plt.plot(x, N_opt, label="N_opt")
    plt.plot(x, y, label="f(x)")
    plt.text(-0.8, -0.2, f'Узлов: {n + 1}')
    plt.legend()

print(table_L)

plt.show()
