import numpy as np
from math import asin, pi, cos
from prettytable import PrettyTable
import matplotlib.pyplot as plt

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

# def optim_points(a, b, n):
#     ans = []
#     for i in range(n + 1):
#         if i == 0:
#             ans.append(b)
#         if i == n:
#             ans.append(a)
#         else:
#             x_i = ((b - a) * cos(pi * (2 * i + 1) / (2 * (n + 1))) + (b + a)) / 2
#             ans.append(x_i)
#     ans.reverse()
#     return ans


def otklon(val_, arg):
    maxx = 0
    y = f(arg)

    for i in range(len(val_)):
        temp = abs(val_[i] - y[i])
        if temp > maxx:
            maxx = temp

    return maxx


# СПЛАЙНЫ
# ЛИНЕЙНЫЙ СПЛАЙН
def k_1(arg, val):
    n = len(arg)
    x = np.eye(2)
    y = np.zeros((2, 1))

    ks = np.zeros((2 * (n - 1), 1))  # вектор коэффициентов
    for i in range(n - 1):
        x[0, 0] = arg[i]
        x[1, 0] = arg[i + 1]
        x[0, 1] = 1
        y[0, 0] = val[i]
        y[1, 0] = val[i + 1]
        solve = np.linalg.solve(x, y)
        ks[2 * i, 0] = solve[0, 0]
        ks[2 * i + 1] = solve[1, 0]

    return ks


def spline_1(x, arg, ks):
    n = len(arg)
    # print("x: ", x)

    for i in range(n - 1):
        if arg[i] <= x <= arg[i + 1]:
            return ks[2 * i, 0] * x + ks[2 * i + 1, 0]

    # for i in range(n - 1):
    #     if arg[i] <= x <= arg[i + 1]:
    #         print("great")
    #         return ks[2 * i, 0] * x + ks[2 * i + 1, 0]
    #     else:
    #         print("arg[i]: ", arg[i])
    #         print("arg[i+1]: ", arg[i+1])
    #         break


# КВАДРАТИЧНЫЙ СПЛАЙН
def k_2(arg, val):
    n = len(arg)

    x = np.zeros((3 * (n - 1), 3 * (n - 1)))
    y = np.zeros((3 * (n - 1), 1))

    y[-1] = 0  # естественный сплайн
    for i in range(n - 1):
        x[3 * i, 3 * i] = arg[i] ** 2
        x[3 * i, 3 * i + 1] = arg[i]
        x[3 * i, 3 * i + 2] = 1
        x[3 * i + 1, 3 * i] = arg[i + 1] ** 2
        x[3 * i + 1, 3 * i + 1] = arg[i + 1]
        x[3 * i + 1, 3 * i + 2] = 1
        x[3 * i + 2, 3 * i] = 2 * arg[i + 1]
        x[3 * i + 2, 3 * i + 1] = 1

        if i != n - 2:
            x[3 * i + 2, 3 * i + 3] = -2 * arg[i + 1]
            x[3 * i + 2, 3 * i + 4] = -1

        y[3 * i] = val[i]
        y[3 * i + 1] = val[i + 1]

    ks = np.linalg.solve(x, y)
    return ks


def spline_2(x, arg, ks):
    n = len(arg)

    for i in range(n - 1):
        if arg[i] <= x <= arg[i + 1]:
            return ks[3 * i, 0] * x ** 2 + ks[3 * i + 1, 0] * x + ks[3 * i + 2, 0]


# КУБИЧЕСКИЙ СПЛАЙН
def k_3(arg, val):
    n = len(arg)
    d = [0, 0]
    h_arg = []
    h_val = []
    for i in range(n - 1):
        h_arg.append(arg[i + 1] - arg[i])
        h_val.append(val[i + 1] - val[i])

    H = np.zeros((n - 2, n - 2))
    gamma = np.zeros((n - 2, 1))

    for i in range(n - 2):
        H[i, i] = 2 * (h_arg[i + 1] + h_arg[i])
        if i > 0:
            H[i, i - 1] = h_arg[i]
            H[i - 1, i] = h_arg[i]

        gamma[i] = 6 * (h_val[i + 1] / h_arg[i + 1] - h_val[i] / h_arg[i])

    solve = np.linalg.solve(H, gamma)

    val_pr2 = []  # при вторых производных
    val_pr2.append(d[0])
    for i in range(n - 2):
        val_pr2.append(solve[i])
    val_pr2.append(d[1])

    val_pr1 = []  # при первых произв
    for i in range(n - 1):
        val_pr1.append(h_val[i] / h_arg[i] - val_pr2[i + 1] * h_arg[i] / 6 - val_pr2[i] * h_arg[i] / 3)
    ks = np.zeros((4 * (n - 1), 1))
    for i in range(n - 1):
        j = 4 * i
        ks[j] = val[i]
        ks[j + 1] = val_pr1[i]
        ks[j + 2] = val_pr2[i] / 2
        ks[j + 3] = (val_pr2[i + 1] - val_pr2[i]) / (6 * h_arg[i])

    return ks


def spline_3(x, arg, ks):
    n = len(arg)
    for i in range(n - 1):
        if arg[i] <= x <= arg[i + 1]:
            j = 4 * i
            a = ks[j] + ks[j + 1] * (x - arg[i]) + ks[j + 2] * (x - arg[i]) ** 2 + ks[j + 3] * (x - arg[i]) ** 3
            return a[0]


table = PrettyTable()
table.field_names = ["Количество узлов (n)", "Количество проверочных точек (m)", "Максимальное отклонение (RS1)",
                     "Максимальное отклонение (RS2)", "Максимальное отклонение (RS3)"]

for j in range(len(num)):
    n = num[j]

    arg = ravn_points(a, b, n)
    val = f(arg)

    arg_opt = optim_points(a, b, n)
    val_opt = f(arg_opt)

    for m in [n * 3, n * 6, n * 9]:
        arg_m = ravn_points(a, b, m)

        S1 = []
        S2 = []
        S3 = []

        S1opt = []
        S2opt = []
        S3opt = []
        k1 = k_1(arg, val)
        k2 = k_2(arg, val)
        k3 = k_3(arg, val)
        for i in range(m + 1):
            # при равномерном
            S1.append(spline_1(arg_m[i], arg, k1))
            S2.append(spline_2(arg_m[i], arg, k2))
            S3.append(spline_3(arg_m[i], arg, k3))

        RS1 = otklon(S1, arg_m)
        RS2 = otklon(S2, arg_m)
        RS3 = otklon(S3, arg_m)

        table.add_row([n + 1, m, RS1, RS2, RS3])

print(table)

###ГРАФИКИ
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


def otklon_list(val_, arg):
    y = f(arg)
    temp = []
    for i in range(len(val_)):
        temp.append(abs(val_[i] - y[i]))

    return temp



plt.figure(figsize=(12, 6))
for i in range(len(num)):
    k = 200
    n = num[i]

    arg = ravn_points(a, b, n)
    val = f(arg)
    arg_opt = optim_points(a, b, n)
    val_opt = f(arg_opt)

    x = ravn_points(a, b, k)
    y = f(x)

    L = []
    L_opt = []
    S3 = []

    k3 = k_3(arg, val)
    for j in range(k + 1):
        L.append(LIP(x[j], arg, val))
        S3.append(spline_3(x[i], arg, k3))

    RL = otklon_list(L, x)
    RS3 = otklon_list(S3, x)

    # 321
    plt.subplot(131 + i)
    plt.grid(True)
    plt.plot(x, RL, label='RL')
    plt.plot(x, RS3, label='RS3')
    plt.text(0.75, 1, f'Узлов: {n + 1}')
    plt.legend()

plt.subplot(131)
plt.show()
