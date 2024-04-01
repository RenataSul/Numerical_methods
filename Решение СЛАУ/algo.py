import numpy as np


def sqrt(x):
    if x == 0:
        return 0
    eps = 0.000001
    a = x / 2
    while 1:
        answ = 1 / 2 * (a + x / a)
        if abs(answ - a) < eps:
            return answ
        a = answ


def E(n):
    m = []
    for i in range(n):
        r = []
        for j in range(n):
            if i == j:
                r.append(1)
            else:
                r.append(0)
        m.append(r)
    return Matrix(m)


def O(n):
    m = []
    for i in range(n):
        r = []
        for j in range(n):
            r.append(0)
        m.append(r)
    return Matrix(m)


class Matrix:
    def __init__(self, m):
        self.data = [[m[i][j] for j in range(len(m[0]))] for i in range(len(m))]
        self.rows = len(m)
        self.cols = len(m[0])

    def __str__(self):
        matrix_str = ""
        for row in self.data:
            for element in row:
                matrix_str += str(element) + " "
            matrix_str += "\n"
        return matrix_str

    def __getitem__(self, index):
        if index[1] is None:  # [int, None]
            return Matrix([self.data[index[0]]])
        elif index[0] is None:  # [None, int]
            return Matrix([[self.data[i][index[1]]] for i in range(self.rows)])

        return self.data[index[0]][index[1]]

    def __setitem__(self, index, value):
        if index[1] is None:  # [int, None]
            for j in range(self.cols):
                self.data[index[0]][j] = value[0, j]
        elif index[0] is None:  # [None, int]
            for i in range(self.rows):
                self.data[i][index[1]] = value[i]
        else:
            self.data[index[0]][index[1]] = value

    def shape(self):
        s = []
        s.append(self.rows)
        s.append(self.cols)

        return s

    def transpose(self):
        transposed_matrix = []
        for i in range(self.cols):
            row = []
            for j in range(self.rows):
                row.append(self.data[j][i])
            transposed_matrix.append(row)
        return Matrix(transposed_matrix)

    def __add__(self, other):
        result = []
        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                el = self.data[i][j] + other.data[i][j]
                row.append(el)
            result.append(row)
        return Matrix(result)

    def __sub__(self, other):
        return self + -1 * other

    def __mul__(self, other):
        result = []
        if isinstance(other, Matrix):
            for i in range(self.rows):
                row = []
                for j in range(other.cols):
                    el = 0
                    for k in range(self.cols):
                        el += self.data[i][k] * other.data[k][j]
                    row.append(el)
                result.append(row)
            return Matrix(result)

        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                el = self.data[i][j] * other
                row.append(el)
            result.append(row)
        return Matrix(result)

    def __rmul__(self, other):
        return self.__mul__(other)

    # def __abs__(self):
    #    result = []
    #    for i in range(self.rows):
    #        row = []
    #        for j in range(self.cols):
    #            el = abs(self.data[i][j])
    #            row.append(el)
    #        result.append(row)
    #    return Matrix(result)

    def norm(self, size_of_norm):
        norm = 0
        if size_of_norm == 1:
            for j in range(self.cols):
                summ = 0
                for i in range(self.rows):
                    summ += abs(self.data[i][j])
                norm = max(norm, summ)

        elif size_of_norm == "inf":
            for i in range(self.rows):
                summ = 0
                for j in range(self.cols):
                    summ += abs(self.data[i][j])
                norm = max(norm, summ)

        elif size_of_norm == 2:
            summ = 0
            for i in range(self.rows):
                for j in range(self.cols):
                    summ += abs(self.data[i][j]) ** 2
            norm = sqrt(summ)

        return norm

    def copy(self):
        result = []

        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                el = self.data[i][j]
                row.append(el)
            result.append(row)
        return Matrix(result)


def swap1(M: Matrix, max_i, i):
    temp = [0 for i in range(M.cols)]  # хранит i-ую строку
    for k in range(M.cols):
        temp[k] = (M.data[i][k])

    for k in range(M.cols):
        M[i, k] = M[max_i, k]
        M[max_i, k] = temp[k]


def LU(A, b, solution):
    b = Matrix(b)
    M = Matrix(A)
    n = M.rows
    P = E(n)

    # 1.1
    for i in range(n):

        max_el = 0
        max_i = i

        for j in range(i, n):
            if abs((M.data[j][i])) > max_el:
                max_el = abs((M.data[j][i]))
                max_i = j  # строка с ведущим элементом в i-ом столбце
        #  меняем i-ую строку с j-ой строкой в матрице М
        swap1(M, max_i, i)
        swap1(P, max_i, i)

        #  Преобразование матрицы M
        for j in range(i + 1, n):
            M[j, i] = M[j, i] / M[i, i]
            for k in range(i + 1, n):
                M[j, k] = M[j, k] - M[j, i] * M[i, k]

    # M = L + U - E

    L = E(n)
    for i in range(1, n):
        for j in range(i):
            L[i, j] = M[i, j]

    U = O(n)
    for i in range(n):
        for j in range(i, n):
            U[i, j] = M[i, j]

    #  Решение
    Pb = P * b  # после перестановки
    x = []
    y = []
    for i in range(n):
        row = []
        row.append(0)
        x.append(row)
        y.append(row)

    x = Matrix(x)
    y = Matrix(y)

    for k in range(n):
        y[k, 0] = Pb[k, 0]
        for i in range(k):
            y[k, 0] -= L[k, i] * y[i, 0]

    for k in reversed(range(n)):
        x[k, 0] = y[k, 0] / U[k, k]
        for i in range(k + 1, n):
            x[k, 0] -= U[k, i] * x[i, 0] / U[k, k]


    solution = Matrix(solution)
    answ_x = []
    answ_pogr = []
    for i in range(n):
        row1 = []
        row2 = []

        row1.append(x[i, 0])
        row2.append((solution - x)[i, 0])
        answ_x.append(row1)
        answ_pogr.append(row2)
    print("РЕШЕНИЕ LU:\t\t\t\tПогрешность LU")
    for t in zip(answ_x, answ_pogr):
        print(t[0], "\t", t[1])


def QR(A, b, solution):
    b = Matrix(b)
    M = Matrix(A)
    n = M.rows
    En = E(n)
    Q = E(n)
    R = M

    for i in range(n - 1):
        y = Matrix([[R[j, i]] for j in range(i, n)])
        z = Matrix([[En[j, i]] for j in range(i, n)])

        a = y[None, 0].norm(2)
        p = (y - z * a).norm(2)
        w = (y - z * a) * (1 / p)

        H = E(n)
        H_ = (E(n - i) - 2 * w * w.transpose())
        for j in range(i, n):
            for k in range(i, n):
                H[j, k] = H_[j - i, k - i]

        Q = Q * H.transpose()
        R = H * R

    y = Q.transpose() * b
    x = []
    for i in range(n):
        row = []
        row.append(0)
        x.append(row)
    x = Matrix(x)

    for k in reversed(range(n)):
        x[k, 0] = y[k, 0] / R[k, k]
        for i in range(k + 1, n):
            x[k, 0] = x[k, 0] - R[k, i] * x[i, 0] / R[k, k]

    solution = Matrix(solution)
    answ_x = []
    answ_pogr = []
    for i in range(n):
        row1 = []
        row2 = []

        row1.append(x[i, 0])
        row2.append((solution - x)[i, 0])
        answ_x.append(row1)
        answ_pogr.append(row2)
    print("РЕШЕНИЕ QR:\t\t\t\tПогрешность QR")
    for t in zip(answ_x, answ_pogr):
        print(t[0], "\t", t[1])


def Zeidel(A, b, solution):
    A = Matrix(A)
    b = Matrix(b)
    n = A.cols
    eps = 0.00001
    # print("A: \n", A)
    mu = 1 / A.norm("inf")
    B = E(n) - mu * A

    if B.norm("inf") >= 1:
        T = A.transpose()
        A, b = T * A, T * b

        mu = 1 / A.norm("inf")
        B = E(n) - mu * A
    c = mu * b
    x = c.copy()
    k = 0

    #  проверка на сходимость
    maxx = 0
    for i in range(n):
        maxx = A[i, i]
        summ = 0
        for j in range(n):
            el = A.data[i][j]
            if j != i:
                summ += abs(el)
        # print(summ)
        if maxx < summ:
            b = A.transpose() * b
            A = A.transpose() * A
            break

    # строим С и d
    B = E(n)
    C = []
    for i in range(n):
        row = []
        row.append(0)
        C.append(row)
    C = Matrix(C)

    for i in range(n):
        C[i, 0] = b[i, 0] / A[i, i]
        for j in range(n):
            if i != j:
                B[i, j] = -1 * A[i, j] / A[i, i]
            else:
                B[i, j] = 0

    x = c.copy()
    k = 0

    while 1:
        x_ = x.copy()
        k += 1

        for i in range(n):
            sum_row = 0
            for j in range(n):
                if i != j:
                    sum_row += B[i, j] * x_[j, 0]
                else:
                    sum_row += C[i, 0]
            x_[i, 0] = sum_row

        if abs((A * x - b).norm(1)) <= eps:
            x = x_.copy()
            break
        x = x_.copy()

    solution = Matrix(solution)

    answ_x = []
    answ_pogr = []
    for i in range(n):
        row1 = []
        row2 = []

        row1.append(x[i, 0])
        row2.append((solution - x)[i, 0])
        answ_x.append(row1)
        answ_pogr.append(row2)
    print("РЕШЕНИЕ Зейдель:\t\t\tПогрешность Зейдель")
    for t in zip(answ_x, answ_pogr):
        print(t[0], "\t", t[1])


def MPI(A, b, solution):
    A1 = Matrix(A).copy()
    b = Matrix(b).copy()
    n = A1.cols
    eps = 0.0001

    m = 1 / A1.norm("inf")
    B = E(n) - m * A1

    if B.norm("inf") >= 1:
        b = A1.transpose() * b
        A1 = A1.transpose() * A1
        m = 1 / A1.norm("inf")
        B = E(n) - m * A1

    c = m * b

    x = c.copy()

    k = 0

    while 1:
        x_ = B * x + c
        k += 1

        if abs(B.norm("inf") * (1 / (1 - B.norm("inf"))) * (x_ - x).norm("inf")) < eps:
            x = x_.copy()
            break

        x = x_.copy()


    solution = Matrix(solution)

    answ_x = []
    answ_pogr = []
    for i in range(n):
        row1 = []
        row2 = []

        row1.append(x[i, 0])
        row2.append((solution - x)[i, 0])
        answ_x.append(row1)
        answ_pogr.append(row2)
    print("РЕШЕНИЕ МПИ:\t\t\tПогрешность МПИ")
    for t in zip(answ_x, answ_pogr):
        print(t[0], "\t", t[1])


def Print(A, b, solution):
    LU(A, b, solution)
    print("")
    QR(A, b, solution)
    print("")
    MPI(A, b, solution)
    print("")
    Zeidel(A, b, solution)
    print("")


def test0():
    print("\nTEST 0")
    A = [[0, 2, 3], [1, 2, 4], [4, 5, 6]]
    b = [[13], [17], [32]]
    solution = [[1], [2], [3]]
    Print(A, b, solution)


def test1():
    print("ТЕСТ 1")
    A = [[11, 1, 1], [1, 13, 1], [1, 1, 15]]
    b = [[13], [15], [17]]
    M = np.array(A)
    v = np.array(b)
    solution = np.linalg.solve(M, v)
    # solution = [[1], [1], [1]]
    Print(A, b, solution)


def test2():
    print("ТЕСТ 2")
    A = [[-11, 1, 1], [1, -13, 1], [1, 1, -15]]
    b = [[-13], [-15], [-17]]
    M = np.array(A)
    v = np.array(b)
    solution = np.linalg.solve(M, v)
    # solution = [[375 / 263], [359 / 263], [347 / 263]]
    Print(A, b, solution)


def test3():
    print("ТЕСТ 3")
    A = [[-11, 12, 13], [14, -13, 10], [13, 14, -15]]
    b = [[13], [15], [17]]
    M = np.array(A)
    v = np.array(b)
    solution = np.linalg.solve(M, v)
    # solution = [[601 / 411], [2479 / 2055], [1928 / 2055]]
    Print(A, b, solution)


def test4():
    print("ТЕСТ 4")
    A = [[11, 10, 10], [10, 13, 10], [10, 10, 15]]
    b = [[13], [15], [17]]
    M = np.array(A)
    v = np.array(b)
    solution = np.linalg.solve(M, v)
    # solution = [[-5 / 49], [31 / 49], [191 / 245]]
    Print(A, b, solution)


def test5(n, eps):
    print("ТЕСТ 5")

    A = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(1)
            if i < j:
                row.append(-1)
            if i > j:
                row.append(0)
        A.append(row)


    B = []
    for i in range(n):
        row = []
        for j in range(n):
            if i <= j:
                row.append(1)

            if i > j:
                row.append(-1)
        B.append(row)

    b = []
    for i in range(n):
        row = []
        if i == (n - 1):
            row.append(1)
        else:
            row.append(-1)
        b.append(row)


    for i in range(n):
        for j in range(n):
            A[i][j] = A[i][j] + eps * B[i][j]

    M = np.array(A)
    v = np.array(b)
    solution = np.linalg.solve(M, v)

    Print(A, b, solution)


test0()
print("--------------------------")
test1()
print("--------------------------")
test2()
print("--------------------------")
test3()
print("--------------------------")
test4()
print("--------------------------")
test5(3, 0.0001)
