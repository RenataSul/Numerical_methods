import math
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


class Matrix(object):
    discard_index = 3  # Количество цифр после запятой

    def __new__(cls, matrix):
        """Объект матрицы можно создать, только если аргументом является либо файлом с матрицей,
        либо двумерным массивом"""
        try:
            if isinstance(matrix, list) and isinstance(matrix[0], list):
                return super().__new__(cls)
        except IndexError:
            pass

    def __init__(self, matrix):
        self.__matrix = matrix

        self.__r_count = len(self.__matrix)
        self.__c_count = len(self.__matrix[0])

    @property
    def dim(self):
        """Возвращает размерность матрицы в виде tuple (A.dim = (n, m))"""
        return self.__r_count, self.__c_count

    def __str__(self):
        """Вывод матрицы (print(A))"""
        str_matrix = []

        for row in self.__matrix:
            str_row = []
            for elem in row:
                format_elem = f"{elem:.{self.discard_index}f}" \
                    if f"{elem:.{self.discard_index}f}" != f'-{0:.{self.discard_index}f}' \
                    else f'{0:.{self.discard_index}f}'

                str_row.append(format_elem)
            str_matrix.append('   '.join(str_row))

        str_matrix = '\n'.join(str_matrix)

        return str_matrix

    def __eq__(self, other):
        """Проверка матриц на равенство (A == B)"""
        for i in range(self.__r_count):
            for j in range(self.__c_count):
                if self.__matrix[i][j] != other.__matrix[i][j]:
                    return False

        return True

    def __ne__(self, other):
        """Проверка матриц на неравенство (A != B)"""
        for i in range(self.__r_count):
            for j in range(self.__c_count):
                if self.__matrix[i][j] != other.__matrix[i][j]:
                    return True

        return False

    def __bool__(self):
        """Функция возвращает True, если матрица ненулевая"""
        for i in range(self.__r_count):
            for j in range(self.__c_count):
                if self.__matrix[i][j] != 0:
                    return True

        return False

    @property
    def is_pd(self):
        """Проверка матрицы на положительную определенность критерием Сильвестра"""
        for i in range(self.__r_count):
            minor = self.copy
            for _ in range(i + 1, self.__c_count):
                m = minor.dim[0] - 1
                minor = minor[m: m]
            elem = minor.det // abs(minor.det) if minor.det != 0 else 0

            if elem <= 0:
                return False

        return True

    @property
    def is_sym(self):
        """Проверка матрицы на симметричность"""
        for i in range(self.__r_count):
            for j in range(i + 1, self.__c_count):
                if self.__matrix[i][j] != self.__matrix[j][i]:
                    return False

        return True

    @property
    def is_dd(self):
        """Проверка матрицы на диагональное преобладание (только для квадратных матрицы)"""
        for i in range(self.__r_count):
            row_sum = 0
            for j in range(self.__c_count):
                row_sum = row_sum + abs(self.__matrix[i][j]) if i != j else row_sum

            if abs(self.__matrix[i][i]) <= row_sum:
                return False

        return True

    def __setitem__(self, index, value):
        """Функция задает значение элементу матрицы (A[i,j] = α),
        либо строке(столбцу) (A(i, None) = a, A(None, j) = a"""

        if isinstance(index, tuple):
            row, column = index

            if None not in index:
                self.__matrix[row][column] = value

            elif column is None:
                for i in range(self.__c_count):
                    self.__matrix[row][i] = value.__matrix[0][i]

            elif row is None:
                for i in range(self.__r_count):
                    self.__matrix[i][column] = value.__matrix[i][0]

    def __getitem__(self, index):
        """Функция возвращает либо обрезанную матрицу (A[i:j], например,
        дана матрица A = [[1, 0], [0, 1]], тогда A[1: 1] = [[1]], к тому же,
        вырезать только строку или столбец можно таким образом: A[i: None],
        A[None: j]), либо один элемент (A[i,j]), либо вектор из матрицы
        (A[i, None] или A[None, i])"""
        if isinstance(index, tuple):
            row, column = index

            if None not in index:
                return self.__matrix[row][column]

            if row is None:
                matrix = []
                for i in range(self.__r_count):
                    matrix.append([self.__matrix[i][column]])

                return Matrix(matrix)

            elif column is None:
                matrix = []

                row_list = []
                for i in range(self.__c_count):
                    row_list.append(self.__matrix[row][i])
                matrix.append(row_list)

                return Matrix(matrix)

        elif isinstance(index, slice):
            row_index = index.start
            column_index = index.stop

            matrix = []

            for i in range(self.__r_count):
                if i != row_index:
                    row = []
                    for j in range(self.__c_count):
                        if j != column_index:
                            row.append(self.__matrix[i][j])
                    matrix.append(row)

            return Matrix(matrix)

    def __mul__(self, other):
        """Функция умножает матрицу либо на число (A * α), либо на матрицу (A * B)"""
        if isinstance(other, float) or isinstance(other, int):
            matrix = []

            for i in range(self.__r_count):
                row = []
                for j in range(self.__c_count):
                    elem = self.__matrix[i][j] * other
                    row.append(elem)
                matrix.append(row)

            return Matrix(matrix)

        else:
            matrix = []

            for i in range(self.__r_count):
                row = []
                for j in range(other.__c_count):
                    elem = 0
                    for k in range(self.__c_count):
                        elem += self.__matrix[i][k] * other.__matrix[k][j]
                    row.append(elem)
                matrix.append(row)

            return Matrix(matrix)

    def __rmul__(self, number):
        """Функция позволяет коммутативно умножать матрицу на число"""
        return self * number

    def __imul__(self, number):
        """Перегрузка оператора *= для чисел"""
        return self * number

    def __truediv__(self, number):
        """Деление матрицы на число"""
        matrix = []

        for i in range(self.__r_count):
            row = []
            for j in range(self.__c_count):
                elem = self.__matrix[i][j] / number
                row.append(elem)
            matrix.append(row)

        return Matrix(matrix)

    def __idiv__(self, number):
        """Перегрузка оператора /= для чисел"""
        return self / number

    def __matmul__(self, other):
        """Скалярное произведение (a @ b)"""
        vec1 = ~self if self.__c_count == 1 else self
        vec2 = ~other if other.__r_count == 1 else other

        in_prod = 0

        for i in range(vec1.__c_count):
            in_prod += vec1.__matrix[0][i] * vec2.__matrix[i][0]

        return in_prod

    def __add__(self, other):
        """Сложение матриц (A + B)"""
        matrix = []

        for i in range(self.__r_count):
            row = []
            for j in range(self.__c_count):
                elem = self.__matrix[i][j] + other.__matrix[i][j]
                row.append(elem)
            matrix.append(row)

        return Matrix(matrix)

    def __iadd__(self, other):
        """Перегрузка оператора +="""
        return self + other

    def __sub__(self, other):
        """Разность матриц (A - B)"""
        matrix = []

        for i in range(self.__r_count):
            row = []
            for j in range(self.__c_count):
                elem = self.__matrix[i][j] - other.__matrix[i][j]
                row.append(elem)
            matrix.append(row)

        return Matrix(matrix)

    def __isub__(self, other):
        """Перегрузка оператора -="""
        return self - other

    def __pow__(self, power, modulo=None):
        """Бинарное возведение матрицы в степень (A ** α)"""
        if power == 1:
            return self.copy
        elif power < -1:
            return (self ** -power) ** -1
        elif power == 0:
            return Matrix.E(self.__r_count)
        elif power == 2:
            return self * self
        elif power == -1:  # Вычисляем обратную матрицу методом Гаусса
            n = self.dim[0]
            A = self.copy
            inv_A = Matrix.E(n)
            for i in range(self.__r_count):
                if A[i, i] == 0:  # Если число в диагонали равно нулю, тогда делаем перестановку
                    max_i = i
                    for j in range(i + 1, self.__r_count):  # Ищем индекс строки, где число в данном столбце не равно 0
                        if A[j, i] != 0:
                            max_i = j
                            break
                    for M in [A, inv_A]:  # Меняем строки в матрицах
                        add_vector = M[max_i, None]
                        M[max_i, None] = M[i, None]
                        M[i, None] = add_vector

                # Преобразуем обратную матрицу
                for j in range(self.__r_count):
                    if j != i:
                        factor = -A[j, i] / A[i, i]
                        for k in range(self.__c_count):
                            inv_A[j, k] += factor * inv_A[i, k]
                for k in range(self.__c_count):
                    inv_A[i, k] /= A[i, i]

                # Преобразуем исходную матрицу
                for j in range(self.__r_count):
                    if j != i:
                        factor = -A[j, i] / A[i, i]
                        for k in range(i + 1, self.__c_count):
                            A[j, k] += factor * A[i, k]
                for k in range(i + 1, self.__c_count):
                    A[i, k] /= A[i, i]
                for j in range(self.__r_count):
                    A[j, i] = 0 if j != i else 1

            return inv_A

        # Непосредственно биномиальное возведение в степень
        if power % 2 == 0:
            return (self ** 2) * (self ** (power // 2))
        else:
            return self * (self ** (power - 1))

    def __ipow__(self, power):
        """Перегрузка оператора **="""
        return self ** power

    def __invert__(self):
        """Транспонированние матрицы (~A)"""
        matrix = []

        for i in range(self.__c_count):
            row = []
            for j in range(self.__r_count):
                row.append(self.__matrix[j][i])
            matrix.append(row)

        return Matrix(matrix)

    @property
    def det(self):
        """Определитель матрицы, высчитываемый методом Гаусса (A.det)"""
        A = self.copy
        n = self.__r_count

        mult_num = 1  # На это число домножим итоговый ответ. Оно принимает значение либо 1, либо -1
        for i in range(n - 1):
            if A[i, i] == 0:  # Если число в диагонали равно нулю, тогда делаем перестановку
                max_i = i
                for j in range(i + 1, n):  # Ищем индекс строки, где число в данном столбце не равно 0
                    if A[j, i] != 0:
                        max_i = j
                        break
                # Меняем строки матрицы
                add_vector = A[max_i, None]
                A[max_i, None] = A[i, None]
                A[i, None] = add_vector

                mult_num *= -1  # При каждой перестановке домножаем на -1

            # Видоизменяем матрицу, зануляя колонку
            for j in range(i + 1, n):
                factor = -A[j, i] / A[i, i]
                for k in range(i + 1, n):
                    A[j, k] += A[i, k] * factor
            for j in range(i + 1, n):
                A[j, i] = 0

        d_elems = [A[i, i] for i in range(n)]  # Находим определитель, перемножив диагональные элементы матрицы
        det = mult_num * d_elems[0]
        for i in range(1, len(d_elems)):
            det *= d_elems[i]

        return det

    def __abs__(self):
        """Нахождение евклидовой нормы (abs(a))"""
        return sqrt(self @ self)

    def norm(self, norm_type):
        """Нахождение единичной нормы (A.norm(1)) и бесконечную норму (A.norm('inf'))"""
        if norm_type == 1:
            norm = 0

            for i in range(self.__c_count):
                column_sum = 0
                for j in range(self.__r_count):
                    column_sum += abs(self.__matrix[j][i])
                norm = column_sum if norm < column_sum else norm

            return norm

        if norm_type == 'inf':
            norm = 0

            for i in range(self.__r_count):
                row_sum = 0
                for j in range(self.__c_count):
                    row_sum += abs(self.__matrix[i][j])
                norm = row_sum if norm < row_sum else norm

            return norm

    @property
    def copy(self):
        """Копирует матрицу. Это необходимо, так как, например, для B = A ссылки на все элементы останутся теми
        же, то есть при изменении B будет изменяться и A"""
        matrix = []

        for i in range(self.__r_count):
            row = []
            for j in range(self.__c_count):
                row.append(self.__matrix[i][j])
            matrix.append(row)

        return Matrix(matrix)

    @staticmethod
    def E(n):
        """Возвращает единичную матрицу той же размерности, что и у матрицы"""
        matrix = []

        for i in range(n):
            row = []
            for j in range(n):
                elem = 0 if i != j else 1
                row.append(elem)
            matrix.append(row)

        return Matrix(matrix)

    @staticmethod
    def NULL_VECTOR(n):
        """Возвращает нулевой вектор"""
        return Matrix([[0] for _ in range(n)])


def LU(A: Matrix, b: Matrix) -> Matrix:
    n = A.dim[0]
    E = Matrix.E(n)

    M = A.copy
    P = E.copy

    for i in range(n - 1):
        # Определяем строку с максимальным по модулю элементом
        max_i = i

        for j in range(i, n):
            max_i = j if abs(M[max_i, i]) < abs(M[j, i]) else max_i

        # Выполняем перестановку в матрицах P и M, меняя строки местами
        for matrix in [P, M]:
            add_vector = matrix[max_i, None]
            matrix[max_i, None] = matrix[i, None]
            matrix[i, None] = add_vector

        # Преобразуем матрицу M
        for j in range(i + 1, n):
            M[j, i] = M[j, i] / M[i, i]
            for k in range(i + 1, n):
                M[j, k] = M[j, k] - M[j, i] * M[i, k]

    # Проводим LU-разложение
    L = E.copy
    U = E.copy

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = M[i, j]
            elif i < j:
                U[i, j] = M[i, j]
            else:
                L[i, i] = 1
                U[i, i] = M[i, i]

    # Решаем систему Ly = Pb и Ux = y
    Pb = P * b
    y = Matrix.NULL_VECTOR(n)
    for k in range(n):
        y[k, 0] = Pb[k, 0]
        for i in range(k):
            y[k, 0] -= L[k, i] * y[i, 0]

    x = Matrix.NULL_VECTOR(n)
    for k in reversed(range(n)):
        x[k, 0] = y[k, 0] / U[k, k]
        for i in range(k + 1, n):
            x[k, 0] -= U[k, i] * x[i, 0] / U[k, k]

    return x


# sin(x+0.5) - y = 1
# x+cos(y-2) = 2
def f(x, y, lmbd):
    return Matrix([[lmbd * math.sin(x + 0.5) - y - 1], [x + lmbd * math.cos(y - 2) - 2]])


def f_p(x, y):
    return Matrix([[math.cos(x + 0.5), -1], [1, -math.sin(y - 2)]])


def approx(st):
    # ищем начальное приближение
    n = 10
    val = st
    x = val[0, 0]
    y = val[1, 0]

    for i in range(1, 10 + 1):
        d_x = LU(f_p(x, y), f(x, y, i / n))
        # print(d_x)

        val = d_x + val
        x = val[0, 0]
        y = val[1, 0]

    return val


def Newton2(st):
    eps = 0.0001

    val = approx(st)
    
    x = val[0, 0]
    y = val[1, 0]

    while 1:
        d_x = LU(f_p(x, y), (-1) * f(x, y, 1))
        val_ = d_x + val

        if abs(val_ - val) < eps:
            return val_

        val = val_
        x = val[0, 0]
        y = val[1, 0]


# примерное расположение решения по графику
st = Matrix([[3.0], [-1.3]])
print("Решение")
k = Newton2(st)

print(k)

print("Подставили")
print(f(k[0, 0], k[1, 0], 1))



