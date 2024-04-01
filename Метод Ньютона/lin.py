import math

def f(x):
    return x ** 2 - math.asin(x - 0.2)


def f_p(x):
    return 2 * x - 1 / (math.sqrt(1 - (x - 0.2) ** 2))


def posled_perebor():
    eps = 0.0001
    a = -0.8
    b = 10
    
    # Последовательный перебор
    N = 10
    h = 1
    k = 0
    otrezok = []
    while (k < N):
        h = (b - a) / N
        x_k = a + k * h
        x_k1 = a + (k + 1) * h
        k += 1
        if abs(x_k - x_k1) < 0.1:
            otrezok.append(a)
            otrezok.append(b)
            return otrezok
        if f(x_k) * f(x_k1) < 0:
            a = x_k
            b = x_k1
            k = 0


def Newton1():
    eps = 0.0001
    otr = posled_perebor()
    print("отрезок")
    print(posled_perebor())
    a = otr[0]
    b = otr[1]
    x = a
    while 1:
        x_ = x - f(x) / f_p(x)

        if not (a <= x_ <= b):  # если x_k+1 вылетел за пределы промежутка, то возвращаем его обратно
            x_ = (a + b) / 2

        if f(a) * f(x_) < 0:
            b = x_
        else:
            a = x_

        if abs(x_ - x) < eps:
            return x_
        x = x_

print("Границы, найденные с помощью последовательного перебора:")
print(posled_perebor())
print("Решение методом Ньютона")
print(Newton1())
#Проверка
print("Подставим решение в уравнение")
print(f(Newton1()))
