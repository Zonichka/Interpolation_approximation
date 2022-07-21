import sympy as sym
import numpy as np
import matplotlib.pyplot as plt


# Область определения (и непрерывности) f(x)
a, b = -1, 1
x = sym.Symbol('x')


# Задание функций, которые будем интерполировать
def f(xj):
    return xj**2 * sym.cos(xj)


# Узлы и значения в узлах
xi = np.linspace(a, b, 5)
yi = list(map(f, xi))


# Метод наименьших квадратов на заданной таблице xi
def MNS(xi, yi):
    A = []
    for i in range(4):
        a = []
        for j in range(5):
            aij = 0
            for k in range(len(xi)):
                aij += xi[k]**(i+j)
            a.append(aij)
        bij = 0
        for k in range(len(xi)):
            bij += xi[k]**i * yi[k]
        a.append(bij)
        A.append(a)
    A = sym.Matrix(A)

    c0, c1, c2, c3 = sym.symbols('c0 c1 c2 c3')
    system = H, b = A[:, :-1], A[:, -1]
    ans = sym.linsolve(system, c0, c1, c2, c3)
    (a0, a1, a2, a3) = next(iter(ans))
    a = [a0, a1, a2, a3]

    l = a[3]*x**3 + a[2]*x**2 + a[1]*x + + a[0]

    return l


# Функция для подсчета многочленов Лежандра
def p(n):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return sym.simplify((2*n - 1)/n * x * p(n-1) - (n-1)/n * p(n-2))


# Коэффициенты Фурье функции f(x)
def c(n):
    I1 = sym.simplify((2*n+1) / 2 * sym.integrate(f(x) * p(n), (x, a, b)))

    return I1


# Полином Лежандра
def pL(n):
    pL = 0
    for i in range(n):
        pL += c(i) * p(i)
        pL = sym.simplify(pL)
    pL = sym.expand(pL)
    return pL


y1 = MNS(xi, yi)
print(y1)
y2 = pL(3)
print(y2)

lam_y1 = sym.lambdify(x, y1, modules=['numpy'])
lam_y2 = sym.lambdify(x, y2, modules=['numpy'])

x_vals = np.linspace(-1, 1, 100)
y1_vals = lam_y1(x_vals)
y2_vals = lam_y2(x_vals)

fig, ax = plt.subplots()
plt.grid(color='#666666', linestyle="--", linewidth=1)
ax.plot(x_vals, y1_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'darkmagenta',
        label = 'P by MSI')

ax.scatter(xi, yi,
        linewidth = 2,
        color = 'darkmagenta')

ax.plot(x_vals, y2_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'darkblue',
        label = 'P by Legandre polinom')

fi = list(map(f, x_vals))
ax.plot(x_vals, fi,
        linestyle='-',
        linewidth = 2,
        color = 'blue',
        label = 'f(x)')

ax.legend()

fig.set_figwidth(12)
fig.set_figheight(6)
plt.show()


f1_vals = []
for i in range(len(y1_vals)):
    f1_vals.append(fi[i]-y1_vals[i])

f2_vals = []
for i in range(len(y2_vals)):
    f2_vals.append(fi[i]-y2_vals[i])

fig, ax = plt.subplots()
plt.grid(color='#666666', linestyle="--", linewidth=1)

ax.plot(x_vals, f1_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'darkblue',
        label = 'f - f1')

ax.plot(x_vals, f2_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'magenta',
        label = 'f - f2')

ax.legend()

fig.set_figwidth(12)
fig.set_figheight(6)
plt.show()
