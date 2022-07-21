import sympy as sym
import math
import numpy as np
import matplotlib.pyplot as plt

# Область определения (и непрерывности) f(x)
ab = [-1, 1]
x = sym.Symbol('x')


# Задание функций, которые будем интерполировать
def f(xi):
    return xi**2 + 1 - math.acos(xi)


def h(xi):
    return f(xi) * abs(xi)


# Функция для построения полинома Лагранжа на заданной таблице xi
def lagrange(xi, g):
    l = 0
    for i in range(len(xi)):
        li = 0
        for j in range(len(xi)):
            if j == 0 or j == 1 and i == 0:
                li = (x - xi[j]) * g(xi[i])
            elif j != i:
                li *= (x - xi[j]) / (xi[i] - xi[j])
                li = sym.expand_multinomial(li)
        l += li
        l = sym.simplify(l)

    return l


# Функция для иного способа задания узлов. Нельзя пользоваться этой функцией  при [a, b] = [-1, 1], поэтому этот способ
# задания узлов не будет рассматриваться для заданной функции (полученные узлы будут вне области определения [-1, 1])
def another_nodes(n):
    xj = []
    for i in range(n):
        xj.append(1/2 * ((ab[1] - ab[0])*sym.cos((2*i + 1)/(2*(n + 1))*sym.pi)+(ab[1]-ab[0])))
    return xj


# Пусть "некоторое число полиномов" будет 4
x1i = np.linspace(ab[0], ab[1], 10)
y1f = lagrange(x1i, f)
y1h = lagrange(x1i, h)
print("f(x), n=10")
print(y1f)
print("h(x), n=10")
print(y1h)
x2i = np.linspace(ab[0], ab[1], 6)
print("f(x), n=6")
y2f = lagrange(x2i, f)
print("h(x), n=6")
y2h = lagrange(x2i, h)
print(y2f)
print(y2h)
x3i = np.linspace(ab[0], ab[1], 3)
print("f(x), n=3")
y3f = lagrange(x3i, f)
print("h(x), n=3")
y3h = lagrange(x3i, h)
print(y3f)
print(y3h)
x4i = np.linspace(ab[0], ab[1], 35)
print("f(x), n=35")
y4f = lagrange(x4i, f)
print("h(x), n=35")
y4h = lagrange(x4i, h)
print(y4f)
print(y4h)
x5i = np.linspace(ab[0], ab[1], 34)
print("f(x), n=34")
y5f = lagrange(x5i, f)
print("h(x), n=34")
y5h = lagrange(x5i, h)
print(y5f)
print(y5h)


lam_y1f = sym.lambdify(x, y1f, modules=['numpy'])
lam_y2f = sym.lambdify(x, y2f, modules=['numpy'])
lam_y3f = sym.lambdify(x, y3f, modules=['numpy'])
lam_y1h = sym.lambdify(x, y1h, modules=['numpy'])
lam_y4f = sym.lambdify(x, y4f, modules=['numpy'])
lam_y5f = sym.lambdify(x, y5f, modules=['numpy'])

x_vals = np.linspace(-1, 1, 100)
y1f_vals = lam_y1f(x_vals)
y2f_vals = lam_y2f(x_vals)
y3f_vals = lam_y3f(x_vals)
y4f_vals = lam_y4f(x_vals)
y5f_vals = lam_y5f(x_vals)
y1h_vals = lam_y1h(x_vals)

fig, ax = plt.subplots()
plt.grid(color='#666666', linestyle="--", linewidth=1)
ax.plot(x_vals, y1f_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'darkmagenta',
        label = 'f(x): n=10')

f1i = list(map(f, x1i))
ax.scatter(x1i, f1i,
        linewidth = 2,
        color = 'darkmagenta')

ax.plot(x_vals, y2f_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'darkblue',
        label = 'f(x): n=5')

f2i = list(map(f, x2i))
ax.scatter(x2i, f2i,
        linewidth = 2,
        color = 'darkblue')

ax.plot(x_vals, y3f_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'indigo',
        label = 'f(x): n=4')

f3i = list(map(f, x3i))
ax.scatter(x3i, f3i,
        linewidth = 2,
        color = 'indigo')

ax.plot(x_vals, y4f_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'red',
        label = 'f(x): n=35')

f4i = list(map(f, x4i))
ax.scatter(x4i, f4i,
        linewidth = 2,
        color = 'red')

ax.plot(x_vals, y5f_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'orange',
        label = 'f(x): n=34')

f5i = list(map(f, x5i))
ax.scatter(x5i, f5i,
        linewidth = 2,
        color = 'orange')


ax.legend()

fig.set_figwidth(12)
fig.set_figheight(6)
plt.show()

fig, ax = plt.subplots()

ax.plot(x_vals, y1f_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'darkmagenta',
        label = 'f(x): n=10')

fx = list(map(f, x1i))
ax.scatter(x1i, fx,
        linewidth = 2,
        color = 'darkmagenta')

ax.plot(x_vals, y1h_vals,
        linestyle = '-',
        linewidth = 2,
        color = 'magenta',
        label = 'h(x): n=10')

hx = list(map(h, x1i))
ax.scatter(x1i, hx,
        linewidth = 2,
        color = 'magenta')

ax.legend()

fig.set_figwidth(12)
fig.set_figheight(6)
plt.show()
