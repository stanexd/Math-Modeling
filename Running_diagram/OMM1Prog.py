import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from numpy import log, arctan, arccos, cos, pi, zeros, linspace

# Параметры сетки
T = 1 
Nx, Nt = 501, 501 
x = linspace(0, -1, Nx) 
t = linspace(0, T, Nt) 
h_x, tau = x[0] - x[1], t[1] - t[0] 
eps = 0.001

u_4 = zeros((Nx, Nt))

u_4[:, 0] = cos(pi * x / 2)
u_4[0, :] = 1 + 0.5 * arctan(t)

#Блок для 4-х точечной схемы

def f(i, j):
    return log(1 + u_4[i, j] ** 2)

def F_proizv(i, j):
    return 1 / (2 * tau) - u_4[i, j] / ((1 + u_4[i, j] ** 2) * h_x) 

def F(i, j):
    return (u_4[i - 1, j] - u_4[i - 1, j - 1] + u_4[i, j] - u_4[i - 1, j]) / (2 * tau) - (f(i, j - 1) - f(i - 1, j - 1) + f(i, j) - f(i, j - 1)) / (2 * h_x)

def Newton(i, j):
    y_0 = u_4[i-1, j-1]
    y = y_0 - F(i, j) / F_proizv(i, j)
    while (abs(y_0 - y) > eps):
        y_0 = y
        u_4[i, j] = y_0
        y = y_0 - F(i, j) / F_proizv(i, j)
    u_4[i, j] = y

for j in range(1, Nt):
    for i in range(1, Nx):
        Newton(i, j)

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# t, x = np.meshgrid(t, x)
# surf = ax.plot_surface(x, t, u_4, cmap='viridis')
# plt.xlabel('x')
# plt.ylabel('t')
# ax.set_zlabel('u')
# ax.set_title('График численного решения для 4-х точечной схемы')
# plt.show()


# plt.plot(t, u_4[0, :])
# plt.xlabel('t')
# plt.ylabel('u')
# plt.grid('both')
# plt.title('График u(x,t) при x = 0')
# plt.show()

# t_i = int((Nt - 1) * 0.2)
# x_i = int((Nx - 1) * 0.7)

# print(f'Шаг: {h_x}, x: {x[x_i]}, t: {t[t_i]}, u_4(-{x_i/(Nx - 1)}, {t_i/(Nt - 1)}) = {u_4[x_i, t_i]}')


# Шаг: 0.01, x: -0.7, t: 0.2, u_4(-0.7, 0.2) = 0.707058442260344 u_an = 0.692800686414118
# Шаг: 0.005, x: -0.7, t: 0.2, u_4(-0.7, 0.2) = 0.7070163049033238 u_an = 0.692800686414118
# Шаг: 0.0033, x: -0.7, t: 0.2, u_4(-0.7, 0.2) = 0.7069791468141973 u_an = 0.692800686414118
# Шаг: 0.0025, x: -0.7, t: 0.2, u_4(-0.7, 0.2) = 0.7069460469234494 u_an = 0.692800686414118
# Шаг: 0.002, x: -0.7, t: 0.2, u_4(-0.7, 0.2) = 0.706916300542111 u_an = 0.692800686414118


steps = [0.01, 0.005, 0.0033, 0.0025, 0.002]
u_an = 0.692800686414118
u_4_values = [0.707058442260344, 0.7070163049033238, 0.7069791468141973, 0.7069460469234494, 0.706916300542111]

diffs = [abs(u_an - u_4) for u_4 in u_4_values]

plt.figure(figsize=(10, 6))
plt.plot(steps, diffs, marker='o')
plt.xlabel('h_x')
plt.ylabel('σ')
plt.title('Зависимость ошибки σ от величины шага h_x')
plt.grid(True)
plt.show()
