import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
from matplotlib import animation
from numpy import exp, sin, cos, pi, outer, zeros, linspace

# Параметры сетки
T = 1 
Nx, Ny, Nt = 30, 30, 100 
x = linspace(0, pi, Nx) 
y = linspace(0, pi, Ny) 
t = linspace(0, T, Nt)  
h_x, h_y, tau = x[1] - x[0], y[1] - y[0], t[1] - t[0] 
gamma_x, gamma_y = tau / (h_x ** 2), tau / (h_y ** 2) 

# Инициализация массива для решения
u = zeros((Nx, Ny, 2 * Nt + 1)) 

# Инициализация массива для аналитического решения
u_an = zeros((Nx, Ny, 2 * Nt + 1))

# Начальные условия
u[:, :, 0] = outer(sin(2 * x), cos(y))

#аналитическое решение
for i1 in range(Nx):
    for i2 in range(Ny):
        for j in range(0, 2 * Nt, 2):
            u_an[i1, i2, j] = ((exp(t[j // 2]) - exp(-2 * t[j // 2])) / 3 ) * sin(x[i1]) * cos(y[i2]) + sin(2 * x[i1]) * cos(y[i2]) * exp(-5 * t[j // 2])

# Функции для расчета значений
def F_1(i1, i2, j):
    return 0.5 * gamma_y * (u[i1, i2 - 1, j - 1] + u[i1, i2 + 1, j - 1]) + (1 - gamma_y) * u[i1, i2, j - 1] + 0.5 * tau * exp(tau * (j + 1) / 2) * sin(x[i1]) * cos(y[i2])

def F_2(i1, i2, j):
    return 0.5 * gamma_x * (u[i1 - 1, i2, j - 1] + u[i1 + 1, i2, j - 1]) + (1 - gamma_x) * u[i1, i2, j - 1] + 0.5 * tau * exp(tau * (j - 1) / 2) * sin(x[i1]) * cos(y[i2])

#прогоночные функции 
def progonka_x(i2, j):
    d = zeros(Nx)
    sigma = zeros(Nx)
    d[1], sigma[1] = 0, 0
    A, B, C = 0.5 * gamma_y, 1 + gamma_y, 0.5 * gamma_y
    for m in range(1, Nx - 1):
        Fm = - F_1(m, i2, j)
        d[m + 1] = C / (B - A * d[m])
        sigma[m + 1] = (Fm - A * sigma[m]) / (A * d[m] - B)
    u[Nx - 1, i2, j] = 0

    for m in range(Nx - 1, 0, -1):
        u[m - 1, i2, j] = d[m] * u[m, i2, j] + sigma[m]

def progonka_y(i1, j):
    d = zeros(Ny)
    sigma = zeros(Ny)
    d[1], sigma[1] = 1, 0
    A, B, C = 0.5 * gamma_y, 1 + gamma_y, 0.5 * gamma_y
    for m in range(1, Ny - 1):
        Fm = - F_2(i1, m, j)
        d[m + 1] = C / (B - A * d[m])
        sigma[m + 1] = (Fm - A * sigma[m]) / (A * d[m] - B)
    u[i1, Ny - 1, j] = sigma[-1] / (1 - d[-1])

    for m in range(Ny - 1, 0, -1):
        u[i1, m - 1, j] = d[m] * u[i1, m, j] + sigma[m]

# Расчет значений сетки
for j in range(1, 2 * Nt, 2):
    for i2 in range(1, Ny - 1):
        progonka_x(i2, j)
    for i1 in range(1, Nx - 1):
        progonka_y(i1, j + 1)

# # Функция создание 3D графика
# def graph(u_test, z_min, z_max):
#     X, Y = np.meshgrid(x, y)  # Создаем сетку
#     fig = plt.figure()  # Создаем фигуру для графика
#     ax = fig.add_subplot(111, projection='3d')  # Добавляем оси для 3D графика

#     # Функция для обновления графика
#     def update_graph(num):
#         ax.clear()  # Очищаем оси перед каждым обновлением
#         ax.set_zlim([z_min, z_max]) #устанавливаем ограничение оси z
#         Z = u_test[:, :,2 * num]  # Выбираем момент времени для визуализации
#         ax.plot_surface(X, Y, Z, cmap='viridis')  # Строим поверхность
#         ax.set_xlabel('X')
#         ax.set_ylabel('Y')
#         ax.set_zlabel('U')
#         ax.set_title(f'График в момент времени {tau * num:.2f}')

#     # Создание анимации
#     ani = animation.FuncAnimation(fig, update_graph, frames=range(Nt), interval=1)

#     # Показываем график
#     return plt.show()
        
# graph(u - u_an, -0.02, 0.02)

mom = 99
T_for_graph, X = np.meshgrid(t, x)  # Создаем сетку
Z = u[:, 20, :2 * Nt:2]
fig = plt.figure()  # Создаем фигуру для графика
ax = fig.add_subplot(111, projection='3d')  # Добавляем оси для 3D графика
ax.set_xlabel('T')
ax.set_ylabel('X')
ax.set_zlabel('U')
ax.set_title(f'График численного решения при y = 2 * π / 3')
ax.plot_surface(T_for_graph, X, Z, cmap = 'viridis')  # Строим 3D-график
plt.show()

