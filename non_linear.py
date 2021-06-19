from settings import *
from math_model import *

from math import cos, sin, atan2, sqrt, pi, factorial
import csv
import numpy as np
from scipy.optimize import minimize, root_scalar, fsolve
from scipy.linalg import lstsq
import matplotlib.pyplot as plt

# === Numeric methods settings === #

# MINIMIZE_METHODS = ['nelder-mead', 'powell', 'cg', 'bfgs', 'newton-cg',
#                     'l-bfgs-b', 'tnc', 'cobyla', 'slsqp', 'trust-constr',
#                     'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov']

# ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'TNC', 'trust-constr']

minimize_method = 'l-bfgs-b'
# 'cg' 4, 'l-bfgs-b' 4+, 'tnc' 4, 'cobyla' 3, 'slsqp' strange 4+

minimize_tolerance = 1e-8
minimize_max_iteration = 5000

overdetermined_method_starts_from_previous_result = False
overdetermined_method_max_iterations = 100
overdetermined_method_precision = 1E-15

# === End of Numeric methods settings === #

points = []
with open(points_file) as csv_file:
    reader = csv.reader(csv_file, quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        points.append(row)


r = [sqrt((point[0]-center[0])**2 + (point[1]-center[1])**2) for point in points]
# th = [abs(atan2(-point[1]+center[1], point[0]-center[0])) for point in points]
th = [(atan2(-point[1]+center[1], point[0]-center[0])) for point in points]


def squared_error(a):
    return sum([err(r[i], th[i], a[:K], a[K:])**2 for i in range(len(r))])


def print_results(name, a):
    print(name)
    print("a1 = ", a[:K])
    print("a2 = ", a[K:])
    print("MSE = ", squared_error(a) / len(r), "\n")


# Начальное приближение
a0 = theoretical_solution(sigma_22_inf=sigma_22_inf, crack_width=crack_width, alpha=0)

print_results("Теоретическое решение для бесконечной пластины", a0)

res = minimize(squared_error, a0, method=minimize_method, tol=minimize_tolerance, options={"maxiter" : minimize_max_iteration})

print_results("Минимизация квадратичного отклонения", res.x)

a = res.x if overdetermined_method_starts_from_previous_result else a0[:]

for iteration in range(overdetermined_method_max_iterations):
    B = [-err(r[i], th[i], a[:K], a[K:]) for i in range(len(r))]
    A = [[derr_da1(r[i], th[i], k, a[:K], a[K:]) if k<=K else derr_da2(r[i], th[i], k-K, a[:K], a[K:])
          for k in range(1, K+M+1)] for i in range(len(r))]

    a_prev = a
    delta_a, _, _, _ = lstsq(A, B)
    a = a+delta_a

    if sum((a-a_prev)**2)/len(a) < overdetermined_method_precision:
        break

print_results("Переопределённый метод", a)


steps_th = 100
theor_points = []

for i in range(1, steps_th):
    theor_th = pi * i / steps_th

    def theor_err(r):
        global theor_th, a0
        return err(r, theor_th, a0[:K], a0[K:])

    sol = root_scalar(theor_err, bracket=[0.0001, 5], method='brentq')
    # sol = root_scalar(theor_err, x0 = 0.1, fprime = True, fprime2 = True, method = 'halley')
    # sol = root_scalar(theor_err, x0=0.1, method = 'newton')
    # theor_r = sol.root[0]

    # sol = minimize(theor_err()**2)
    # sol = fsolve(theor_err(theor_th)**2, np.array([0.1]))

    # x1 = [0.01 + 7*i/1000 for i in range(1000)]
    # y1 = [theor_err(x1[i]) for i in range(1000)]
    # plt.plot(x1, y1, label="line 1")
    # plt.show()
    #
    # break

    theor_r = sol.root

    theor_points.append([theor_r, theor_th])
    print(f"{theor_r}, {theor_th}")




x1 = [point[0]*cos(point[1]) for point in theor_points]
y1 = [point[0]*sin(point[1]) for point in theor_points]
plt.plot(x1, y1, label="line 1")
plt.show()

