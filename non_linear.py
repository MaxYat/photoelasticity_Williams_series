from math import cos, sin, atan2, sqrt
import csv
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import lstsq

# === Settings === #

example = 0

N = 3
h = [5.7, 5.0][example]  # mm
fs = 1.823 # kg/mm из статьи Л. В. Степановой и В. С. Долгих
# fs = 1.945659 # из программы тарировки


# K, M = 5, 3
K, M = 7, 5
# K, M = 15, 9

points_file = [
    f"data/One_crack/points_Р01-100кг_N_{N}_of_5.txt",
    f"data/Two_cracks/Points_Р 11-100кг_N_{N}_of_4.txt"
    ][example]

center = [
    [34.29787, 51.604862],
    [38.214905, 31.455807]
    ][example]

# center = [34.3769, 51.604862]

# Начальное приближение
a0 = [0.1 for k in range(K+M)]

# MINIMIZE_METHODS = ['nelder-mead', 'powell', 'cg', 'bfgs', 'newton-cg',
#                     'l-bfgs-b', 'tnc', 'cobyla', 'slsqp', 'trust-constr',
#                     'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov']

# ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'TNC', 'trust-constr']

minimize_method = 'nelder-mead'
minimize_tolerance = 1e-8
minimize_max_iteration = 5000

overdetermined_method_starts_from_previous_result = True
overdetermined_method_max_iterations = 100
overdetermined_method_precision = 1E-15

# === End of settings === #

points = []
with open(points_file) as csv_file:
    reader = csv.reader(csv_file, quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        points.append(row)


r = [sqrt((point[0]-center[0])**2 + (point[1]-center[1])**2) for point in points]
th = [(atan2(-point[1]+center[1], point[0]-center[0])) for point in points]


def f_1_11(k, th):
    return 0.5 * k * (
            (2+k/2+(-1)**k)*cos((k/2-1)*th)
            - (k/2 - 1)*cos((k/2-3)*th)
    )

def f_1_22(k, th):
    return 0.5 * k * (
            (2-k/2-(-1)**k)*cos((k/2-1)*th)
            + (k/2 - 1)*cos((k/2-3)*th)
    )

def f_1_12(k, th):
    return 0.5 * k * (
            -(k/2+(-1)**k)*sin((k/2-1)*th)
            + (k/2 - 1)*sin((k/2-3)*th)
    )

def f_2_11(k, th):
    return -0.5 * k * (
            (2+k/2-(-1)**k)*sin((k/2-1)*th)
            - (k/2 - 1)*sin((k/2-3)*th)
    )

def f_2_22(k, th):
    return -0.5 * k * (
            (2-k/2+(-1)**k)*sin((k/2-1)*th)
            + (k/2 - 1)*sin((k/2-3)*th)
    )

def f_2_12(k, th):
    return 0.5 * k * (
            -(k/2-(-1)**k)*cos((k/2-1)*th)
            + (k/2 - 1)*cos((k/2-3)*th)
    )

def sigma11_minus_sigma22(r, th, a1, a2):
    return sum([a1[k-1]*r**(k/2-1)*(f_1_11(k, th) - f_1_22(k, th)) for k in range(1,K+1)]) \
        + sum([a2[k - 1] * r ** (k / 2 - 1) * (f_2_11(k, th) - f_2_22(k, th)) for k in range(1, M + 1)])

def sigma12(r, th, a1, a2):
    return sum([a1[k-1]*r**(k/2-1)*f_1_12(k, th) for k in range(1,K+1)]) \
               +sum([a2[k-1]*r**(k/2-1)*f_2_12(k, th) for k in range(1,M+1)])

left = (N*fs/h) ** 2

def err(r, th, a1, a2):
    return sigma11_minus_sigma22(r, th, a1, a2) ** 2 \
        + 4 * sigma12(r, th, a1, a2) ** 2 - left


def objective(a):
    return sum([err(r[i], th[i], a[:K], a[K:])**2 for i in range(len(r))])


res = minimize(objective, a0, method=minimize_method, tol=minimize_tolerance, options={"maxiter" : minimize_max_iteration})

print(res.x)
print(objective(res.x))

def derr_da1(r, th, p, a1, a2):
     return 2 * r**(p/2-1) * (
        sigma11_minus_sigma22(r, th, a1, a2) * (f_1_11(p, th)-f_1_22(p, th)) +
        4 * sigma12(r, th, a1, a2) * f_1_12(p, th)
     )


def derr_da2(r, th, p, a1, a2):
    return 2 * r**(p/2-1) * (
        sigma11_minus_sigma22(r, th, a1, a2) * (f_2_11(p, th)-f_2_22(p, th)) +
        4 * sigma12(r, th, a1, a2) * f_2_12(p, th)
    )


a = res.x if overdetermined_method_starts_from_previous_result else a0[:]

for iteration in range(overdetermined_method_max_iterations):
    B = [-err(r[i], th[i], a[:K], a[K:]) for i in range(len(r))]
    A = [[derr_da1(r[i], th[i], k, a[:K], a[K:]) if k<=K else derr_da2(r[i], th[i], k-K, a[:K], a[K:])
          for k in range(1, K+M+1)] for i in range(len(r))]

    a_prev = a
    delta_a, _, _, _ = lstsq(A, B)
    a = a+delta_a

    if sum((a-a_prev)**2) < overdetermined_method_precision:
        break

print(a)
print(objective(a))

