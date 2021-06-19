from settings import *
from math_model import *

from math import cos, sin, atan2, sqrt, pi, factorial
import numpy as np
from scipy.optimize import minimize, root_scalar, fsolve
import matplotlib.pyplot as plt
import matplotlib.image as mp_img
import csv

# Начальное приближение
a0 = theoretical_solution(sigma_22_inf=sigma_22_inf, crack_width=crack_width, alpha=0)

steps_th = 50
theor_points = []

for i in range(1, steps_th):
    theor_th = pi * i / steps_th

    def theor_err(r):
        global theor_th, a0
        return err(r, theor_th, a0[:K], a0[K:])

    # left_border, right_border = 0.00001, 5
    left_border, right_border = 0.00001, 5

    max_steps = 10
    while theor_err(left_border) < 0 and max_steps > 0:
        left_border /= 10
        max_steps -= 1

    if max_steps == 0:
        continue

    max_steps = 10
    while theor_err(right_border) > 0 and max_steps > 0:
        left_border *= 10
        max_steps -= 1

    if max_steps == 0:
        continue

    # x1 = [left_border + (right_border-left_border)*i/1000 for i in range(1000)]
    # y1 = [theor_err(x1[i]) for i in range(1000)]
    # plt.plot(x1, y1, label="line 1")
    # plt.show()
    # exit()

    sol = root_scalar(theor_err, bracket=[left_border, right_border], method='brentq')
    theor_r = sol.root

    theor_points.append([theor_r, theor_th])
    # print(f"{theor_r}, {theor_th}")


x1 = [point[0]*cos(point[1]) for point in theor_points]
y1 = [point[0]*sin(point[1]) for point in theor_points]
plt.plot(x1, y1, label="line 1")
plt.show()

img = mp_img.imread(image_file)
img = np.array(img)

for i in range(len(x1)):
    x = round((x1[i]+center[0]) * img.shape[1] / w)
    y = round((-y1[i]+center[1]) * img.shape[1] / w)
    img[y, x, :] = 255

points = []
with open(points_file) as csv_file:
    reader = csv.reader(csv_file, quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        points.append(row)

for i in range(len(points)):
    x = round(points[i][0] * img.shape[1] / w)
    y = round(points[i][1] * img.shape[1] / w)
    img[y, x, 0] = 255


fig = plt.figure(figsize=(10, 10 / img.shape[1] * img.shape[0]), dpi=300)

plt.imshow(img)
plt.show()
