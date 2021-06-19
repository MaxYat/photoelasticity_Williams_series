from settings import *
from math_model import *

from math import cos, sin, atan2, sqrt, pi, factorial
import numpy as np
from scipy.optimize import minimize, root_scalar, fsolve
import matplotlib.pyplot as plt
import matplotlib.image as mp_img


steps_th = 100
theor_points = []

for i in range(1, steps_th):
    theor_th = pi * i / steps_th

    def theor_err(r):
        global theor_th, a0
        return err(r, theor_th, a0[:K], a0[K:])

    sol = root_scalar(theor_err, bracket=[0.0001, 5], method='brentq')

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

img = mp_img.imread(image_file)
plt.imshow(img)
