from random import random

from settings import *
from math_model import *

from math import cos, sin, atan2, sqrt, pi, factorial
import numpy as np
from scipy.optimize import minimize, root_scalar, fsolve
import matplotlib.pyplot as plt
import matplotlib.image as mp_img
import csv


def get_points_by_solution(a, steps_th = 250):
    r, th = [], []

    for i in range(1, steps_th):
        theor_th = pi * i / steps_th

        def theor_err(r):
            return err(r, theor_th, a[:K], a[K:])

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

        r.append(theor_r)
        th.append(theor_th)
        # print(f"{theor_r}, {theor_th}")
    return r, th


def not_closer_than(r, th, min_r):
    result_r, result_th = [], []
    for i in range(len(r)):
        if r[i] >= min_r:
            result_r.append(r[i])
            result_th.append(th[i])
    return result_r, result_th


def show_artificial_isochrome(a, steps_th = 250, points = None, points_file = None,
                              show_graphics = True, show_photo = True):

    theor_r, theor_th = get_points_by_solution(a, steps_th)

    x1 = [theor_r[i] * cos(theor_th[i]) for i in range(len(theor_r))]
    y1 = [theor_r[i] * sin(theor_th[i]) for i in range(len(theor_r))]

    if show_graphics:
        plt.plot(x1, y1, label="line 1")
        plt.show()

    if show_photo:
        img = mp_img.imread(image_file)
        img = np.array(img)

        for i in range(len(x1)):
            x = round((x1[i]+center[0]) * img.shape[1] / w)
            y = round((-y1[i]+center[1]) * img.shape[1] / w)
            img[y, x, :] = 255

        if points_file is not None:
            points = []
            with open(points_file) as csv_file:
                reader = csv.reader(csv_file, quoting=csv.QUOTE_NONNUMERIC)
                for row in reader:
                    points.append(row)

        if points is not None:
            for i in range(len(points)):
                x = round(points[i][0] * img.shape[1] / w)
                y = round(points[i][1] * img.shape[1] / w)
                img[y, x, 0] = 255

        fig = plt.figure(figsize=(10, 10 / img.shape[1] * img.shape[0]), dpi=300)

        plt.imshow(img)
        plt.show()


def add_noise_polar(r, th, noise_max_r):
    for i in range(len(r)):
        noise_th = random()*2*pi
        noise_r = random()*noise_max_r
        x = r[i] * cos(th[i]) + noise_r * cos(noise_th)
        y = r[i] * sin(th[i]) + noise_r * sin(noise_th)
        r[i] = sqrt(x**2 + y**2)
        th[i] = atan2(y, x)
