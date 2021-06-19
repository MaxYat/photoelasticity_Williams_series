from settings import *
from math import cos, sin, atan2, sqrt, pi, factorial


def theoretical_solution(sigma_22_inf, crack_width, alpha):
    a = np.zeros((K+M))
    for n in range(0, (K - 1) // 2):
        a[2*n] = (-1)**(n+1) * factorial(2*n) * sigma_22_inf \
            / (2**(3*n+1/2) * (factorial(n)**2) *(2*n-1) * (crack_width / 2) ** (n-1/2))

    a[1] = sigma_22_inf * (alpha - 1) / 4

    return a


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


def derr_da1(r, th, p, a1, a2):
    return 2 * r ** (p / 2 - 1) * (
            sigma11_minus_sigma22(r, th, a1, a2) * (f_1_11(p, th) - f_1_22(p, th)) +
            4 * sigma12(r, th, a1, a2) * f_1_12(p, th)
    )


def derr_da2(r, th, p, a1, a2):
    return 2 * r ** (p / 2 - 1) * (
            sigma11_minus_sigma22(r, th, a1, a2) * (f_2_11(p, th) - f_2_22(p, th)) +
            4 * sigma12(r, th, a1, a2) * f_2_12(p, th)
    )
