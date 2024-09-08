"""
    问题3代码
    取标准单位,长度为m
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

plt.rcParams["font.sans-serif"] = ["SimHei"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题

theta_all = np.zeros(224)
x_all = np.zeros(224)
y_all = np.zeros(224)

Collision = 0.15

R = 4.5  # 回转半径

r_count = 0  # 螺距二分次数

r_0 = 0.46  # 螺距上限初始化
r_1 = 0.43  # 螺距下限初始化

LEGNTH = 1000  # 角度遍历步进


def theta_solve_func(x, thetai, l, b):
    """
    计算龙身位置函数
    """
    return (
        ((b * thetai) * (b * thetai))
        + ((b * x) * (b * x))
        - 2 * (b * thetai) * (b * x) * np.cos(x - thetai)
        - l * l
    )


def find_range(arr, target):
    """
    检查 target 在 arr 内的哪两个值之间
    """
    for i in range(len(arr) - 1):
        if arr[i] <= target < arr[i + 1]:
            return (i, i + 1)

    return None


while True:
    # 二分更新
    r_temp = (r_0 + r_1) / 2
    b_0 = r_temp / (2 * np.pi)

    print("------%d次二分计算,螺距%.2f-------" % (r_count, r_temp))

    # 角度遍历初始化
    theta_start = R / b_0
    theta_fail = 0
    # 角度遍历
    for j in range(0, LEGNTH):
        # theta_all[0] = R / b_0 + 0.59 * np.pi
        theta_all[0] = theta_start + (j * (np.pi) * 2 / (LEGNTH))
        x_all[0] = b_0 * theta_all[0] * np.cos(theta_all[0])
        y_all[0] = b_0 * theta_all[0] * np.sin(theta_all[0])
        print("------%d次角度计算,螺距%f,初始角度%f-------" % (j, r_temp, theta_all[0]))
        for i in range(1, 224):
            if i == 1:
                result = opt.root(
                    lambda x: theta_solve_func(x, theta_all[0], 2.86, b_0),
                    (theta_all[0] + 1),
                )
                x_all[i] = b_0 * (result.x[0]) * np.cos(result.x[0])
                y_all[i] = b_0 * (result.x[0]) * np.sin(result.x[0])
                theta_all[i] = result.x[0]
            else:
                result = opt.root(
                    lambda x: theta_solve_func(x, theta_all[i - 1], 1.65, b_0),
                    (theta_all[i - 1] + 1),
                )
                x_all[i] = b_0 * (result.x[0]) * np.cos(result.x[0])
                y_all[i] = b_0 * (result.x[0]) * np.sin(result.x[0])
                theta_all[i] = result.x[0]

        theta_d = theta_all[0] + 2 * (np.pi)

        # 寻找索引
        (index_a1, index_b1) = find_range(theta_all, theta_d)
        x_a1 = b_0 * theta_all[index_a1] * np.cos(theta_all[index_a1])
        x_b1 = b_0 * theta_all[index_b1] * np.cos(theta_all[index_b1])
        y_a1 = b_0 * theta_all[index_a1] * np.sin(theta_all[index_a1])
        y_b1 = b_0 * theta_all[index_b1] * np.sin(theta_all[index_b1])

        # 计算此两点直角坐标系直线
        _k1 = (y_b1 - y_a1) / (x_b1 - x_a1)
        cos_beta1 = (np.abs(x_b1 - x_a1)) / (
            np.sqrt((y_b1 - y_a1) * (y_b1 - y_a1) + (x_b1 - x_a1) * (x_b1 - x_a1))
        )
        _b1 = -_k1 * x_a1 + y_a1

        # 计算外点坐标
        phi1 = np.arcsin(
            (b_0 * theta_all[1] * np.sin(theta_all[1] - theta_all[0])) / 2.86
        )
        psi1 = 2 * np.pi - (np.pi) / 2 - np.arctan(27.5 / 15) - phi1
        r_c1 = np.sqrt(
            (b_0 * theta_all[0]) * (b_0 * theta_all[0])
            + 0.098125
            - 2 * (b_0 * theta_all[0]) * np.sqrt(0.098125) * np.cos(psi1)
        )
        theta_c1 = theta_all[0] - np.arcsin(np.sqrt(0.098125) * np.sin(psi1) / r_c1)
        y_c1 = r_c1 * np.sin(theta_c1)
        x_c1 = r_c1 * np.cos(theta_c1)

        d_1 = np.abs(_k1 * x_c1 - y_c1 + _b1) / np.sqrt(_k1 * _k1 + 1)

        theta_d = theta_all[0] + 2 * (np.pi)

        # 寻找索引
        (index_a2, index_b2) = find_range(theta_all, theta_d)
        x_a2 = b_0 * theta_all[index_a2] * np.cos(theta_all[index_a2])
        x_b2 = b_0 * theta_all[index_b2] * np.cos(theta_all[index_b2])
        y_a2 = b_0 * theta_all[index_a2] * np.sin(theta_all[index_a2])
        y_b2 = b_0 * theta_all[index_b2] * np.sin(theta_all[index_b2])

        # 计算此两点直角坐标系直线
        _k2 = (y_b2 - y_a2) / (x_b2 - x_a2)
        cos_beta2 = (np.abs(x_b2 - x_a2)) / (
            np.sqrt((y_b2 - y_a2) * (y_b2 - y_a2) + (x_b2 - x_a2) * (x_b2 - x_a2))
        )
        _b2 = -_k2 * x_a2 + y_a2

        # 计算外点坐标
        phi2 = np.arcsin(
            (b_0 * theta_all[0] * np.sin(theta_all[1] - theta_all[0])) / 2.86
        )
        psi2 = 2 * np.pi - (np.pi) / 2 - np.arctan(27.5 / 15) - phi2
        r_c2 = np.sqrt(
            (b_0 * theta_all[1]) * (b_0 * theta_all[1])
            + 0.098125
            - 2 * (b_0 * theta_all[1]) * np.sqrt(0.098125) * np.cos(psi2)
        )
        theta_c2 = theta_all[1] + np.arcsin(np.sqrt(0.098125) * np.sin(psi2) / r_c2)
        y_c2 = r_c2 * np.sin(theta_c2)
        x_c2 = r_c2 * np.cos(theta_c2)

        d_2 = np.abs(_k2 * x_c2 - y_c2 + _b2) / np.sqrt(_k2 * _k2 + 1)

        if d_1 < Collision or d_2 < Collision:
            theta_fail = 1
            break

        # __r = np.arange(0, 12, 0.0001)
        # __theta = __r / b_0
        # __ax = plt.subplot(111, projection="polar")
        # __ax.plot(__theta, __r, linewidth=1, color="red")
        # __ax.scatter(theta_c1, r_c1, c="g", marker="o", s=4)
        # __ax.scatter(theta_c1, r_test_c1, c="purple", marker="o", s=4)
        # __ax.scatter(theta_all, b_0 * theta_all, c="b", marker="o", s=4)
        # __ax.plot(theta_all, b_0 * theta_all, c="b", linewidth=0.8)
        # __ax.grid(True)
        # plt.show()
        # break

    if theta_fail == 1:
        r_1 = r_temp
        theta_fail = 0
    else:
        r_0 = r_temp
        theta_fail = 0

    if r_count == 10:
        break

    r_count = r_count + 1

print(r_0)

b = r_0 / (2 * np.pi)

plt.figure(figsize=(5.5, 5.5))
__theta = np.arange(0, (21 * np.pi * 2), 0.0001)
__x = b * (__theta) * np.cos(__theta)
__y = b * (__theta) * np.sin(__theta)
plt.plot(__x, __y, linewidth=1, color="red", label="螺线")
__theta_all = np.zeros(224)

for i in range(0, 224):
    __theta_all[i] = theta_all[i]

__xi = b * (__theta_all) * np.cos(__theta_all)
__yi = b * (__theta_all) * np.sin(__theta_all)

plt.plot(__xi, __yi, linewidth=1, color="blue", label="龙身")
plt.scatter(__xi, __yi, linewidth=1, color="blue", marker="o", s=4)

__theta_cir = np.arange(0, (np.pi * 2), 0.0001)
__x_cir = 4.5 * np.cos(__theta_cir)
__y_cir = 4.5 * np.sin(__theta_cir)

plt.plot(__x_cir, __y_cir, linewidth=1, color="purple", label="掉头空间")

plt.legend()
plt.show()
