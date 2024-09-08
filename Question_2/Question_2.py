"""
    问题2代码
    取标准单位,长度为m
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

plt.rcParams["font.sans-serif"] = ["SimHei"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题

# 等距螺线方程为：r = bθ
b = 0.55 / (2 * np.pi)  # 螺距参数


theta_0 = 8.5 * (np.pi)  # 龙头角度初始估计值
## 遍历步长
D = 3000
delta_theta = 0.5 * (np.pi) / D

theta_all = np.zeros(224)
x_all = np.zeros(224)
y_all = np.zeros(224)

theta_all[0] = theta_0
x_all[0] = b * theta_0 * np.cos(theta_0)
y_all[0] = b * theta_0 * np.sin(theta_0)

flag = 0  # 遍历终止标志位
count = 0  # 遍历次数标志位


def theta_solve_func(x, thetai, l):
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


while flag == 0:
    print("%d次遍历计算" % (count))

    print("---------1. 更新龙身位置---------")

    if count != 0:
        theta_all[0] = theta_0
        x_all[0] = b * theta_0 * np.cos(theta_0)
        y_all[0] = b * theta_0 * np.sin(theta_0)

    for i in range(1, 224):
        if i == 1:
            result = opt.root(
                lambda x: theta_solve_func(x, theta_0, 2.86), (theta_0 + 1)
            )
            x_all[i] = b * (result.x[0]) * np.cos(result.x[0])
            y_all[i] = b * (result.x[0]) * np.sin(result.x[0])
            theta_all[i] = result.x[0]
        else:
            result = opt.root(
                lambda x: theta_solve_func(x, theta_all[i - 1], 1.65),
                (theta_all[i - 1] + 1),
            )
            x_all[i] = b * (result.x[0]) * np.cos(result.x[0])
            y_all[i] = b * (result.x[0]) * np.sin(result.x[0])
            theta_all[i] = result.x[0]

    print("---------2. 寻找龙头外圈区域的两个龙身---------")

    theta_d = theta_0 + 2 * (np.pi)
    # 寻找索引
    (index_a, index_b) = find_range(theta_all, theta_d)
    x_a = b * theta_all[index_a] * np.cos(theta_all[index_a])
    x_b = b * theta_all[index_b] * np.cos(theta_all[index_b])
    y_a = b * theta_all[index_a] * np.sin(theta_all[index_a])
    y_b = b * theta_all[index_b] * np.sin(theta_all[index_b])

    # 计算此两点向下15cm的直角坐标系直线
    _k = (y_b - y_a) / (x_b - x_a)
    cos_beta = (np.abs(x_b - x_a)) / (
        np.sqrt((y_b - y_a) * (y_b - y_a) + (x_b - x_a) * (x_b - x_a))
    )
    _b = -_k * x_a + y_a - 0.15 / cos_beta

    print("---------3. 计算外点坐标---------")

    # 计算外点坐标
    phi = np.arcsin((b * theta_all[1] * np.sin(theta_all[1] - theta_all[0])) / 2.86)
    psi = 2 * np.pi - (np.pi) / 2 - np.arctan(27.5 / 15) - phi
    r_c = np.sqrt(
        (b * theta_all[0]) * (b * theta_all[0])
        + 0.098125
        - 2 * (b * theta_all[0]) * np.sqrt(0.098125) * np.cos(psi)
    )
    theta_c = theta_all[0] - np.arcsin(np.sqrt(0.098125) * np.sin(psi) / r_c)
    y_c = r_c * np.sin(theta_c)
    x_c = r_c * np.cos(theta_c)

    y_test_c = _k * x_c + _b

    if y_test_c > y_c:
        flag = 0
        theta_0 = theta_0 - delta_theta
        count = count + 1
    else:
        flag = 1


print("----------finish:theta = ------------")

finish_theta0 = theta_0

print("---------1. 计算龙身位置---------")

theta_all[0] = finish_theta0
x_all[0] = b * finish_theta0 * np.cos(finish_theta0)
y_all[0] = b * finish_theta0 * np.sin(finish_theta0)

for i in range(1, 224):
    if i == 1:
        result = opt.root(
            lambda x: theta_solve_func(x, finish_theta0, 2.86), (finish_theta0 + 1)
        )
        x_all[i] = b * (result.x[0]) * np.cos(result.x[0])
        y_all[i] = b * (result.x[0]) * np.sin(result.x[0])
        theta_all[i] = result.x[0]
    else:
        result = opt.root(
            lambda x: theta_solve_func(x, theta_all[i - 1], 1.65),
            (theta_all[i - 1] + 1),
        )
        x_all[i] = b * (result.x[0]) * np.cos(result.x[0])
        y_all[i] = b * (result.x[0]) * np.sin(result.x[0])
        theta_all[i] = result.x[0]

plt.figure(figsize=(5.5, 5.5))
__theta = np.arange(0, (16 * np.pi * 2), 0.0001)
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

plt.scatter(x_c, y_c, linewidth=1, color="green", marker="o", s=4, label="龙头前进方向边界点")
plt.scatter(x_c, y_test_c, linewidth=1, color="purple", marker="o", s=4, label="龙身边界线上一点")

plt.legend()
plt.show()

print("---------2. 计算龙身速度---------")

v_0 = 1  # 龙头速度
v_i = np.zeros(223)


def alpha_cal_func(theta1, theta0, l):
    """
    中间变量 alpha 计算函数
    """
    return np.arcsin((b * theta1 * np.sin(theta1 - theta0)) / (l))


def beta_cal_func(theta1, theta0, l):
    """
    中间变量 alpha 计算函数
    """
    return np.arcsin((b * theta0 * np.sin(theta1 - theta0)) / (l))


for i in range(0, 223):
    k = i + 1
    p = (int)(i - 1)
    if i == 0:
        v_i[i] = (
            v_0
            * np.abs(
                np.cos(
                    alpha_cal_func(theta_all[k], theta_all[i], 2.86)
                    + (np.pi) / 2
                    - np.arctan(1 / theta_all[i])
                )
            )
            / np.abs(
                np.cos(
                    -beta_cal_func(theta_all[k], theta_all[i], 2.86)
                    + (np.pi) / 2
                    - np.arctan(1 / theta_all[k])
                )
            )
        )
    else:
        v_i[i] = (
            v_i[p]
            * np.abs(
                np.cos(
                    alpha_cal_func(theta_all[k], theta_all[i], 1.65)
                    + (np.pi) / 2
                    - np.arctan(1 / theta_all[i])
                )
            )
            / np.abs(
                np.cos(
                    -beta_cal_func(theta_all[k], theta_all[i], 1.65)
                    + (np.pi) / 2
                    - np.arctan(1 / theta_all[k])
                )
            )
        )

v_all = np.zeros(224)
v_all[0] = v_0
for i in range(0, 223):
    v_all[i + 1] = v_i[i]

theta_init = 16 * 2 * (np.pi)

t = b * (
    0.5 * theta_init * np.sqrt(1 + theta_init * theta_init)
    - 0.5 * np.log(theta_init + np.sqrt(1 + theta_init * theta_init))
) - b * (
    0.5 * finish_theta0 * np.sqrt(1 + finish_theta0 * finish_theta0)
    - 0.5 * np.log(finish_theta0 + np.sqrt(1 + finish_theta0 * finish_theta0))
)

df = pd.DataFrame({"x": x_all, "y": y_all, "v": v_all})
csv_file_path = "x_y_v_data.csv"
df.to_csv(csv_file_path, index=False)

print(t)
