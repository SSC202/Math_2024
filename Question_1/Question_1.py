"""
    问题1代码
    取标准单位,长度为m
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

plt.rcParams["font.sans-serif"] = ["SimHei"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题

T = 1e-7

# 等距螺线方程为：r = bθ

b = 0.55 / (2 * np.pi)  # 解螺距参数
theta_init = 16 * 2 * np.pi  # 初始角度
theta_all = np.zeros((224, 301))

print("---------0. 解龙头位置---------")


def theta_0_solve_func(x, t):
    """
    解龙头位置方程函数
    """
    return (
        b
        * (
            0.5 * theta_init * np.sqrt(1 + theta_init * theta_init)
            - 0.5 * np.log(theta_init + np.sqrt(1 + theta_init * theta_init))
        )
        - b * (0.5 * x * np.sqrt(1 + x * x) - 0.5 * np.log(x + np.sqrt(1 + x * x)))
    ) - t


theta_0 = []
x_0 = []
y_0 = []

# 0. 解龙头位置
for i in range(0, 301):
    result_0 = opt.root(theta_0_solve_func, 0, i)
    x_0.append(b * (result_0.x[0]) * np.cos(result_0.x[0]))
    y_0.append(b * (result_0.x[0]) * np.sin(result_0.x[0]))
    theta_0.append(result_0.x[0])
    # print(result_0.x[0])

# df = pd.DataFrame({'x': x_0, 'y': y_0})
# csv0_file_path = 'x0_y0_data.csv'
# df.to_csv(csv0_file_path, index=False)
for i in range(0, 301):
    theta_all[0][i] = theta_0[i]

print("---------1. 计算龙身位置---------")


# 1. 计算龙身位置
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


theta_1 = []
x_1 = []
y_1 = []

### 第一段，龙头位置较长
for i in range(0, 301):
    result_1 = opt.root(
        lambda x: theta_solve_func(x, theta_0[i], 2.86), (theta_0[i] + 1)
    )
    x_1.append(b * (result_1.x[0]) * np.cos(result_1.x[0]))
    y_1.append(b * (result_1.x[0]) * np.sin(result_1.x[0]))
    theta_1.append(result_1.x[0])
    # print(result_1.x[0])

for i in range(0, 301):
    theta_all[1][i] = theta_1[i]

# df = pd.DataFrame({'x': x_1, 'y': y_1})
# csv1_file_path = 'x1_y1_data.csv'
# df.to_csv(csv1_file_path, index=False)

### 第二段，从第2段龙身到第222段龙身

x_i = np.ones((222, 301))
y_i = np.ones((222, 301))
theta_i = np.ones((222, 301))

for i in range(0, 222):
    for j in range(0, 301):
        if i == 0:
            result = opt.root(
                lambda x: theta_solve_func(x, theta_1[j], 1.65), (theta_1[j] + 1)
            )
            x_i[i][j] = b * (result.x[0]) * np.cos(result.x[0])
            y_i[i][j] = b * (result.x[0]) * np.sin(result.x[0])
            theta_i[i][j] = result.x[0]
        else:
            k = i - 1
            result = opt.root(
                lambda x: theta_solve_func(x, theta_i[k][j], 1.65), (theta_i[k][j] + 1)
            )
            x_i[i][j] = b * (result.x[0]) * np.cos(result.x[0])
            y_i[i][j] = b * (result.x[0]) * np.sin(result.x[0])
            theta_i[i][j] = result.x[0]

body_coordinate = np.ones((444, 301))

for i in range(0, 444):
    for j in range(0, 301):
        if i % 2 == 0:
            k = int(i / 2)
            body_coordinate[i][j] = x_i[k][j]
        else:
            k = int((i - 1) / 2)
            body_coordinate[i][j] = y_i[k][j]

for j in range(0, 222):
    k = j + 2
    for i in range(0, 301):
        theta_all[k][i] = theta_i[j][i]

# df = pd.DataFrame(body_coordinate)
# csvi_file_path = 'xi_yi_data.csv'
# df.to_csv(csvi_file_path, index=False)

# for i in range(0, 222):
#     print(theta_i[i][0])

# df = pd.DataFrame(theta_all)
# csv_theta_file_path = 'theta_data.csv'
# df.to_csv(csv_theta_file_path, index=False)

plt.figure(figsize=(6.5, 6.5))
__theta = np.arange(0, (22 * np.pi * 2), 0.0001)
__x = b * (__theta) * np.cos(__theta)
__y = b * (__theta) * np.sin(__theta)
plt.plot(__x, __y, linewidth=1, color="red", label="螺线")
__theta_all = np.zeros(224)

for i in range(0, 224):
    __theta_all[i] = theta_all[i][300]

print(theta_all)
__xi = b * (__theta_all) * np.cos(__theta_all)
__yi = b * (__theta_all) * np.sin(__theta_all)

plt.plot(__xi, __yi, linewidth=1, color="blue", label="龙身")
plt.scatter(__xi, __yi, linewidth=1, color="blue", marker="o", s=4)

plt.legend()
plt.show()


v_i = np.zeros((223, 301))

print("---------2. 计算龙身速度:公式解法---------")

v_0 = 1  # 龙头速度


def alpha_cal_func(theta1, theta0, l):
    """
    中间变量 alpha 计算函数
    """
    return np.arcsin((b * theta1 * np.sin(theta1 - theta0)) / (l))


def beta_cal_func(theta1, theta0, l):
    """
    中间变量 beta 计算函数
    """
    return np.arcsin((b * theta0 * np.sin(theta1 - theta0)) / (l))


for i in range(0, 223):
    k = i + 1
    p = (int)(i - 1)
    for j in range(0, 301):
        if i == 0:
            v_i[i][j] = (
                v_0
                * np.abs(
                    np.cos(
                        alpha_cal_func(theta_all[k][j], theta_all[i][j], 2.86)
                        + (np.pi) / 2
                        - np.arctan(1 / theta_all[i][j])
                    )
                )
                / np.abs(
                    np.cos(
                        -beta_cal_func(theta_all[k][j], theta_all[i][j], 2.86)
                        + (np.pi) / 2
                        - np.arctan(1 / theta_all[k][j])
                    )
                )
            )
        else:
            v_i[i][j] = (
                v_i[p][j]
                * np.abs(
                    np.cos(
                        alpha_cal_func(theta_all[k][j], theta_all[i][j], 1.65)
                        + (np.pi) / 2
                        - np.arctan(1 / theta_all[i][j])
                    )
                )
                / np.abs(
                    np.cos(
                        -beta_cal_func(theta_all[k][j], theta_all[i][j], 1.65)
                        + (np.pi) / 2
                        - np.arctan(1 / theta_all[k][j])
                    )
                )
            )

df = pd.DataFrame(v_i)
csv_vi_approximate_file_path = "vi_approximate_data.csv"
df.to_csv(csv_vi_approximate_file_path, index=False)

# print("---------3. 计算龙身速度:差分解法---------")

# all_coordinate_0 = np.zeros((448, 301))
# all_coordinate_1 = np.zeros((448, 301))
# delta_coordinate = np.zeros((448, 301))

# for i in range(0, 448):
#     k = i - 4
#     for j in range(0, 301):
#         if i >= 0 and i <= 3:
#             if i == 0:
#                 all_coordinate_0[i][j] = x_0[j]
#             elif i == 1:
#                 all_coordinate_0[i][j] = y_0[j]
#             elif i == 2:
#                 all_coordinate_0[i][j] = x_1[j]
#             elif i == 3:
#                 all_coordinate_0[i][j] = y_1[j]
#         else:
#             all_coordinate_0[i][j] = body_coordinate[k][j]

# # df = pd.DataFrame(all_coordinate_0)
# # csv_coordinate_0_file_path = 'coordinate_0_data.csv'
# # df.to_csv(csv_coordinate_0_file_path, index=False)

# _theta_0 = []
# _x_0 = []
# _y_0 = []
# _theta_all = np.zeros((224, 301))

# # 0. 解龙头位置
# for i in range(0, 301):
#     _result_0 = opt.root(theta_0_solve_func, 0, (i + T))
#     _x_0.append(b * (_result_0.x[0]) * np.cos(_result_0.x[0]))
#     _y_0.append(b * (_result_0.x[0]) * np.sin(_result_0.x[0]))
#     _theta_0.append(_result_0.x[0])

# for i in range(0, 301):
#     _theta_all[0][i] = _theta_0[i]

# _theta_1 = []
# _x_1 = []
# _y_1 = []

# ### 第一段，龙头位置较长
# for i in range(0, 301):
#     _result_1 = opt.root(
#         lambda x: theta_solve_func(x, _theta_0[i], 2.86), (_theta_0[i] + 1)
#     )
#     _x_1.append(b * (_result_1.x[0]) * np.cos(_result_1.x[0]))
#     _y_1.append(b * (_result_1.x[0]) * np.sin(_result_1.x[0]))
#     _theta_1.append(_result_1.x[0])
#     # print(result_1.x[0])

# for i in range(0, 301):
#     _theta_all[1][i] = _theta_1[i]

# # df = pd.DataFrame({'x': x_1, 'y': y_1})
# # csv1_file_path = 'x1_y1_data.csv'
# # df.to_csv(csv1_file_path, index=False)

# ### 第二段，从第2段龙身到第222段龙身

# _x_i = np.ones((222, 301))
# _y_i = np.ones((222, 301))
# _theta_i = np.ones((222, 301))

# for i in range(0, 222):
#     for j in range(0, 301):
#         if i == 0:
#             _result = opt.root(
#                 lambda x: theta_solve_func(x, _theta_1[j], 1.65), (_theta_1[j] + 1)
#             )
#             _x_i[i][j] = b * (_result.x[0]) * np.cos(_result.x[0])
#             _y_i[i][j] = b * (_result.x[0]) * np.sin(_result.x[0])
#             _theta_i[i][j] = _result.x[0]
#         else:
#             k = i - 1
#             _result = opt.root(
#                 lambda x: theta_solve_func(x, _theta_i[k][j], 1.65),
#                 (_theta_i[k][j] + 1),
#             )
#             _x_i[i][j] = b * (_result.x[0]) * np.cos(_result.x[0])
#             _y_i[i][j] = b * (_result.x[0]) * np.sin(_result.x[0])
#             _theta_i[i][j] = _result.x[0]

# _body_coordinate = np.ones((444, 301))

# for i in range(0, 444):
#     for j in range(0, 301):
#         if i % 2 == 0:
#             k = int(i / 2)
#             _body_coordinate[i][j] = _x_i[k][j]
#         else:
#             k = int((i - 1) / 2)
#             _body_coordinate[i][j] = _y_i[k][j]

# for i in range(0, 448):
#     k = i - 4
#     for j in range(0, 301):
#         if i >= 0 and i <= 3:
#             if i == 0:
#                 all_coordinate_1[i][j] = _x_0[j]
#             elif i == 1:
#                 all_coordinate_1[i][j] = _y_0[j]
#             elif i == 2:
#                 all_coordinate_1[i][j] = _x_1[j]
#             elif i == 3:
#                 all_coordinate_1[i][j] = _y_1[j]
#         else:
#             all_coordinate_1[i][j] = _body_coordinate[k][j]

# # df = pd.DataFrame(all_coordinate_1)
# # csv_coordinate_1_file_path = 'coordinate_1_data.csv'
# # df.to_csv(csv_coordinate_1_file_path, index=False)

# ## 两个表进行减法
# for i in range(0, 448):
#     for j in range(0, 301):
#         delta_coordinate[i][j] = all_coordinate_1[i][j] - all_coordinate_0[i][j]

# ## 求速度
# for i in range(0, 224):
#     k = 2 * i
#     p = 2 * i + 1
#     for j in range(0, 301):
#         v_i[i][j] = (
#             np.sqrt(
#                 delta_coordinate[k][j] * delta_coordinate[k][j]
#                 + delta_coordinate[p][j] * delta_coordinate[p][j]
#             )
#             / T
#         )

# df = pd.DataFrame(v_i)
# csv_vi_file_path = 'vi_data.csv'
# df.to_csv(csv_vi_file_path, index=False)
