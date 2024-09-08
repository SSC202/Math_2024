"""
    问题4代码
    取标准单位,长度为m
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

plt.rcParams["font.sans-serif"] = ["SimHei"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题

b = 1.7 / (np.pi * 2)

# plt.figure(figsize=(5.5, 5.5))
# __theta = np.arange(0, (4 * np.pi * 2), 0.0001)
# __x_in = b * (__theta) * np.cos(__theta)
# __y_in = b * (__theta) * np.sin(__theta)
# __x_out = -__x_in
# __y_out = -__y_in
# plt.plot(__x_out, __y_out, linewidth=1, color="red")
# plt.plot(__x_in, __y_in, linewidth=1, color="red", label="螺线")


LENGTH_R0 = int(100)  # R0 遍历采样
LENGTH_R1 = int(100)  # R1 遍历采样
LENGTH_K = int(20000)  # k 遍历采样

path_length = np.zeros((LENGTH_R0, LENGTH_R1))

r0_sample = np.linspace(1.705, 4.5, LENGTH_R0)
k_sample = np.linspace(0, 1, LENGTH_K)

for i_r0 in range(0, LENGTH_R0):

    # 生成 r1 采样点
    r1_sample = np.linspace(
        (r0_sample[i_r0] - (np.pi / 2) * b),
        (r0_sample[i_r0] + (np.pi / 2) * b),
        LENGTH_R1,
    )

    for i_r1 in range(0, LENGTH_R1):
        theta_0 = r0_sample[i_r0] / b
        theta_1 = r1_sample[i_r1] / b
        # 进入点和退出点计算
        x_0 = b * theta_0 * np.cos(theta_0)
        y_0 = b * theta_0 * np.sin(theta_0)
        x_1 = -b * theta_1 * np.cos(theta_1)
        y_1 = -b * theta_1 * np.sin(theta_1)
        # 拟合求圆心
        for i in range(0, LENGTH_K):
            k = k_sample[LENGTH_K - i - 1]
            point_1_x = k * x_1
            point_1_y = k * y_1
            point_0_x = (
                (r0_sample[i_r0] - 2 * (1 - k) * r1_sample[i_r1]) / (r0_sample[i_r0])
            ) * x_0
            point_0_y = (
                (r0_sample[i_r0] - 2 * (1 - k) * r1_sample[i_r1]) / (r0_sample[i_r0])
            ) * y_0
            d = ((point_0_x - point_1_x) ** 2) + ((point_0_y - point_1_y) ** 2)
            d_test = (3 * (1 - k) * r1_sample[i_r1]) ** 2
            if np.abs(d - d_test) < 1e-2:
                # print(k)
                break
        # 求弧长
        out_r = (1 - k) * r1_sample[i_r1]
        in_r = 2 * out_r

        d_in = ((x_0 - point_1_x) ** 2) + (
            (y_0 - point_1_y) ** 2
        )  # 进入点和小圆圆心距离平方
        d_out = ((x_1 - point_0_x) ** 2) + (
            (y_1 - point_0_y) ** 2
        )  # 出点和大圆圆心距离平方

        temp_in = ((4 * out_r * out_r) + d - d_in) / (4 * out_r * np.sqrt(d))
        if temp_in > 1:
            temp_in = 1
        elif temp_in < -1:
            temp_in = -1

        temp_out = ((out_r * out_r) + d - d_out) / (2 * out_r * np.sqrt(d))
        if temp_out > 1:
            temp_out = 1
        elif temp_out < -1:
            temp_out = -1

        in_theta = np.arccos(temp_in)
        out_theta = np.arccos(temp_out)
        if r0_sample[i_r0] > r1_sample[i_r1]:
            _in_theta = np.pi * 2 - in_theta
            l = _in_theta * in_r + out_theta * out_r
        else:
            _out_theta = np.pi * 2 - out_theta
            l = in_theta * in_r + _out_theta * out_r

        # if out_r < 1.43:
        #     l = 0

        # if r0_sample[i_r0] + r1_sample[i_r1] > 9:
        #     l = 0

        path_length[i_r0][i_r1] = l

        print(l)

Z = path_length
X = np.linspace(1.705, 4.5, LENGTH_R0)
Y = np.linspace(
    (-(np.pi / 2) * b),
    (+(np.pi / 2) * b),
    LENGTH_R1,
)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap="viridis", edgecolor="none")
ax.set_xlabel("r0")
ax.set_ylabel("theta")
ax.set_zlabel("path_length")
ax.set_title("未加入约束条件的路径长度，起点坐标和终点坐标关系")

plt.show()


# df = pd.DataFrame(path_length)
# csv_file_path = 'path_data.csv'
# df.to_csv(csv_file_path, index=False)


# theta_0 = r0_sample[94] / b
# r1_sample = np.linspace(
#     (r0_sample[94] - (np.pi / 2) * b),
#     (r0_sample[94] + (np.pi / 2) * b),
#     LENGTH_R1,
# )
# theta_1 = r1_sample[76] / b


# # 进入点和退出点计算
# x_0 = b * theta_0 * np.cos(theta_0)
# y_0 = b * theta_0 * np.sin(theta_0)
# x_1 = -b * theta_1 * np.cos(theta_1)
# y_1 = -b * theta_1 * np.sin(theta_1)
# # 拟合求圆心
# for i in range(0, LENGTH_K):
#     k = k_sample[LENGTH_K - i - 1]
#     point_1_x = k * x_1
#     point_1_y = k * y_1
#     point_0_x = ((r0_sample[94] - 2 * (1 - k) * r1_sample[76]) / (r0_sample[94])) * x_0
#     point_0_y = ((r0_sample[94] - 2 * (1 - k) * r1_sample[76]) / (r0_sample[94])) * y_0
#     d = ((point_0_x - point_1_x) ** 2) + ((point_0_y - point_1_y) ** 2)
#     d_test = (3 * (1 - k) * r1_sample[76]) ** 2
#     if np.abs(d - d_test) < 1e-2:
#         # print(k)
#         break

# # print(r0_sample[96],r1_sample[76])

# # 计算切换点
# point_x_temp = point_1_x + (point_0_x - point_1_x) / 3
# point_y_temp = point_1_y + (point_0_y - point_1_y) / 3

# out_r = (1 - k) * r1_sample[76]
# in_r = 2 * out_r

# d_in = ((x_0 - point_1_x) ** 2) + ((y_0 - point_1_y) ** 2)  # 进入点和小圆圆心距离平方
# d_out = ((x_1 - point_0_x) ** 2) + ((y_1 - point_0_y) ** 2)  # 出点和大圆圆心距离平方

# temp_in = ((4 * out_r * out_r) + d - d_in) / (4 * out_r * np.sqrt(d))
# if temp_in > 1:
#     temp_in = 1
# elif temp_in < -1:
#     temp_in = -1

# temp_out = ((out_r * out_r) + d - d_out) / (2 * out_r * np.sqrt(d))
# if temp_out > 1:
#     temp_out = 1
# elif temp_out < -1:
#     temp_out = -1

# in_theta = np.arccos(temp_in)
# out_theta = np.arccos(temp_out)

# if r0_sample[94] > r1_sample[76]:
#     in_theta = np.pi * 2 - in_theta
# else:
#     out_theta = np.pi * 2 - out_theta

# in_theta_reg = in_theta * 180 / np.pi
# out_theta_reg = out_theta * 180 / np.pi

# in_theta_start = np.arctan2((y_0 - point_0_y), (x_0 - point_0_x)) * 180 / np.pi + 360
# out_theta_start = np.arctan2((y_1 - point_1_y), (x_1 - point_1_x)) * 180 / np.pi
# in_theta_end = (
#     np.arctan2((point_y_temp - point_0_y), (point_x_temp - point_0_x)) * 180 / np.pi
# )
# out_theta_end = (
#     np.arctan2((point_y_temp - point_1_y), (point_x_temp - point_1_x)) * 180 / np.pi
# )

# # print(in_theta_start, out_theta_start)
# # print(in_theta_end, out_theta_end)


# def circle_array(xc, yc, r, start, end):
#     # 根据圆心、半径、起始角度、结束角度，生成圆弧的数据点
#     phi1 = start * np.pi / 180.0
#     phi2 = end * np.pi / 180.0
#     dphi = (phi2 - phi1) / np.ceil(
#         200 * np.pi * r * (phi2 - phi1)
#     )  # 根据圆弧周长设置等间距
#     # array = np.arange(phi1, phi2+dphi, dphi) #生成双闭合区间(Python默认为左闭右开)
#     # 当dphi为无限小数时，上述方法可能由于小数精度问题，导致array多一个元素，将代码优化为下面两行
#     array = np.arange(phi1, phi2, dphi)
#     array = np.append(array, array[-1] + dphi)  # array在结尾等距增加一个元素
#     return xc + r * np.cos(array), yc + r * np.sin(array)


# X0, Y0 = circle_array(point_0_x, point_0_y, in_r, in_theta_end, in_theta_start)
# X1, Y1 = circle_array(point_1_x, point_1_y, out_r, out_theta_end, out_theta_start)
# X2, Y2 = circle_array(0, 0, 4.5, 0, 360)
# plt.plot(X0, Y0, linewidth=1, color="blue")
# plt.plot(X1, Y1, linewidth=1, color="blue")
# plt.plot(X2, Y2, linestyle="dashed", linewidth=1, color="purple")
# plt.scatter(x_0, y_0, s=2)
# plt.scatter(x_1, y_1, s=2)
# plt.scatter(point_0_x, point_0_y, s=2)
# plt.scatter(point_1_x, point_1_y, s=2)
# plt.scatter(point_x_temp, point_y_temp, s=2)

# print(x_0, y_0)
# print(point_0_x, point_0_y)
# print(point_1_x, point_1_y)
# print(point_x_temp, point_y_temp)
# print(in_r, out_r)

# plt.show()
