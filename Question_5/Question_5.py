"""
    问题4_2,问题5代码
    取标准单位,长度为m
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

plt.rcParams["font.sans-serif"] = ["SimHei"]  # 设置字体
plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题

b = 1.7 / (np.pi * 2)

T = 0.00001
v_count = 0  # 速度最大值计数
v0_max = 1  # 速度最大值初始化
v0_min = 1  # 速度最小值初始化

while True:
    # 二分更新
    v_0 = (v0_max + v0_min) / 2

    print("----------------1. 求解圆弧参数-----------------")

    ## 1. 求解圆弧角参数

    LENGTH_K = int(20000)  # k 遍历采样
    LENGTH_R1 = int(100)  # R1 遍历采样
    LENGTH_R0 = int(100)  # R0 遍历采样

    k_sample = np.linspace(0, 1, LENGTH_K)
    r0_sample = np.linspace(1.705, 4.5, LENGTH_R0)

    theta_0 = r0_sample[94] / b  # 进入角
    r1_sample = np.linspace(
        (r0_sample[94] - (np.pi / 2) * b),
        (r0_sample[94] + (np.pi / 2) * b),
        LENGTH_R1,
    )
    theta_1 = r1_sample[76] / b  # 离开角

    ## 2. 进入点和退出点计算
    x_0 = b * theta_0 * np.cos(theta_0)
    y_0 = b * theta_0 * np.sin(theta_0)
    x_1 = -b * theta_1 * np.cos(theta_1)
    y_1 = -b * theta_1 * np.sin(theta_1)

    ## 3. 拟合求圆心
    for i in range(0, LENGTH_K):
        k = k_sample[LENGTH_K - i - 1]
        point_1_x = k * x_1
        point_1_y = k * y_1
        point_0_x = (
            (r0_sample[94] - 2 * (1 - k) * r1_sample[76]) / (r0_sample[94])
        ) * x_0
        point_0_y = (
            (r0_sample[94] - 2 * (1 - k) * r1_sample[76]) / (r0_sample[94])
        ) * y_0
        d = ((point_0_x - point_1_x) ** 2) + ((point_0_y - point_1_y) ** 2)
        d_test = (3 * (1 - k) * r1_sample[76]) ** 2
        if np.abs(d - d_test) < 1e-2:
            # print(k)
            break

    ## 4. 计算切换点
    point_x_temp = point_1_x + (point_0_x - point_1_x) / 3
    point_y_temp = point_1_y + (point_0_y - point_1_y) / 3

    ## 5. 计算圆弧半径
    out_r = (1 - k) * r1_sample[76]
    in_r = 2 * out_r

    ## 6. 计算圆心角
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
    if r0_sample[94] > r1_sample[76]:
        _in_theta = np.pi * 2 - in_theta
        l_in = _in_theta * in_r
        l_out = out_theta * out_r
    else:
        _out_theta = np.pi * 2 - out_theta
        l_in = in_theta * in_r
        l_out = _out_theta * out_r

    if r0_sample[94] > r1_sample[76]:
        in_theta = np.pi * 2 - in_theta
    else:
        out_theta = np.pi * 2 - out_theta

    print("----------------2. 构建t(s)对坐标的映射函数-----------------")

    ## 1. 求首段函数
    def theta_1_solve_func(x, t):
        """
        t对theta的映射,注意时间零点
        """
        return (
            b
            * (
                0.5 * theta_0 * np.sqrt(1 + theta_0 * theta_0)
                - 0.5 * np.log(theta_0 + np.sqrt(1 + theta_0 * theta_0))
            )
            - b * (0.5 * x * np.sqrt(1 + x * x) - 0.5 * np.log(x + np.sqrt(1 + x * x)))
        ) - v_0 * t

    def theta_2_solve_func(x, __s):
        """
        t对theta的映射,注意时间零点
        """
        return (
            b
            * (
                0.5 * theta_1 * np.sqrt(1 + theta_1 * theta_1)
                - 0.5 * np.log(theta_1 + np.sqrt(1 + theta_1 * theta_1))
            )
            - b * (0.5 * x * np.sqrt(1 + x * x) - 0.5 * np.log(x + np.sqrt(1 + x * x)))
        ) + __s

    ## 2. 求分段时间
    t_temp = l_in / v_0
    t_out = (l_in + l_out) / v_0

    def get_coordinate(t):
        """
        t对坐标的映射,返回当前坐标值
        """
        if t <= 0:
            _result = opt.root(theta_1_solve_func, 0, t)
            _theta = _result.x[0]
            _x = b * _theta * np.cos(_theta)
            _y = b * _theta * np.sin(_theta)
            return _x, _y
        elif t > 0 and t <= t_temp:
            _rotate_theta = t * v_0 / in_r
            # 计算旋转初态向量
            _temp_rotate_x = x_0 - point_0_x
            _temp_rotate_y = y_0 - point_0_y
            # 旋转(顺时针)
            _temp_rotated_x = _temp_rotate_x * np.cos(
                _rotate_theta
            ) + _temp_rotate_y * np.sin(_rotate_theta)
            _temp_rotated_y = -_temp_rotate_x * np.sin(
                _rotate_theta
            ) + _temp_rotate_y * np.cos(_rotate_theta)
            # 求坐标值
            _x = point_0_x + _temp_rotated_x
            _y = point_0_y + _temp_rotated_y
            return _x, _y
        elif t > t_temp and t <= t_out:
            _rotate_theta = (t - t_temp) * v_0 / out_r
            # 计算旋转初态向量
            _temp_rotate_x = point_x_temp - point_1_x
            _temp_rotate_y = point_y_temp - point_1_y
            # 旋转(顺时针)
            _temp_rotated_x = _temp_rotate_x * np.cos(
                _rotate_theta
            ) - _temp_rotate_y * np.sin(_rotate_theta)
            _temp_rotated_y = _temp_rotate_x * np.sin(
                _rotate_theta
            ) + _temp_rotate_y * np.cos(_rotate_theta)
            # 求坐标值
            _x = point_1_x + _temp_rotated_x
            _y = point_1_y + _temp_rotated_y
            return _x, _y
        else:
            _s = t * v_0 - l_in - l_out
            _result = opt.root(theta_2_solve_func, 0, _s)
            _theta = _result.x[0]
            _x = -b * _theta * np.cos(_theta)
            _y = -b * _theta * np.sin(_theta)
            return _x, _y

    def get_trail_coordinate(__s):
        return get_coordinate(__s / v_0)

    def get_all_coord(t):
        all_coord = get_coordinate(t)
        all_x = [all_coord[0]]
        all_y = [all_coord[1]]
        all_s = [t * v_0]

        def equation1(__s):
            x_next, y_next = get_trail_coordinate(__s)
            return (x_next - all_x[0]) ** 2 + (y_next - all_y[0]) ** 2 - 2.86**2

        initial_guess = all_s[0] - 2.86

        solution = opt.fsolve(equation1, initial_guess)
        all_s.append(solution[0])
        all_x.append(get_trail_coordinate(all_s[1])[0])
        all_y.append(get_trail_coordinate(all_s[1])[1])
        for i in range(222):

            def equation2(__s):
                x_next, y_next = get_trail_coordinate(__s)
                return (
                    (x_next - all_x[i + 1]) ** 2
                    + (y_next - all_y[i + 1]) ** 2
                    - 1.65**2
                )

            initial_guess = all_s[i + 1] - 1.65
            solution = opt.fsolve(equation2, initial_guess)
            all_s.append(solution[0])
            all_x.append(get_trail_coordinate(all_s[i + 2])[0])
            all_y.append(get_trail_coordinate(all_s[i + 2])[1])

        return [all_x, all_y]

    print("----------------3. 画图-----------------")

    # in_theta_reg = in_theta * 180 / np.pi
    # out_theta_reg = out_theta * 180 / np.pi

    # in_theta_start = (
    #     np.arctan2((y_0 - point_0_y), (x_0 - point_0_x)) * 180 / np.pi + 360
    # )
    # out_theta_start = np.arctan2((y_1 - point_1_y), (x_1 - point_1_x)) * 180 / np.pi
    # in_theta_end = (
    #     np.arctan2((point_y_temp - point_0_y), (point_x_temp - point_0_x)) * 180 / np.pi
    # )
    # out_theta_end = (
    #     np.arctan2((point_y_temp - point_1_y), (point_x_temp - point_1_x)) * 180 / np.pi
    # )

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

    # plt.figure(figsize=(10, 10))
    # __theta = np.arange(0, (10 * np.pi * 2), 0.0001)
    # __x_in = b * (__theta) * np.cos(__theta)
    # __y_in = b * (__theta) * np.sin(__theta)
    # __x_out = -__x_in
    # __y_out = -__y_in
    # plt.plot(__x_out, __y_out, linewidth=1, color="red")
    # plt.plot(__x_in, __y_in, linewidth=1, color="red", label="螺线")
    # X0, Y0 = circle_array(point_0_x, point_0_y, in_r, in_theta_end, in_theta_start)
    # X1, Y1 = circle_array(point_1_x, point_1_y, out_r, out_theta_end, out_theta_start)
    # X2, Y2 = circle_array(0, 0, 4.5, 0, 360)
    # plt.plot(X0, Y0, linewidth=1, color="blue")
    # plt.plot(X1, Y1, linewidth=1, color="blue", label="调头曲线")
    # plt.plot(X2, Y2, linestyle="dashed", linewidth=1, color="purple", label="调头空间")
    # plt.scatter(x_0, y_0, s=2)
    # plt.scatter(x_1, y_1, s=2)
    # plt.scatter(point_0_x, point_0_y, s=2)
    # plt.scatter(point_1_x, point_1_y, s=2)
    # plt.scatter(point_x_temp, point_y_temp, s=2)
    # plt.plot(get_all_coord(100)[0], get_all_coord(100)[1], color="yellow", label="龙身")
    # plt.scatter(
    #     get_all_coord(100)[0],
    #     get_all_coord(100)[1],
    #     linewidth=1,
    #     color="blue",
    #     marker="o",
    #     s=4,
    # )

    # plt.legend()
    # plt.show()

    print("----------------4. 数据处理-----------------")

    all_coordination = np.zeros((448, 201))

    for j in range(0, 201):
        print(j)
        # 时间标记
        m = j - 100
        x_m = get_all_coord(m)[0]
        y_m = get_all_coord(m)[1]
        for i in range(0, 448):
            n1 = (int)(i / 2)
            n2 = (int)((i - 1) / 2)
            if i % 2 == 0:
                all_coordination[i][j] = x_m[n1]
            else:
                all_coordination[i][j] = y_m[n2]

    # df = pd.DataFrame(all_coordination)
    # csv_x_y_file_path = "coordination_data.csv"
    # df.to_csv(csv_x_y_file_path, index=False)

    print("---finish0!---")

    all_coordination1 = np.zeros((448, 201))

    for j in range(0, 201):
        print(j)
        # 时间标记
        m = j - 100 + T
        x_m = get_all_coord(m)[0]
        y_m = get_all_coord(m)[1]
        for i in range(0, 448):
            n1 = (int)(i / 2)
            n2 = (int)((i - 1) / 2)
            if i % 2 == 0:
                all_coordination1[i][j] = x_m[n1]
            else:
                all_coordination1[i][j] = y_m[n2]

    print("---finish1!---")

    delta_coordination = np.zeros((448, 201))

    ## 两个表进行减法
    for i in range(0, 448):
        for j in range(0, 201):
            delta_coordination[i][j] = all_coordination1[i][j] - all_coordination[i][j]

    v_i = np.zeros((224, 201))

    ## 求速度
    for i in range(0, 224):
        k = 2 * i
        p = 2 * i + 1
        for j in range(0, 201):
            v_i[i][j] = (
                np.sqrt(
                    delta_coordination[k][j] * delta_coordination[k][j]
                    + delta_coordination[p][j] * delta_coordination[p][j]
                )
                / T
            )

    X = np.arange(0, 223)
    Y = np.arange(-100, 100)
    Z = v_i

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap="viridis", edgecolor="none")
    ax.set_xlabel("龙身节数")
    ax.set_ylabel("时间")
    ax.set_zlabel("速度")
    ax.set_title("龙身的速度-时间关系")

    plt.show()

    # ## 速度更新
    # if v_i.max() > 2:
    #     v0_max = v_0
    # else:
    #     v0_min = v_0

    v_count = v_count + 1
    # print("%d次二分,v最大值%f", (v_count, v_i.max()))

    if v_count == 1:
        break

# print("ok!")
# print(v_0)

# df = pd.DataFrame(v_i)
# csv_vi_file_path = 'vi_data.csv'
# df.to_csv(csv_vi_file_path, index=False)
