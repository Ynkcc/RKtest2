import numpy as np
import matplotlib.pyplot as plt
# 设置参数
Lx = 1.0   # 空间区间长度
Ly = 1.0
Nx = 101   # 空间网格数
Ny = 101
dx = Lx / (Nx - 1)   # 空间步长
dy = Ly / (Ny - 1)
T0 = 0.0   # 初始温度
T1 = 100.0   # 边界温度
alpha = 1   # 热扩散系数
dt = 0.001   # 时间步长
tmax = 1.0   # 最大时间

# 初始化温度场
T = T0 * np.ones((Nx, Ny))
T[0, :] = T1
T[Nx-1, :] = T1
T[:, 0] = T1
T[:, Ny-1] = T1

# 定义中心差分函数
def center_diff(u, h):
    return (u[2:] - 2*u[1:-1] + u[:-2]) / h**2

# 迭代求解温度场
t = 0.0
while t < tmax:
    # 使用4阶中心差分法计算空间偏导数
    Tx = center_diff(T, dx)
    Ty = center_diff(T.T, dy).T
    
        # 使用龙格库塔法迭代更新
    T1 = T + 0.5 * dt * alpha * (Tx[1:-1, :] + Tx[:-2, :] + Ty[:, 1:-1] + Ty[:, :-2])
    Tx1 = center_diff(T1, dx)
    Ty1 = center_diff(T1.T, dy).T
    T2 = T + 0.5 * dt * alpha * (Tx1[1:-1, :] + Tx1[:-2, :] + Ty1[:, 1:-1] + Ty1[:, :-2])
    Tx2 = center_diff(T2, dx)
    Ty2 = center_diff(T2.T, dy).T
    T3 = T + dt * alpha * (Tx2[1:-1, :] + Tx2[:-2, :] + Ty2[:, 1:-1] + Ty2[:, :-2])

    # 更新温度场
    T = T + dt/6 * (Tx[1:-1, :] + 2*Tx1[1:-1, :] + 2*Tx2[1:-1, :] + center_diff(T3, dx))

    # 更新时间
    t += dt


    # 绘制热图
    plt.imshow(T, cmap='hot', extent=[0, Lx, 0, Ly], origin='lower')
    plt.colorbar()

    # 设置坐标轴标签
    plt.xlabel('x')
    plt.ylabel('y')

    # 显示图形
    plt.show()

print(T)