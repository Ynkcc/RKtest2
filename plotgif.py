import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# 设置参数
num_frames = 999  # 总帧数
fps = 50  # 帧速率

# 循环生成每一帧图像
frames = []
for i in range(num_frames):
    # 读取CSV文件
    filename = f'{i}.csv'
    data = np.genfromtxt(filename, delimiter=',')

    # 生成x轴和y轴的坐标值
    x = np.linspace(0, 10, data.shape[0])
    y = np.linspace(0, 10, data.shape[1])
    X, Y = np.meshgrid(x, y)

    # 绘制等高线图
    plt.contourf(X, Y, data, cmap='coolwarm', levels=60, vmin=0, vmax=100)
    plt.colorbar()

    # 将当前帧的图像保存为PIL图像对象
    fig = plt.gcf()
    fig.canvas.draw()
    frame = Image.frombytes('RGB', fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
    frames.append(frame)

    # 清除当前帧的图像
    plt.clf()

# 保存图像序列为GIF动画文件
frames[0].save('animation.gif', save_all=True, append_images=frames[1:], optimize=False, duration=int(1000 / fps), loop=0)
