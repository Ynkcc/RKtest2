import numpy as np
import matplotlib.pyplot as plt

# 读取CSV文件
data = np.genfromtxt('600.csv', delimiter=',')

# 生成x轴和y轴的坐标值
x = np.linspace(0, 10, data.shape[0])
y = np.linspace(0, 10, data.shape[1])
X, Y = np.meshgrid(x, y)

# 绘制等高线图
plt.contourf(X, Y, data, cmap='coolwarm',levels=60)
plt.colorbar()
plt.show()