import numpy as np
import time

class HeatConductionEquationSolver:
    def __init__(self, L=10, T=2):
        self.Lx = L
        self.Ly = L
        self.T = T
        self.dx = 0.1
        self.dy = 0.1
        self.dt = 0.001
        self.alpha = 1
        self.Nx = int(self.Lx / self.dx)
        self.Ny = int(self.Ly / self.dy)
        self.tn = round(self.T / self.dt)
        self.u = np.zeros((self.Nx, self.Ny))
        for i in range(self.Nx):
            self.u[i, 0] = self.u[i, self.Ny - 1] = 100
        for j in range(self.Ny):
            self.u[0, j] = self.u[self.Nx - 1, j] = 100
    
    def deriv_center(self, f):
        ddf = np.zeros((self.Nx, self.Ny))
        for i in range(1, self.Nx - 1):
            for j in range(1, self.Ny - 1):
                ddf[i, j] = (f[i + 1, j] + f[i - 1, j] - 2 * f[i, j]) / (self.dx * self.dx) + (f[i, j + 1] + f[i, j - 1] - 2 * f[i, j]) / (self.dy * self.dy)
        return self.alpha * ddf
    
    def solve(self):
        for t in range(1, self.tn):
            K1 = self.deriv_center(self.u[:, :])
            K2 = self.deriv_center(self.u[:, :] + 0.5 * self.dt * K1)
            K3 = self.deriv_center(self.u[:, :] + 0.5 * self.dt * K2)
            K4 = self.deriv_center(self.u[:, :] + self.dt * K3)
            self.u[:, :] = self.u[:, :] + self.dt * (1.0 / 6 * K1 + 1.0 / 3 * K2 + 1.0 / 3 * K3 + 1.0 / 6 * K4)

if __name__ == '__main__':
    start_time = time.time()


    simulate = HeatConductionEquationSolver()
    simulate.solve()
    # 程序代码

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"程序运行时间为：{elapsed_time:.2f}秒")
