#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    // 设置参数
    double alpha = 1.0;  // 热扩散系数
    double dt = 0.0001;  // 时间步长
    double dx = 0.1;  // 空间步长
    double L = 1.0;  // 区域长度
    int N = L / dx;  // 空间网格数
    double T = 1.0;  // 模拟总时间
    int M = T / dt;  // 时间网格数
    double r = alpha * dt / dx / dx;  // 参数r

    // 初始化矩阵
    mat U(N + 1, N + 1, fill::zeros);  // 初始温度分布
    U(N / 2, N / 2) = 100.0;  // 在正中央加入热源

    // 迭代计算
    for (int k = 0; k < M; k++)
    {
        // 计算下一个时间步的温度分布
        for (int i = 1; i < N; i++)
        {
            for (int j = 1; j < N; j++)
            {
                U(i, j) = U(i, j) + r * (U(i - 1, j) - 2 * U(i, j) + U(i + 1, j)) + r * (U(i, j - 1) - 2 * U(i, j) + U(i, j + 1));
            }
        }

        // 输出当前温度分布
        //cout << "Time step " << k << ":" << endl;
        //cout << U << endl;
    }
    U.save("u.csv", csv_ascii);
    return 0;
}
