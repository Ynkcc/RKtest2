#include <armadillo>
#include <chrono>
using namespace arma;

class HeatConductionEquationSolver {
public:
	HeatConductionEquationSolver(int L = 10, int T = 20) :Lx(L), Ly(L), T(T) {

		dx = dy = 0.1;
		dt = 0.001;
		alpha = 1;
		Nx = Ny = Lx / dx;
		tn = round(T / dt);
		u = zeros(Nx, Ny);

		for (int i = 0; i < Nx; i++)
		{
			u(i, 0) = u(i, Ny - 1) = 100;
		}
		for (int j = 0; j < Ny; j++)
		{
			u(0, j) = u(Nx - 1, j) = 100;
		}
		u.save(  "0.csv", csv_ascii);
	}





	void solve() {
		for (int t = 1; t < tn; t++) {
			mat K1 = deriv_center(u);
			mat K2 = deriv_center(u + 0.5 * dt * K1);
			mat K3 = deriv_center(u + 0.5 * dt * K2);
			mat K4 = deriv_center(u + dt * K3);
			u +=  dt * (1.0 / 6 * K1 + 1.0 / 3 * K2 + 1.0 / 3 * K3 + 1.0 / 6 * K4);
			if (t%20==0)
			{
				u.save(std::to_string(t / 20)+".csv" , csv_ascii);
			}
		}
	}



private:

	int Nx, Ny, T;
	double Lx, Ly;
	double dt;
	double alpha;
	int tn;
	double dx, dy;
	mat u;
	vec x;
	vec t;
	vec error;

	mat deriv_center(const mat& f) {
		mat ddf(Nx, Ny);
		for (int i = 1; i < Nx - 1; i++) {
			for (int j = 1; j < Ny - 1; j++)
			{
				ddf(i, j) = (f(i + 1, j) + f(i - 1, j) - 2 * f(i, j)) / (dx * dx) + (f(i, j + 1) + f(i, j - 1) - 2 * f(i, j)) / (dy * dy);
			}
		}

		return alpha * ddf;
	}
};
int main() {
	auto start = std::chrono::steady_clock::now(); // 记录开始时间
	auto simulate = HeatConductionEquationSolver();
	simulate.solve();
	auto end = std::chrono::steady_clock::now(); // 记录结束时间
	auto diff = end - start; // 计算时间差
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(diff); // 转换为毫秒
	std::cout << "Time taken: " << duration.count() << " milliseconds." << std::endl;
	system("pause");
	return 0;
}