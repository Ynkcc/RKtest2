#include <armadillo>
using namespace arma;

class HeatConductionEquationSolver {
public:
	HeatConductionEquationSolver(int L = 10, int T = 10) :Lx(L), Ly(L), T(T) {

		dx = dy = 0.4;
		dt = 0.001;
		alpha = 1;
		Nx = Ny = Lx / dx;
		tn = round(T / dt);
		u = zeros(Nx, Ny, tn);
		for (size_t i = 0; i < Nx; i++)
		{
			for (size_t j = 0; j < Ny; j++)
			{
				u(i, j, 0) = sin(i+j);
			}
		}

		//for (int i = 0; i < Nx; i++)
		//{
		//	u(i, 0, 0) = u(i, Ny - 1, 0) = 100;
		//}
		//for (int j = 0; j < Ny; j++)
		//{
		//	u(0, j, 0) = u(Nx - 1, j, 0) = 100;
		//}
	}
	void solve() {
		for (int t = 1; t < tn; t++) {
			mat K1 = deriv_center(u.slice(t - 1));
			mat K2 = deriv_center(u.slice(t - 1) + 0.5 * dt * K1);
			mat K3 = deriv_center(u.slice(t - 1) + 0.5 * dt * K2);
			mat K4 = deriv_center(u.slice(t - 1) + dt * K3);
			u.slice(t) = u.slice(t - 1) + dt * (1.0 / 6 * K1 + 1.0 / 3 * K2 + 1.0 / 3 * K3 + 1.0 / 6 * K4);
		}
		u.slice(tn - 1).save("tn.csv", csv_ascii);
	}

private:

	int Nx, Ny, T;
	double Lx, Ly;
	double dt;
	double alpha;
	int tn;
	double dx, dy;
	cube u;
	vec x;
	vec t;
	vec error;

	mat deriv_center(const mat& f) {
		mat ddf(Nx, Ny);
		for (int i = 1; i < Nx - 1; i++) {
			for (int j = 1; j < Ny - 1; j++)
			{
				if (i == 1 || i == Nx - 2) {
					if (j == 1 || j == Ny - 2) {
						if (i == 1 && j == 1) ddf(i, j) = (f(i, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i + 2, j)) / (12 * dx * dx) + (f(i, j) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j + 2)) / (12 * dy * dy);
						if (i == 1 && j == Ny - 2) ddf(i, j) = (f(i, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i + 2, j)) / (12 * dx * dx) + (f(i, j - 2) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j)) / (12 * dy * dy);
						if (i == Nx - 2 && j == 1) ddf(i, j) = (f(i - 2, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i, j)) / (12 * dx * dx) + (f(i, j) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j + 2)) / (12 * dy * dy);
						if (i == Nx - 2 && j == Ny - 2) ddf(i, j) = (f(i - 2, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i, j)) / (12 * dx * dx) + (f(i, j - 2) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j)) / (12 * dy * dy);
						continue;
					}
					ddf(i, j) = 0;
					if (i == 1) ddf(i, j) = (f(i, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i + 2, j)) / (12 * dx * dx) + (f(i, j - 2) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j + 2)) / (12 * dy * dy);
					if (i == Nx - 2) ddf(i, j) = (f(i - 2, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i, j)) / (12 * dx * dx) + (f(i, j - 2) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j + 2)) / (12 * dy * dy);
					continue;
				}
				if (j == 1 || j == Ny - 2) {
					if (i == 1 || i == Ny - 2) {
						continue;
					}
					if (j == 1) ddf(i, j) = (f(i - 2, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i + 2, j)) / (12 * dx * dx) + (f(i, j) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j + 2)) / (12 * dy * dy);
					if (j == Ny - 2) ddf(i, j) = (f(i - 2, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i + 2, j)) / (12 * dx * dx) + (f(i, j - 2) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j)) / (12 * dy * dy);
					continue;
				}

				ddf(i, j) = (f(i - 2, j) - 16 * f(i - 1, j) + 30 * f(i, j) - 16 * f(i + 1, j) + f(i + 2, j)) / (12 * dx * dx) + (f(i, j - 2) - 16 * f(i, j - 1) + 30 * f(i, j) - 16 * f(i, j + 1) + f(i, j + 2)) / (12 * dy * dy);

			}
		}

		return alpha * ddf;
	}
};
int main() {
	auto simulate = HeatConductionEquationSolver();
	simulate.solve();
	system("pause");
	return 0;
}