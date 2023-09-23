// CalcMath.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <cstring>
#include <omp.h>

using namespace std;

const int n = 2;
const int m = 4;

void PrintMethodResults(double* y, double deltaT, string methodName)
{
	cout << "-=" << methodName << "=-" << endl << " Время выполнения = " << deltaT << endl << " Результаты:" << endl;

	for (int i = 0; i < n; i++)
		cout << '\t' << "y[" << i << "] = " << y[i] << endl;

	cout << endl;
}

double ApproximateDeravitive(double (*func)(double*, double), double* y, double t, short dArg)
{
	const double delta = 0.001;
	double* yy;

	if (dArg == 1)
		yy = new double[n] { y[0] + delta, y[1]};
	else
		yy = new double[n] { y[0], y[1] + delta};

	double res = (func(yy, t) - func(y, t)) / delta;

	delete[] yy;
	return res;
}

double F1(double *y, double time) {
	return y[0] * (1.0 - sqrt(pow(y[0], 2) + pow(y[1], 2))) - y[1];
}

double F2(double* y, double time) {
	return y[1] * (1.0 - sqrt(pow(y[0], 2) + pow(y[1], 2))) - y[0];
}

void ExplicitEuler(double* y, int size, double (*funcs[])(double*, double), double t, double tMax, double tau)
{
	double *yy = new double[size] { 0.0 };

	do {
		for (int i = 0; i < size; i++)
			yy[i] = y[i] + tau * (funcs[i](y, t));

		for (int i = 0; i < size; i++)
			y[i] = yy[i];

		t += tau;
	} while (t <= tMax);

	delete[] yy;
}

void RungeKutta2(double* y, int size, double (*funcs[])(double*, double), double t, double tMax, double tau)
{
	double ff[n] = { 0.0 };
	double *yy = new double[size] { 0.0 };

	do {
		for (int i = 0; i < size; i++)
			yy[i] = y[i] + tau * 0.5 * funcs[i](y, t);

		for (int i = 0; i < size; i++)
			ff[i] = tau * funcs[i](yy, t + 0.5 * tau);

		for (int i = 0; i < size; i++)
			y[i] += ff[i];

		t += tau;
	} while (t <= tMax);

	delete[] yy;
}

void PredictorCorrector(double* y, int size, double (*funcs[])(double*, double), double t, double tMax, double tau)
{
	double ff[n] = { 0.0 };
	double* yy = new double[size] { 0.0 };

	do {
		for (int i = 0; i < size; i++)
			yy[i] = y[i] + tau * funcs[i](y, t);

		for (int i = 0; i < size; i++)
			ff[i] = tau * (funcs[i](y, t) + funcs[i](yy, t + tau)) / 2.0;

		for (int i = 0; i < size; i++)
			y[i] += ff[i];

		t += tau;
	} while (t <= tMax);

	delete[] yy;
}

void RungeKutta4(double* y, int size, double (*funcs[])(double*, double), double t, double tMax, double tau)
{
	double *yy = new double[size] { 0.0 };
	double *R = new double[m, size]{ 0.0 };

	do {
		for (int i = 0; i < size; i++)
			R[0, i] = tau * funcs[i](y, t);

		for (int i = 0; i < size; i++)
			yy[i] = y[i] + 0.5 * R[0, i];

		for (int i = 0; i < size; i++)
			R[1, i] = tau * funcs[i](y, t + 0.5 * tau);

		for (int i = 0; i < size; i++)
			yy[i] = y[i] + 0.5 * R[1, i];

		for (int i = 0; i < size; i++)
			R[2, i] = tau * funcs[i](yy, t + 0.5 * tau);

		for (int i = 0; i < size; i++)
			yy[i] = y[i] + R[2, i];

		for (int i = 0; i < size; i++)
			R[3, i] = tau * funcs[i](yy, t + tau);

		for (int i = 0; i < size; i++)
			y[i] += (R[0 , i] + 2 * R[1, i] + 2 * R[2, i] + R[3, i]) / 6.0;

		t += tau;
	} while (t <= tMax);

	delete[] yy, R;
}

void ImplicitEuler(double* y, int size, double (*funcs[])(double*, double), double t, double tMax, double tau)
{
	double bi[n] = { 0.0 };
	double p[n] = { 0.0 };

	do {
		for (int i = 0; i < size; i++)
			bi[i] = -funcs[i](y, t);

		double R[2][2] = { { ApproximateDeravitive(funcs[0], y, t, 1) - 1.0 / tau, ApproximateDeravitive(funcs[0], y, t, 2)},
		   { ApproximateDeravitive(funcs[1], y, t, 1), ApproximateDeravitive(funcs[1], y, t, 2) - 1.0 / tau} };

		double det = R[0][0] * R[1][1] - R[0][1] * R[1][0];

		p[0] = R[1][1] * bi[0] - R[0][1] * bi[1];
		p[1] = R[0][0] * bi[1] - R[1][0] * bi[0];

		for (int i = 0; i < size; i++)
			y[i] += p[i] / det;

		t += tau;
	} while (t <= tMax);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	double t0 = 0.0, t = t0, tMax = 100, tau = 0.0001, y[n] = { 1.0, 1.0 }, tStart, tEnd, deltaT;

	double (*funcs[])(double*, double) = { F1, F2 };

	tStart = omp_get_wtime();
	ExplicitEuler(y, n, funcs, t, tMax, tau);
	tEnd = omp_get_wtime();

	deltaT = tEnd - tStart;
	
	PrintMethodResults(y, deltaT, "Метод Эйлера (Явный)");

	t = t0;
	y[0] = 1.0; y[1] = 1.0;

	tStart = omp_get_wtime();

	RungeKutta2(y, n, funcs, t, tMax, tau);

	tEnd = omp_get_wtime();

	deltaT = tEnd - tStart;

	PrintMethodResults(y, deltaT, "Метод Рунге-Кутта 2");

	t = t0;
	y[0] = 1.0; y[1] = 1.0;

	tStart = omp_get_wtime();

	PredictorCorrector(y, n, funcs, t, tMax, tau);

	tEnd = omp_get_wtime();

	deltaT = tEnd - tStart;

	PrintMethodResults(y, deltaT, "Метод Предиктор-Корретор");

	t = t0;
	y[0] = 1.0; y[1] = 1.0;

	tStart = omp_get_wtime();

	RungeKutta4(y, n, funcs, t, tMax, tau);

	tEnd = omp_get_wtime();

	deltaT = tEnd - tStart;

	PrintMethodResults(y, deltaT, "Метод Рунге-Кутта 4");

	t = t0;
	y[0] = 1.0; y[1] = 1.0;

	tStart = omp_get_wtime();

	ImplicitEuler(y, n, funcs, t, tMax, tau);

	tEnd = omp_get_wtime();

	deltaT = tEnd - tStart;

	PrintMethodResults(y, deltaT, "Метод Эйлера (Неявный)");

	system("Pause");

	cout << "ТЕСТ МАКС 1 залив" << endl;

	cout << "Я ненавижу вычислительную математику" << endl;

	return 0;
}