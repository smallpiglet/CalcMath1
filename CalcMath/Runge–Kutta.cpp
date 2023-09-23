// CalcMath.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <cstring>
#include <omp.h>

//using namespace std;
//
//const int n = 2;


//double func(double* y, int i)
//{
//    double w;
//    switch (i)
//    {
//    case 0:
//        w = y[0] * (1.0 - sqrt(pow(y[0], 2) + pow(y[1], 2))) - y[1];
//        break;
//    case 1:
//        w = y[1] * (1.0 - sqrt(pow(y[0], 2) + pow(y[1], 2))) - y[0];
//        break;
//    }
//
//    return w;
//}
//
//int main()
//{
//    setlocale(LC_ALL, "Russian");
//    double t0 = 0.0, t = t0, tMax = 60, tau = 0.01, y[2] = { 1.0, 1.0 }, yy[2] = { 0.0 }, ff[n] = { 0.0 }, tStart, tEnd, deltaT;
//
//    tStart = omp_get_wtime();
//
//    do {
//        for (int i = 0; i < n; i++)
//            yy[i] = y[i] + tau * 0.5 *func(yy, i + 0.5 * tau);
//
//        for (int i = 0; i < n; i++)
//            ff[i] = y[i] + tau * 0.5 * func(yy, i + 0.5 * tau);
//
//        for (int i = 0; i < n; i++)
//            y[i] = ff[i];
//
//        t += tau;
//    } while (t <= tMax);
//
//    tEnd = omp_get_wtime();
//    deltaT = tEnd - tStart;
//    cout << "Время выполнения методом Рунгк-Кутта 2 = " << deltaT << endl << "Результаты:" << endl;
//
//    for (int i = 0; i < n; i++)
//        cout << "y[" << i << "] = " << y[i] << endl;
//
//    system("Pause");
//
//    return 0;
//}