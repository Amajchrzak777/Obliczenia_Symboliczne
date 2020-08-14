// Adam Majchrzak s176708 02.06.2020

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

double j = 1;
static const int ilosc_sum = 100;

double f(double x) {
    return std::exp(-std::pow(x - 5, 2));
}

double calcLambda(double k) {
    return std::pow((k * M_PI) / 10, 2);
}

double Cm(double x, double k) {
    return f(x) * std::sin(sqrt(calcLambda(k)) * x);
}

double metodaProstokatow() {
    const int N = 100; //liczba punktów/prostokątów podziałowych
    double xp, xk, s, dx;
    xp = 0;
    xk = 10;
    s = 0;
    int i;
    dx = (xk - xp) / N;
    for (i = 1; i <= N; i++) s += Cm(xp + i * dx, j);
    s *= dx;
    j++;
    if (j == ilosc_sum) {
        j = 1;
    }

    double score1 = (2 / xk) * s;
    return score1;
}

double funcU(double x, double t)
{
    double u = 0.0;
    for (int i = 1; i < ilosc_sum; i++) {
        u += metodaProstokatow() * std::exp(-calcLambda(i) * t) * std::sin(std::sqrt(calcLambda(i)) * x);
    }
    //std::cout << u << std::endl;
    return u;
}

int main()
{
    double j = 1;
    std::cout << "Rownanie przewodnictwa cieplnego w precie ograniczonym jednorodnymi warunkami brzegowymi. \n";
    std::cout << "--------------------------------------------------------------------------------------------";
    std::cout << "\nWarunki brzegowe: \n";
    std::cout << "u(t,10) = 0\n";
    std::cout << "u(t,0) = 0\n";
    std::cout << "Nasze rownanie: \n";
    std::cout << "u(t,x) = X(x)*T(t)\n";
    std::cout << "Warunek poczatkowy: \nu(0,x) = f(x) \t dla 0 < x < l";
    std::cout << "f(x) = exp(-(x-5)^2)";

    std::cout << "Na początku rozdzielam zmienne: \n";
    std::cout << "d^2X/dx^2 + lambda*X(x) = 0\n";
    std::cout << "dT/dt + lambda*T(t) = 0\n";

    std::cout << "T(t) = Ce^(-lambda*t)\n";
    std::cout << "X(x) = A*cos(sqrt(lambda)*x) + B*sin(sqrt(lambda)*x)";
    std::cout << "u(t,x) = Ce^(-lambda*t)*A*cos(sqrt(lambda)*x) + B*sin(sqrt(lambda)*x)";
    std::cout << "\nKorzystajac z warunkow brzegowych obliczam lambda \n";
    std::cout << "Wstawiam do wzoru na Cm, obliczone wstawiam do Sumy i otrzymuje wyniki: \n";
    std::cout << "Przyjete granice calkowania 0 - 10, ilosc sum = 100, \nliczba prostokatow do przyblizenia na przedziale - 100";
    std::cout << "\nPrzyjeta precyzja do 15 liczb po przecinku.";
    std::cout << std::endl << std::fixed << std::setprecision(15);
    std::cout << "t = 1\n" <<funcU(0, 1) << '\n';
    std::cout << funcU(1, 1) << '\n';
    std::cout << funcU(2, 1) << '\n';
    std::cout << funcU(3, 1) << '\n';
    std::cout << funcU(4, 1) << '\n';
    std::cout << funcU(5, 1) << '\n';
    std::cout << funcU(6, 1) << '\n';
    std::cout << funcU(7, 1) << '\n';
    std::cout << funcU(8, 1) << '\n';
    std::cout << funcU(9, 1) << '\n';
    std::cout << funcU(10,1) << '\n';

    std::cout << std::endl << std::endl;
    std::cout << "t = 5\n" << funcU(0, 5) << '\n';
    std::cout << funcU(1, 5) << '\n';
    std::cout << funcU(2, 5) << '\n';
    std::cout << funcU(3, 5) << '\n';
    std::cout << funcU(4, 5) << '\n';
    std::cout << funcU(5, 5) << '\n';
    std::cout << funcU(6, 5) << '\n';
    std::cout << funcU(7, 5) << '\n';
    std::cout << funcU(8, 5) << '\n';
    std::cout << funcU(9, 5) << '\n';
    std::cout << funcU(10, 5) << '\n';

    std::cout << std::endl << std::endl;
    std::cout << "t = 10\n" << funcU(0, 10) << '\n';
    std::cout << funcU(1, 10) << '\n';
    std::cout << funcU(2, 10) << '\n';
    std::cout << funcU(3, 10) << '\n';
    std::cout << funcU(4, 10) << '\n';
    std::cout << funcU(5, 10) << '\n';
    std::cout << funcU(6, 10) << '\n';
    std::cout << funcU(7, 10) << '\n';
    std::cout << funcU(8, 10) << '\n';
    std::cout << funcU(9, 10) << '\n';
    std::cout << funcU(10, 10) << '\n';

    std::cout << std::endl << std::endl;
    std::cout << "t = 15\n" << funcU(0, 15) << '\n';
    std::cout << funcU(1, 15) << '\n';
    std::cout << funcU(2, 15) << '\n';
    std::cout << funcU(3, 15) << '\n';
    std::cout << funcU(4, 15) << '\n';
    std::cout << funcU(5, 15) << '\n';
    std::cout << funcU(6, 15) << '\n';
    std::cout << funcU(7, 15) << '\n';
    std::cout << funcU(8, 15) << '\n';
    std::cout << funcU(9, 15) << '\n';
    std::cout << funcU(10, 15) << '\n';
}

