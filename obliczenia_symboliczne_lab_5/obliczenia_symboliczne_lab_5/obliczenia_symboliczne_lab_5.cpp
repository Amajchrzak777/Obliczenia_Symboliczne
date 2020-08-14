// Adam Majchrzak s176708 02.06.2020

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>


//deklaracje funkcji
double f(double x);
double gnti(double& x, double t);
double mf(double x, double t);
double tntm(double& x, double t);
double calcUs(double x, double tm);
double funcU(double x, double t);
double metodaProstokatow();
void vector_ti();

std::vector<double> t1_ti;
double j = 1.0;
static const int ilosc_sum = 15;
double tau = 1.0;

//double f(double x) {
//    return std::exp(-std::pow(x - 5, 2));
//}
//
double calcLambda(double k) {
    return (k * M_PI) / 10;
}

double Cm(double x, double k) {
    return f(x) * std::sin(sqrt(calcLambda(k)) * x);
}

double tntm(double& x, double t) {
    double Tn = 0.0;
    for (int i = 0; i < ilosc_sum; i++) {
        Tn += std::exp(pow(calcLambda(i),2) * (t1_ti[i] - t1_ti[ilosc_sum])) * gnti(x, t1_ti[i]) * tau;
    }
    return Tn;
}
//
double calcUs(double x, double tm) {
    double us = 0.0;
    for (int i = 0; i < ilosc_sum; i++) {
        us += tntm(x, tm) * sin(calcLambda(i) * x);
    }
    return us;
}

double metodaProstokatow() {
    const int N = 100; //liczba punktów/prostok¹tów podzia³owych
    double xp, xk, s, dx;
    xp = 0;
    xk = 2*M_PI;
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
        u += metodaProstokatow() * std::exp(pow(calcLambda(i),2) * t) * std::sin(calcLambda(i) * x);
    }
    //std::cout << u << std::endl;
    return u;
}

double f(double x) {
    return sin(x);
}

double gnti(double& x, double t) {
    //deklaruje tablice, obliczam hx, xm, a nastêpnie sumuje g(xm, ti)d
    //dla czasu t = 1 zmienna do zmieniania dla ró¿nych czasów
    double xm[1000];
    xm[0] = x;
    double hx = ilosc_sum/1000.0;
    //std::cout << "hx: "<<hx << std::endl ;
    double L = 2 * M_PI;
    //std::cout << "l/hx"  <<L / hx << "\n";
    double nx = L / hx;
    for (int i = 1; i < nx; i++) {
        xm[i] = xm[0] + (i * hx);
    }
    //std::cout << "nx = " << nx << std::endl;
    double sum = 0.0;
    for (int i = 0; i < nx; i++) {
        sum += mf(xm[i], t)* sin(calcLambda(i)*xm[i]);
    }
    return (2/nx)*sum;
}

double mf(double x, double t) {
    return 2 * std::exp(-t) * std::sin(x);
}

void vector_ti() {
    double zmienna = 0.0;
    for (int i = 0; i <= 420; ++i) {
        zmienna = tau * i;
        t1_ti.push_back(zmienna);
    }

}

int main()
{
   
    double j = 1;
    
    
    vector_ti();
    //for (auto n : t1_ti) std::cout << n << std::endl;
    

   
    std::cout << "Warunki brzegowe dla Us1: " << std::endl;
    std::cout << "Us1: x = 0, t = 0: " << calcUs(0, 0) << std::endl;
    std::cout << "Us1: x = 2*pi, t = 0: " << calcUs(2*M_PI, 0) << std:: endl;

    std::cout << "Warunki brzegowe dla Us1: " << std::endl;
    std::cout << "Us2: x = 0, t = 0: " << funcU(0, 0) << std::endl;
    std::cout << "Us2: x = 2*pi, t = 0: " << funcU(2 * M_PI, 0) << std::endl;
    double x = funcU(2 * M_PI, 0);
    double y = calcUs(2 * M_PI, 0);
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << x + y << std::endl;

    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = 1, t = 0: " << funcU(1, 0) << std::endl;
    std::cout << "Us2: x = 1, t = 0: " << calcUs(1, 0) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(1, 0) + calcUs(1, 0);

    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = 2, t = 0: " << funcU(2, 0) << std::endl;
    std::cout << "Us2: x = 2, t = 0: " << calcUs(2, 0) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(2, 0) + calcUs(2, 0);


    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = pi, t = 5: " << funcU(M_PI, 5) << std::endl;
    std::cout << "Us2: x = pi, t = 5: " << calcUs(M_PI, 5) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(M_PI, 5) + calcUs(M_PI, 5);




    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = 1, t = 5: " << funcU(1, 5) << std::endl;
    std::cout << "Us2: x = 1, t = 5: " << calcUs(1, 5) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(1, 5) + calcUs(1, 5);

    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = 2, t = 5: " << funcU(2, 5) << std::endl;
    std::cout << "Us2: x = 2, t = 5: " << calcUs(2, 5) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(2, 5) + calcUs(2, 5);

    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = 3, t = 5: " << funcU(3, 5) << std::endl;
    std::cout << "Us2: x = 3, t = 5: " << calcUs(3, 5) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(3, 5) + calcUs(3, 5);

    std::cout << "\nTEST: " << std::endl;
    std::cout << "Us2: x = 4, t = 5: " << funcU(4, 5) << std::endl;
    std::cout << "Us2: x = 4, t = 5: " << calcUs(4, 5) << std::endl;
    std::cout << "Us = Us1 + Us2 = ";
    std::cout << funcU(4, 5) + calcUs(4, 5);

}

