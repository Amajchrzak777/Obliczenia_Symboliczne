// Adam Majchrzak s176708 16.04.2020 
//do zadania uzylem wielomiany z maksymalna iloscia pierwsiatkow rzeczywistych
//w funkcji main() szukam wywoluje szukanie kazdego pierwiastka z osobna dla wielomianu
// 3 i 5 stopnia w momencie kiedy pierwiastek potrzbuje więcej iteracji szukam az do znalezienia,

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
void print_wielomian_5() {
    std::cout << "x ^ 5 - 5 * x ^ 4 - 70 * x ^ 3 + 230 * x ^ 2 + 789 * x - 945" << std::endl;
}
void print_wielomian_5p() {
    std::cout << "Pochodnia wielomianu:\n5 * x ^ 4 - 20 * x ^ 3 - 210 * x ^ 2 + 460 * x + 789" << std::endl;
}
void print_wielomian() {
    std::cout << "x ^ 3 + 2 * x ^ 2 - 5 * x - 6" << std::endl;
}
void pochodna_print_wielomian() {
    std::cout << "Pochodna wielomianu:\n3 * x ^ 2 + 4 * x - 5" << std::endl;
}
//Wielomian 5 stopnia
double f5(double x) {
    return (pow(x, 5) - 5 * pow(x, 4) - 70 * pow(x, 3) + 230 * pow(x, 2) + 789 * x - 945);
}
//Pochodna wielomianiu
double f5_p(double x) {
    return (5 * pow(x, 4) - 20 * pow(x, 3) - 210 * pow(x, 2) + 460 * x + 789);
}
//Wielomian 3 stopnia
double f(double x) {
    return (pow(x,3) + 2*pow(x,2)-5*x-6);
}
//Pochodna wielomianu
double f_p(double x) {
    return (3*pow(x,2) + 4*x - 5);
}


//Metoda Newtona
void Newton_Method(double x, int n) {
    double x2;
    for (int i = 0; i < n; ++i) {
        x2 = x - (f(x) / f_p(x));
        std::cout << "i= " << i << " x0= " << std::setprecision(5) << x2 << std::endl;
        x = x2;
    }
    std::cout << std::endl;
}

//Metoda siecznych
void Metoda_siecznych(double x1, double x2, int n) {
    double X;
    for (int i = 0; i < n; ++i) {
        X = x1 - f(x1) * (x1 - x2) / (f(x1) - f(x2));
        std::cout << "i= " << i << " x0= " << std::setprecision(5) << X << std::endl;
        x1 = X;
    }
    std::cout << std::endl;
}

//Metoda bisekcji
void Bisection_method(double x1, double x2, int n) {
    double x;
    if (f(x1) * f(x2) > 0.0) std::cout << "Wrong scope";
    for (int i = 0; i < n; ++i) {
        x = (x1 + x2) / 2;
        if ((f(x1) * f(x)) > 0.0) {
            x1 = x;
        }
        else {
            x2 = x;
        }
        std::cout << "i= " << i << " x0= " << std::setprecision(5) << x << std::endl;
    }
    std::cout << std::endl;
}

//Metoda Newtona
void Newton_Method5(double x, int n) {
    double x2;
    for (int i = 0; i < n; ++i) {
        x2 = x - (f5(x) / f5_p(x));
        std::cout << "i= " << i << " x0= " << std::setprecision(5) << x2 << std::endl;
        x = x2;
    }
    std::cout << std::endl;
}

//Metoda siecznych
void Metoda_siecznych5(double x1, double x2, int n) {
    double X;
    for (int i = 0; i < n; ++i) {
        X = x1 - f5(x1) * (x1 - x2) / (f5(x1) - f5(x2));
        std::cout << "i= " << i << " x0= " << std::setprecision(5) << X << std::endl;
        x1 = X;
    }
    std::cout << std::endl;
}

//Metoda bisekcji
void Bisection_method5(double x1, double x2, int n) {
    double x;
    if (f5(x1) * f5(x2) > 0.0) std::cout << "Wrong scope";
    for (int i = 0; i < n; ++i) {
        x = (x1 + x2) / 2;
        if ((f5(x1) * f5(x)) > 0.0) {
            x1 = x;
        }
        else {
            x2 = x;
        }
        std::cout << "i= " << i << " x0= " << std::setprecision(5) << x << std::endl;
    }
    std::cout << std::endl;
}
int main()
{

    //Dla wielomianów 3 stopnia:
    //Newton method
    std::cout << "Metoda Newtona-Raphsona: " << std::endl;
    std::cout << "Wielomian 3 stopnia, ktorego uzywam do zwizualizowania rozwiazanego zadania: " << std::endl;
    std::cout << "Wielomian posiada 3 pierwastki. " << std::endl;
    print_wielomian();
    pochodna_print_wielomian();
    Newton_Method(2.75, 5);
    Newton_Method(-0.2, 5);
    Newton_Method(-5, 5);
    //Metoda siecznych
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "Metoda Siecznych: " << std::endl;
    std::cout << "Wielomian 3 stopnia, ktorego uzywam do zwizualizowania rozwiazanego zadania: " << std::endl;
    std::cout << "Wielomian posiada 3 pierwastki. " << std::endl;
    print_wielomian();
    pochodna_print_wielomian();
    Metoda_siecznych(1, 5, 28);
    Metoda_siecznych(-2.9, 0.9, 15);
    Metoda_siecznych(-10, -5, 25);
    //Metoda bisekcji
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "Metoda Bisekcji: " << std::endl;
    std::cout << "Wielomian 3 stopnia, ktorego uzywam do zwizualizowania rozwiazanego zadania: " << std::endl;
    std::cout << "Wielomian posiada 3 pierwastki. " << std::endl;
    print_wielomian();
    pochodna_print_wielomian();
    Bisection_method(1, 5, 25);
    Bisection_method(-5, -2.5, 25);
    Bisection_method(-1.5, 0, 25);
    std::cout << "Jak widac wielomian w postaci iloczynowej wyraza sie wzorem (x+1)(x-2)(x+3)";
    std::cout << std::endl<< "-----------------------------------------------------------------------------" << std::endl << std::endl;
    //Dla wielomianów 5 stopnia:
    //Newton method
    std::cout << "Metoda Newtona-Raphsona: " << std::endl;
    std::cout << "Wielomian 5 stopnia, ktorego uzywam do zwizualizowania rozwiazanego zadania: " << std::endl;
    std::cout << "Wielomian posiada 5 pierwastkow. " << std::endl;
    print_wielomian_5();
    print_wielomian_5p();
    Newton_Method5(13,10);
    Newton_Method5(-11,10);
    Newton_Method5(7,10);
    Newton_Method5(-5,10);
    Newton_Method5(0,10);
    //Metoda siecznych
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "Metoda Siecznych: " << std::endl;
    std::cout << "Wielomian 5 stopnia, ktorego uzywam do zwizualizowania rozwiazanego zadania: " << std::endl;
    std::cout << "Wielomian posiada 5 pierwastkow. " << std::endl;
    print_wielomian_5();
    print_wielomian_5p();
    Metoda_siecznych5(10, 8, 45);
    Metoda_siecznych5(-9, -6, 45);
    Metoda_siecznych5(6, 4, 25);
    Metoda_siecznych5(-5, -2, 25);
    Metoda_siecznych5(-2, 2.5, 25);
    //Metoda bisekcji
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "Metoda Bisekcji: " << std::endl;
    std::cout << "Wielomian 5 stopnia, ktorego uzywam do zwizualizowania rozwiazanego zadania: " << std::endl;
    std::cout << "Wielomian posiada 5 pierwastkow. " << std::endl;
    print_wielomian_5();
    print_wielomian_5p();
    Bisection_method5(10, 8, 28);
    Bisection_method5(-9, -6, 25);
    Bisection_method5(6, 3, 25);
    Bisection_method5(-5, -2, 25);
    Bisection_method5(-0.5, 2.5, 25);
    std::cout << "Jak widac wielomian w postaci iloczynowej wyraza sie wzorem (x-9)(x+7)(x-5)(x+3)(x-1)";
    std::cout << std::endl << "-----------------------------------------------------------------------------" << std::endl << std::endl;
}

