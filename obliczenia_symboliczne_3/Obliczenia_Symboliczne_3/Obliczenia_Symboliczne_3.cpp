//Adam Majchrzak s176708 21.05.2020

#include <iostream>
#include <vector>

using namespace std;

const int n = 11;

//funkcja testowa
double f(double x)
{
	return (1.0 / (1.0 + x * x));
}

double f_1(double x)
{
	return sin(x);
}

//funkcja testowa, pokazuje dane w wektorze
void print_vector(std::vector<double> z) {
	for (int i = 0; i < z.size(); ++i) {
		std::cout << z[i] << " ";
	}
	std::cout << std::endl;
}

void enter_step(double zmienna, std::vector<double>& h_w) {
	h_w[0] = 0.0;
	for (int i = 1; i < h_w.size(); ++i) {
		std::cin >> zmienna;
		h_w[i] = zmienna;
	}
}

void enter_first_value(double zmienna, std::vector<double>& h_w, double xj[], double pj[]) {
	double jakas_zmienna = 0.0;
	std::cout << "\nwprowadz pierwsza wartosc zgodnie z ustalonym przez ciebie krokiem: ";
	std::cin >> zmienna;
	xj[0] = zmienna;
	std::cout << "\n";
	int wybor = 0;
	std::cout << "Do przetestowania programu przygotowalem 2 funkcje \n*(1) f(x) = 1.0 / (1.0 + x * x)\n"
		<< " druga z nich to \n*(2) - f(x) = sin(x)\nwybierz na ktorej przetestowac program 1/2 jesli wpiszesz "
		<< "inna liczbe calkowita program wykona sie dla f(x) = sin(x)\n";
	std::cin >> wybor;
	if (wybor == 1) {
		for (int j = 0; j < n; j++)
		{
			if (j > 0) xj[j] = xj[j - 1] + h_w[j];
			pj[j] = f(xj[j]);
			std::cout << xj[j] << " ";
		}
		std::cout << "\n";
	}
	else if (wybor == 2) {
		for (int j = 0; j < n; j++)
		{
			if (j > 0) xj[j] = xj[j - 1] + h_w[j];
			pj[j] = f_1(xj[j]);
			std::cout << xj[j] << " ";
		}
		std::cout << "\n";
	}
	else {
		for (int j = 0; j < n; j++)
		{
			if (j > 0) xj[j] = xj[j - 1] + h_w[j];
			pj[j] = f_1(xj[j]);
			std::cout << xj[j] << " ";
		}
		std::cout << "\n";
	}
}

void enter_data(vector<double>& d_w, vector<double>& a_w, vector<double>& b_w, vector<double>& c_w, vector<double>& h_w,
	double r[], double pj[], double a[], double b[], double c[]) {
	r[0] = 0.0;
	r[n - 1] = 0.0;
	d_w[0] = r[0];
	d_w[n - 1] = r[n - 1];
	for (int j = 1; j < (n - 1); j++)
	{
		r[j] = (6.0 / h_w[j]) * (pj[j + 1] - 2.0 * pj[j] + pj[j - 1]);
		d_w[j] = r[j];
	}
	print_vector(d_w);

	////Initial values for a,b,c
	a[0] = 0.0; //not used
	a[1] = 0.0;
	a[n - 1] = 0.0;
	a_w[0] = a[0];
	a_w[1] = a[1];
	a_w[n - 1] = a[n - 1];
	for (int j = 2; j < (n - 1); j++)
	{
		a[j] = h_w[j];
		a_w[j] = a[j];
	}
	print_vector(a_w);

	b[0] = 1.0;
	b[n - 1] = 1.0;
	b_w[0] = b[0];
	b_w[n - 1] = b[n - 1];
	for (int j = 1; j < (n - 1); j++)
	{
		b[j] = 4.0 * h_w[j];
		b_w[j] = b[j];
	}
	print_vector(b_w);

	c[0] = 0.0;
	c[n - 2] = 0.0;
	c[n - 1] = 0.0; //not used
	c_w[0] = c[0];
	c_w[n - 1] = c[n - 1];
	c_w[n - 2] = c[n - 2];
	for (int j = 1; j < (n - 2); j++)
	{
		c[j] = h_w[j];
		c_w[j] = c[j];
	}
	print_vector(c_w);
}

void enter_beta(vector<double>& beta, vector<double>& a_w, vector<double>& b_w, vector<double>& c_w) {
	double val = 0;
	val = -c_w[0] / b_w[0];
	beta.push_back(val);
	for (int i = 1; i < n; ++i) {
		val = -(c_w[i] / ((a_w[i] * beta[i - 1]) + b_w[i]));
		if (val == -0) val = 0; //warunek dla poprawnosci zapisu (nie może być -0)
		beta.push_back(val);
	}
	std::cout << "\nObliczona beta z naszych danych 1 do n: \n";
	print_vector(beta);
}

void enter_gamma(vector<double>& gamma,vector<double>& beta, vector<double>& a_w, vector<double>& b_w, vector<double>& d_w) {
	double val2 = 0;
	val2 = d_w[0] / b_w[0];
	gamma.push_back(val2);

	for (int i = 1; i < n; ++i) {
		val2 = (d_w[i] - a_w[i] * gamma[i - 1]) / (a_w[i] * beta[i - 1] + b_w[i]);
		gamma.push_back(val2);
	}
	std::cout << "\nObliczona gamma z naszych danych od 1 do n: \n";
	print_vector(gamma);
}

void score_of_this_program(vector<double>& gamma, vector<double>& x_w, vector<double>& beta) {
	double val1 = 0;
	x_w.push_back(gamma[n - 1]);
	int j = 0;
	for (int i = n; i > 1; --i) {
		val1 = (beta[i - 2] * x_w[j]) + (gamma[i - 2]);
		x_w.push_back(val1);
		j++;
	}
	std::cout << "\nNasze x w kolejnosci: pj,pj-1, ... pj3,pj2,pj1: \n";
	print_vector(x_w);
}
//
void wzor(std::vector<double>& x_w, double pj[], std::vector<double>& h_w, double xj[]) {
	std::vector<double> score_a(n);
	std::vector<double> score_b(n);
	std::vector<double> score_c(n);
	double zmienna_przechowuje_wynik_a = 0.0;
	for (int i = 1; i < 10; ++i) {
		std::cout << i << std::endl;
		zmienna_przechowuje_wynik_a = pj[i] + (	(pj[i + 1] - pj[i]) / h_w[i] - (1.0 / 6.0 * h_w[i] * x_w[i + 1]) - (1.0 / 3.0 * h_w[i] * x_w[i]));
		score_a[i] = zmienna_przechowuje_wynik_a;
	}		
	std::cout << "\n Dane do 'a'";
	print_vector(score_a);
	double zmienna_przechowuje_wynik_b = 0.0;
	for (int i = 1; i < 10; ++i) {
		zmienna_przechowuje_wynik_b = x_w[i]/2.0;
		score_b[i] = zmienna_przechowuje_wynik_b;
	}
	std::cout << "\n Dane do 'b'";
	print_vector(score_b);
	double zmienna_przechowuje_wynik_c = 0.0;
	for (int i = 1; i < 10; ++i) {
		zmienna_przechowuje_wynik_c = (x_w[i+1] - x_w[i]) / (6.0*h_w[i]);
		score_c[i] = zmienna_przechowuje_wynik_c;
	}
	std::cout << "\n Dane do 'c'";
	print_vector(score_c);
	std::cout << std::endl;
	for (int i = 1; i < 12; i++) {
		std::cout << i << ". p(x) = " << score_a[i - 1] << "(x + (" << xj[i - 1] << ")) + " << score_b[i - 1] << "(x + (" << xj[i - 1] << "))^2 + "
			<< score_c[i - 1] << "(x + (" << xj[i - 1] << "))^3" << std::endl;
	}
}



int main() {
	std::vector<double> a_w(n); //wektor przechowuje wartosci nad diagonala
	std::vector<double> b_w(n); //wektor przechowuje wartosci diagonali
	std::vector<double> c_w(n); //wektor przechowuje wartosci pod diagonala
	std::vector<double> d_w(n); //wektor przechowuje wartosci wynikowe
	std::vector<double> beta;
	std::vector<double> gamma;
	std::vector<double> x_w;
	std::vector<double> h_w(n); //wektor przechowujje wartości kroku !
	double xj[n], a[n], b[n], c[n], r[n], pjPP[n], pj[n];
	double h, p;
	
	h = 1.0;
	std::cout << "1. Dla zestawu wezlow x={-6, -5, -3, -2, 0, 1.5, 3, 4, 6, 7, 8}\n"
		//		<< "Testuje dla zestawu: \nZatem wprowadzam odgornie "
		//		<< "narzucone 'h' czyli krok, ktory prezentuje sie nastepujaco:\n"
		<< "krok h wynosi: \n1, 2, 1, 2, 1.5, 1.5, 1, 2, 1, 1\n"
		<< "\n2 . Drugi zestaw dla ktorego bedziemy testowac program to x={-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5}\n"
		<< "krok h wynosi: \n1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1"
		<< "\n\n\nProsze wybrac zestaw testowy. a nastepnie wpisac jego krok z klawiatury\n";

	double jakas_zmienna = 0.0;
	enter_step(jakas_zmienna, h_w);
	enter_first_value(jakas_zmienna, h_w, xj, pj);
	enter_data(d_w, a_w, b_w, c_w, h_w, r,  pj,  a,  b,  c);
	enter_beta(beta, a_w, b_w, c_w);
	enter_gamma(gamma, beta, a_w, b_w, d_w);
	score_of_this_program(gamma, x_w, beta);
	wzor(x_w, pj, h_w, xj);

	for (int i = 0; i < 10; i++)std::cout << pj[i] << " ";
	return 0;
}