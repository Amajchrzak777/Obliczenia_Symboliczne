// Adam Majchrzak s176708 05.05.2020
//
//
//#include <iostream>
//#include <vector>
//#include <fstream>
//
//void fill_vector(std::vector<double>& z, int n) {
//    double value;
//
//    for (int i = 0; i < n; ++i) {
//        std::cin >> value;
//        z.push_back(value);
//    }
//}
//
//void print_vector(std::vector<double> z) {
//    for (int i = 0; i < z.size(); ++i) {
//        std::cout << z[i] << " ";
//    }
//    std::cout << std::endl;
//}
//
//
//int main()
//{
//    std::cout << "Mamy macierz trojdiagonalna\n";
//    int n;
//    std::cout << "Wpisz ile wymiarow ma miec macierz (rozwazamy macierze kwadratowe wiec wystarczy 1 cyfra))?\n";
//    std::cin >> n;
//    std::vector<double> a;
//    std::vector<double> b;
//    std::vector<double> c;
//    std::vector<double> d;
//    std::vector<double> beta;
//    std::vector<double> gamma;
//    std::vector<double> x;
//    //uzupelnij vectory
//    std::cout << "Wprowadz dane do przekatnej a ktora znajduje sie pod diagonala.\n";
//    fill_vector(a,n);
//    std::cout << "Wprowadz dane do przekatnej b ktora znajduje sie nad diagonala.\n";
//    fill_vector(b,n);
//    std::cout << "Wprowadz dane do przekatnej c ktora jest diagonala.\n";
//    fill_vector(c,n);
//    std::cout << "Wprowadz dane do macierzy wynikowej d\n";
//    fill_vector(d,n);
//    //pokaz vectory
//    a[0] = 0;
//    c[n-1] = 0;
//    std::cout << "Nasze wektory: \n";
//    std::cout << "a: ";
//    print_vector(a);
//    std::cout << "\nb: ";
//    print_vector(b);
//    std::cout << "\nc: ";
//    print_vector(c);
//    std::cout << "\nd: ";
//    print_vector(d);
//    double val = 0;
//    //double val1 = 0;
//    val = -c[0] / b[0];
//    beta.push_back(val);
//    for (int i = 1; i < n; ++i) {
//        val = -(c[i] / ((a[i] * beta[i-1]) + b[i]));
//        if (val == -0) val = 0; //warunek dla poprawnosci zapisu (nie może być -0)
//        beta.push_back(val);
//    }
//    std::cout << "\nObliczona beta z naszych danych 1 do n: ";
//    print_vector(beta);
//    
//    double val2 = 0;
//    val2 = d[0] / b[0];
//    gamma.push_back(val2);
//    for (int i = 1; i < n; ++i) {
//        val2 = (d[i] - a[i] * gamma[i - 1]) / (a[i] * beta[i - 1] + b[i]);
//        gamma.push_back(val2);
//    }
//    std::cout << "\nObliczona gamma z naszych danych od 1 do n: ";
//    print_vector(gamma);
//    
//
//    double val1 = 0;
//    x.push_back(gamma[n-1]);
//    int j = 0;
//    for (int i = n; i > 1; --i) {
//        val1 = (beta[i - 2] * x[j]) + (gamma[i - 2]);
//        x.push_back(val1);
//        j++;
//    }
//    std::cout << "\nNasze x w kolejnosci: xn,xn-1, ... x3,x2,x1: ";
//    print_vector(x);
//}


#include <iostream>
#include <vector>

using namespace std;

//funkcja testowa
double f(double x)
{
	return (1.0 / (1.0 + x * x));
}

//funkcja testowa, pokazuje dane w wektorze
void print_vector(std::vector<double> z) {
	for (int i = 0; i < z.size(); ++i) {
		std::cout << z[i] << " ";
	}
	std::cout << std::endl;
}

int main() {
	std::vector<double> a_w(11); //wektor przechowuje wartosci nad diagonala
	std::vector<double> b_w(11); //wektor przechowuje wartosci diagonali
	std::vector<double> c_w(11); //wektor przechowuje wartosci pod diagonala
	std::vector<double> d_w(11); //wektor przechowuje wartosci wynikowe
	std::vector<double> beta;
	std::vector<double> gamma;
	std::vector<double> x_w;

	const int n = 11;
	double xj[n], a[n], b[n], c[n], r[n], pjPP[n], pj[n], Beta[n], Rho[n];
	double x, h, p, xR;

	//Initial values	
	h = 1.0;

	//Initial values for pj and xj
	xj[0] = -5.0;
	for (int j = 0; j < n; j++)
	{
		if (j > 0) xj[j] = xj[j - 1] + h;
		pj[j] = f(xj[j]);
	}

	//Initial values for r	
	r[0] = 0.0;
	r[n - 1] = 0.0;
	d_w[0] = r[0];
	d_w[n - 1] = r[n - 1];
	for (int j = 1; j < (n - 1); j++)
	{
		r[j] = (6.0 / h) * (pj[j + 1] - 2.0 * pj[j] + pj[j - 1]);
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
		a[j] = h;
		//std::cout << a[j] << '\n';
		a_w[j] = a[j];
	}
	print_vector(a_w);

	b[0] = 1.0;
	b[n - 1] = 1.0;
	b_w[0] = b[0];
	b_w[n - 1] = b[n - 1];
	for (int j = 1; j < (n - 1); j++)
	{
		b[j] = 4.0 * h;
		//std::cout << b[j] << '\n';
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
		c[j] = h;
		//std::cout << c[j] << '\n';
		c_w[j] = c[j];
	}
	print_vector(c_w);

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

	double val2 = 0;
	val2 = d_w[0] / b_w[0];
	//std::cout << val2;
	gamma.push_back(val2);

	for (int i = 1; i < n; ++i) {
		val2 = (d_w[i] - a_w[i] * gamma[i - 1]) / (a_w[i] * beta[i - 1] + b_w[i]);
		gamma.push_back(val2);
	}
	std::cout << "\nObliczona gamma z naszych danych od 1 do n: \n";
	print_vector(gamma);


	double val1 = 0;
	x_w.push_back(gamma[n - 1]);
	int j = 0;
	for (int i = n; i > 1; --i) {
		val1 = (beta[i - 2] * x_w[j]) + (gamma[i - 2]);
		x_w.push_back(val1);
		j++;
	}
	std::cout << "\nNasze x w kolejnosci: xn,xn-1, ... x3,x2,x1: ";
	print_vector(x_w);


	std::cout << "\n------------------------------------------------------\n";
	//Gaussian Elimination
	Beta[0] = b[0];
	Rho[0] = r[0];
	for (int j = 1; j < n; j++)
	{
		Beta[j] = b[j] - (a[j] / Beta[j - 1]) * c[j - 1];
		Rho[j] = r[j] - (a[j] / Beta[j - 1]) * Rho[j - 1];
	}
	pjPP[n - 1] = Rho[n - 1] / Beta[n - 1];
	for (int j = 2; j <= n; j++)
	{
		pjPP[n - j] = (Rho[n - j] - c[n - j] * pjPP[n - j + 1]) / Beta[n - j];
	}

	//Print pjPP
	for (int j = 0; j < n; j++)
	{
		cout << "p" << j + 1 << "PP= " << pjPP[j] << endl;
	}

	////Calculation of the polynomial p
	//x = -5.0;
	//xR = xj[1]; //initial interval xj[0] xj[1]
	//int j = 0;
	//do
	//{
	//	if (x > xR) //jump to next interval
	//	{
	//		j++;
	//		xR = xj[j + 1];
	//	}
	//	p = pj[j] + ((pj[j + 1] - pj[j]) / h - (h * pjPP[j + 1]) / 6.0 - (h * pjPP[j]) / 3.0) * (x - xj[j]);
	//	p = p + (pjPP[j] / 2.0) * (x - xj[j]) * (x - xj[j]) + ((pjPP[j + 1] - pjPP[j]) / (6.0 * h)) * (x - xj[j]) * (x - xj[j]) * (x - xj[j]);

	//	cout << x << "  " << p << "  " << f(x) << endl;
	//	x = x + 0.1;
	//} while (x < 5.0001);

	return 0;
}