

//https://www.bragitoff.com/2017/12/legendre-polynomial-c-program/

namespace ps
{


	template <typename T>
	T L0
	(
		const T& x
	)
	{
		return 1;
	}

	//Roots[LegendreP[2, cos(x)] == 0, x]
	//x = π n - 1/2 cos^(-1)(-1/3) and n element Z
	//54.7356°

	template <typename T>
	T L1
	(
		const T& x
	)
	{
		return x;
	}
	//The following is a general functoin that returns the value of the Legendre Polynomial for any given x and n=0,1,2,3,...

	template <typename T>
	T Legendre
	(
		unsigned int n,
		const T& x
	)
	{
		if (n == 0) {
			return L0(x);
		}
		else if (n == 1) {
			return L1(x);
		}
		else {
			return (double(2 * n - 1) * x * Legendre(n - 1, x) - double(n - 1) * Legendre(n - 2, x)) / double(n);
		}
	}
	////////

	////////

	template<typename T>
	T H0(T x) {
		return 1;
	}

	template<typename T>
	T H1(T x) {
		return 2 * x;
	}

	//The following is a general function that returns the value of the Hermite Polynomial for any given x and n=0,1,2,3,...
	template<typename T>
	T Hermite(int n, T x) {
		if (n == 0) {
			return H0(x);
		}
		else if (n == 1) {
			return H1(x);
		}
		else {
			return 2 * x * Hermite(n - 1, x) - 2 * (n - 1) * Hermite(n - 2, x);
		}
	}

	template<typename T>
	T La0
	(
		const T& x
	)
	{
		return 1;
	}

	template<typename T>
	T La1
	(
		const T& x
	)
	{
		return -x + 1;
	}

	//The following is a general function that returns the value of the Laguerre Polynomial for any given x and n=0,1,2,3,...
	template<typename T>
	T Laguerre
	(
		unsigned int n,
		const T& x
	)
	{
		if (n == 0) {
			return La0(x);
		}
		else if (n == 1) {
			return La1(x);
		}
		else {
			return ((2 * n - x - 1) * Laguerre(n - 1, x) - ((n - 1) * Laguerre(n - 2, x))) / n;
		}
	}

	template<typename T>
	T aLa0
	(
		unsigned int m,
		const T& x
	)
	{
		return 1;
	}

	template<typename T>
	T aLa1
	(
		unsigned int m,
		const T& x
	)
	{
		return -x + 1 + m;
	}

	//The following is a general function that returns the value of the assoc_Laguerre Polynomial for any given x and n=0,1,2,3,...
	template<typename T>
	T assoc_Laguerre
	(
		unsigned int n,
		unsigned int m,
		const T& x
	)
	{
		if (n == 0) {
			return aLa0(m, x);
		}
		else if (n == 1) {
			return aLa1(m, x);
		}
		else {
			return ((2 * n + m - x - 1) * assoc_Laguerre(n - 1, m, x) - ((n - 1 + m) * assoc_Laguerre(n - 2, m, x))) / n;
		}
	}

	template<typename T>
	T C0
	(
		const T& x
	)
	{
		return 1;
	}

	template<typename T>
	T C1
	(
		const T& x
	)
	{
		return x;
	}

	template<typename T>
	T Chebyshev
	(
		unsigned int n,
		const T& x
	)
	{
		if (n == 0) {
			return C0(x);
		}
		else if (n == 1) {
			return C1(x);
		}
		else {
			return 2 * x * Chebyshev(n - 1, x) - Chebyshev(n - 2, x);
		}
	}

	template<typename T>
	T Hp0
	(
		const T& x
	)
	{
		return 1;
	}

	template<typename T>
	T Hp1
	(
		const T& x
	)
	{
		return 2 * x;
	}

	template<typename T>
	T falling_factorial
	(
		const T& x,
		const unsigned int& n
	)
	{
		T c = 1;
		for (unsigned int k = 0; k < n; k++)
			c *= x - k;
		return c;
	}

	template<typename T>
	T rising_factorial
	(
		const T& x,
		const unsigned int& n
	)
	{
		T c = 1;
		for (unsigned int k = 0; k < n; k++)
			c *= x + k;
		return c;
	}

	template<typename T>
	T lower_gamma
	(
		const T& s,
		const T& z
	)
	{
		T c = 0;
		for (unsigned int k = 0; k < 30; k++)
			c += pow(z, k) / rising_factorial(s, k + 1);
		return pow(z, s) * exp(-z) * c;
	}

	template<typename T>
	T upper_gamma
	(
		const T& s,
		const T& z
	)
	{
		return tgamma(s) - lower_gamma(s, z);
	}

	template<typename T>
	void driver2
	(
		const T& d1,
		const T& x
	)
	{
		std::cout << "\n\n Polynomial test \n\n";

		std::cout << std::legendre(7, d1) << std::endl << std::endl;
		std::cout << Legendre(7, x) << std::endl;
		// cout << std::legendre(7, x) << endl << endl;

		std::cout << Hermite(7, x) << std::endl;
		//cout << std::hermite(7, x) << endl << endl;

		std::cout << Laguerre(7, x) << std::endl;
		//cout << std::laguerre(7, x) << endl << endl;

		std::cout << assoc_Laguerre(7, 2, x) << std::endl;
		//cout << std::assoc_laguerre(7, 2, x) << endl << endl;

		std::cout << Chebyshev(7, x) << std::endl;
		//cout << boost::math::chebyshev_t(7, x) << endl << endl;
	}

	template<typename T>
	T sph_bessel
	(
		const unsigned& v,
		const T& x
	)
	{
		T c = 0;
		const T pi = 3.141592653589793238462643383279502884197169399375105820974;

		for (unsigned k = 0; k < 20; k++) {
			c += (T(pow(-1, k)) * pow(x / 2., T(2 * k + v))) /
				factorial<T>(k) * tgamma(v + k + 1 + 0.5);
		}

		return sqrt(pi / 4) * c;
	}

	template<typename T>
	T cyl_bessel_i
	(
		const T& v,
		const T& x
	)
	{
		T c = 0;

		for (unsigned k = 0; k < 20; k++) {
			c += pow(x / 2., 2. * k + v) / factorial<T>(k) * tgamma(v + k + 1);
		}

		return c;
	}

	template<typename T>
	T cyl_bessel_j
	(
		const T& v,
		const T& x
	)
	{
		T c = 0;

		for (unsigned k = 0; k < 20; k++) {
			c += (T(pow(-1, k)) * pow(x / 2., 2 * k + v)) /
				factorial<T>(k) * tgamma(v + k + 1);
		}

		return c;
	}

	template<typename T>
	const T pi = 3.141592653589793238462643383279502884197169399375105820974;

	template<typename T>
	T cyl_bessel_k
	(
		const T& v,
		const T& x
	)
	{
		return (pi<T> / 2) * ((cyl_bessel_i(-v, x) - cyl_bessel_i(v, x)) / sin(v * pi<T>));
	}

	template<typename T>
	T cyl_neumann
	(
		const T& v,
		const T& x
	)
	{
		return (cyl_bessel_j(v, x) * cos(v * pi<T>) - cyl_bessel_j(-v, x)) / sin(v * pi<T>);
	}

	//////////////////

	template <typename T>
	void bessel_eval()
	{

		double xb = 1.2345;
		std::cout << "j_1(" << xb << ") = " << std::sph_bessel(1, xb) << '\n';

		std::cout << "j_1(" << xb << ") = " << sph_bessel(1, xb) << "\n\n";


		std::cout << "J_0(" << xb << ") = " << std::cyl_bessel_j(0, xb) << '\n';
		std::cout << "J_0(" << xb << ") = " << cyl_bessel_j(0., xb) << "\n\n";

		std::cout << "I_0(" << xb << ") = " << std::cyl_bessel_i(0, xb) << '\n';
		std::cout << "I_0(" << xb << ") = " << cyl_bessel_i(0., xb) << "\n\n";

		std::cout << "N_.5(" << xb << ") = " << std::cyl_neumann(0.5, xb) << '\n';
		std::cout << "N_.5(" << xb << ") = " << cyl_neumann(0.5, xb) << "\n\n";

		std::cout << "K_.5(" << xb << ") = " << std::cyl_bessel_k(.5, xb) << '\n';
		std::cout << "K_.5(" << xb << ") = " << cyl_bessel_k(0.5, xb) << "\n\n";
	}

}