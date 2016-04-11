//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include <iostream>
#include <time.h>
#include <conio.h>
#include <sstream>
//---------------------------------------------------------------------------
#include "lazymatr.h"
//---------------------------------------------------------------------------
//#define GLM_TEST
//---------------------------------------------------------------------------
#ifdef GLM_TEST
#	include "glm\glm\glm.hpp"
//---------------------------------------------------------------------------
std::ostream& operator << (std::ostream& s, const glm::mat4x4& m) {

	s << '[' << 4 << 'x' << 4 << ']' << '\n';

	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j)
			s << m[i][j] << " ";

		s << '\n';
	}

	return s;
}
#endif
//---------------------------------------------------------------------------
namespace tst {

	template<class T>
	struct mat4x4 {

		mat4x4() {}

		mat4x4(
			T a00, T a01, T a02, T a03,
			T a10, T a11, T a12, T a13,
			T a20, T a21, T a22, T a23,
			T a30, T a31, T a32, T a33)
			: m{ a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33 }
		{}	

		mat4x4 transpose() const {

			return mat4x4(
				m[0], m[4], m[8 ], m[12],
				m[1], m[5], m[9 ], m[13], 
				m[2], m[6], m[10], m[14], 
				m[3], m[7], m[11], m[15]
			);
		}
		
		T m[16];
	};

	template<class T>
	mat4x4<T> operator*(mat4x4<T> const& a, mat4x4<T> const& b) {

		return mat4x4<T>(
			a.m[ 0] * b.m[ 0] + a.m[ 1] * b.m[ 4] + a.m[ 2] * b.m[ 8] + a.m[ 3] * b.m[12],
			a.m[ 0] * b.m[ 1] + a.m[ 1] * b.m[ 5] + a.m[ 2] * b.m[ 9] + a.m[ 3] * b.m[13],
			a.m[ 0] * b.m[ 2] + a.m[ 1] * b.m[ 6] + a.m[ 2] * b.m[10] + a.m[ 3] * b.m[14],
			a.m[ 0] * b.m[ 3] + a.m[ 1] * b.m[ 7] + a.m[ 2] * b.m[11] + a.m[ 3] * b.m[15],

			a.m[ 4] * b.m[ 0] + a.m[ 5] * b.m[ 4] + a.m[ 6] * b.m[ 8] + a.m[ 7] * b.m[12],
			a.m[ 4] * b.m[ 1] + a.m[ 5] * b.m[ 5] + a.m[ 6] * b.m[ 9] + a.m[ 7] * b.m[13],
			a.m[ 4] * b.m[ 2] + a.m[ 5] * b.m[ 6] + a.m[ 6] * b.m[10] + a.m[ 7] * b.m[14],
			a.m[ 4] * b.m[ 3] + a.m[ 5] * b.m[ 7] + a.m[ 6] * b.m[11] + a.m[ 7] * b.m[15],

			a.m[ 8] * b.m[ 0] + a.m[ 9] * b.m[ 4] + a.m[10] * b.m[ 8] + a.m[11] * b.m[12],
			a.m[ 8] * b.m[ 1] + a.m[ 9] * b.m[ 5] + a.m[10] * b.m[ 9] + a.m[11] * b.m[13],
			a.m[ 8] * b.m[ 2] + a.m[ 9] * b.m[ 6] + a.m[10] * b.m[10] + a.m[11] * b.m[14],
			a.m[ 8] * b.m[ 3] + a.m[ 9] * b.m[ 7] + a.m[10] * b.m[11] + a.m[11] * b.m[15],

			a.m[12] * b.m[ 0] + a.m[13] * b.m[ 4] + a.m[14] * b.m[ 8] + a.m[15] * b.m[12],
			a.m[12] * b.m[ 1] + a.m[13] * b.m[ 5] + a.m[14] * b.m[ 9] + a.m[15] * b.m[13],
			a.m[12] * b.m[ 2] + a.m[13] * b.m[ 6] + a.m[14] * b.m[10] + a.m[15] * b.m[14],
			a.m[12] * b.m[ 3] + a.m[13] * b.m[ 7] + a.m[14] * b.m[11] + a.m[15] * b.m[15]
		);
	}

	template<class T>
	mat4x4<T> operator-(mat4x4<T> const& a, mat4x4<T> const& b) {

		return mat4x4<T>(
			a.m[ 0] - b.m[0],
			a.m[ 1] - b.m[1],
			a.m[ 2] - b.m[2],
			a.m[ 3] - b.m[3],
			a.m[ 4] - b.m[4],
			a.m[ 5] - b.m[5],
			a.m[ 6] - b.m[6],
			a.m[ 7] - b.m[7],
			a.m[ 8] - b.m[8],
			a.m[ 9] - b.m[9],
			a.m[10] - b.m[10],
			a.m[11] - b.m[11],
			a.m[12] - b.m[12],
			a.m[13] - b.m[13],
			a.m[14] - b.m[14],
			a.m[15] - b.m[15]
		);
	}

	template<class T>
	mat4x4<T> operator+(mat4x4<T> const& a, mat4x4<T> const& b) {

		return mat4x4<T>(
			a.m[0] + b.m[0],
			a.m[1] + b.m[1],
			a.m[2] + b.m[2],
			a.m[3] + b.m[3],
			a.m[4] + b.m[4],
			a.m[5] + b.m[5],
			a.m[6] + b.m[6],
			a.m[7] + b.m[7],
			a.m[8] + b.m[8],
			a.m[9] + b.m[9],
			a.m[10] + b.m[10],
			a.m[11] + b.m[11],
			a.m[12] + b.m[12],
			a.m[13] + b.m[13],
			a.m[14] + b.m[14],
			a.m[15] + b.m[15]
		);
	}

	template<class T>
	std::ostream& operator << (std::ostream& s, mat4x4<T> const& m) {

		s << '[' << 4 << 'x' << 4 << ']' << '\n';

		for (size_t i = 0; i < 4; ++i) {
			for (size_t j = 0; j < 4; ++j)
				s << m.m[i * 4 + j] << " ";

			s << '\n';
		}

		return s;
	}
}
//---------------------------------------------------------------------------
template<class M>
double performance_test() {

	const int size = 4 * 1024;
	M m[size], m1[size], m2;

	for (int i = 0; i < size; ++i)
		m1[i] = { 0.04f*float(rand() % size), 0.01f, 0.2f, 0.39f, -0.09f, -0.03f, 0.1f, 0.17f, -0.32f, -0.02f, 0.07f, 0.1f, 0.21f, -0.23f, -0.13f, -0.05f };

	m2 = { 0.02f, 0.01f, 0.2f, 0.39f, -0.09f, -0.03f, 0.1f, 0.17f, -0.32f, -0.02f, 0.07f, 0.1f, 0.21f, -0.23f, -0.13f, -0.05f };

	clock_t cl1, cl2;

	cl1 = clock();

	for (int k = 0; k < size; ++k)
		for (int i = 0; i < size; ++i)
			m[i] = m1[i] * m2 - m1[i] + m2;

	cl2 = clock();

	double sec = double(cl2 - cl1) / CLOCKS_PER_SEC;

	// Nflop = 2*(N^3) + 16 + 16
	double gflops = (size*size*(4.0*4.0*4.0*2.0 + 32)) / (sec * 1e+9);

	std::cout << m[rand() % size];

	return gflops;
}
//---------------------------------------------------------------------------
template<class M>
void mul_test() {

	M m;

	M m1{ 0.04f, 0.01f, 0.2f, 0.39f, -0.09f, -0.03f, 0.1f, 0.17f, -0.32f, -0.02f, 0.07f, 0.1f, 0.21f, -0.23f, -0.13f, -0.05f };
	M m2{ 0.02f, 0.01f, 0.2f, 0.39f, -0.09f, -0.03f, 0.1f, 0.17f, -0.32f, -0.02f, 0.07f, 0.1f, 0.21f, -0.23f, -0.13f, -0.05f };
	M res{ 0.0178f, -0.0936f, -0.0277f, 0.0178f, 0.0046f, -0.0411f, -0.0361f, -0.0387f, -0.006f, -0.027f, -0.0741f, -0.1262f, 0.056f, 0.0231f, 0.0164f, 0.0323f };

	m = m1*m2;

	//glm не прав ;) см. https://www.wolframalpha.com/input/?i=%7B%7B0.04,0.01,0.2,0.39%7D,%7B-0.09,-0.03,0.1,0.17%7D,%7B-0.32,-0.02,0.07,0.1%7D,%7B0.21,-0.23,-0.13,-0.05%7D%7D*%7B%7B0.02,0.01,0.2%09,0.39%7D,%7B-0.09,-0.03,0.1,0.17%7D,%7B-0.32,-0.02,0.07,0.1%7D,%7B0.21,-0.23,-0.13,-0.05%7D%7D
	std::cout << m << std::endl;
}
//---------------------------------------------------------------------------
///векторные и матричные операции на основе отложенных вычислений
void usage() {
	
	using namespace mtr;

	///1. инициализация, вывод в поток, сравнение
	{
		///если кол-во значений не совпадает с размерностю матрицы или вектора
		///будет ошибка компиляции
		mat3x3f m = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f };

		///можно вытащить строку или столбец
		std::cout << m.get_row(0) << m.get_col(1) << std::endl;

		mat1x3f v1 = m.get_row(0);
		mat1x3f v2 = m.get_row(1);

		// инициализация из вектров-строк
		mat2x3f mm{v1, v2};

		vec2f v3 = mm.get_col(0);
		vec2f v4 = mm.get_col(1);
		vec2f v5 = mm.get_col(2);

		// инициализация из вектров-столбцов
		mat2x3f mmm{v3, v4, v5 };

		// range based for
		for (auto& val : m)
			val *= 1.0f;

		mat3x3f m2; vec4i v{1, 2, 3, 4};

		std::cout << m << v << std::endl;
		///заполнить нулями все элементы
		m.zero();	v.zero();
		std::cout << m << v << std::endl;
		///еденичная матрица
		m.identity();
		std::cout << m << std::endl;

		///определены операции для поэлементного сравнения матриц == и != 
		m2 = { 1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.5f, 8.6f, 9.0f };

		if (m != m2)
			std::cout << "m != m2" << std::endl;

		///определены операции для поэлементного сравнения с константой > < и ф-ия abs_lt/abs_gt (меньше/больше по абсолютному значение)
		if (m > 0.0f && m2 < 15.0f)
			std::cout << "m > 0.0f && m2 < 15.0f" << std::endl;

	}

	///2. Векторные операции
	///Соответствие размерностей всех векторов в выражениях проверяются в compile time
	{
		vec3f c, d;

		auto a = vec(1.0f, 2.0f, 3.0f);
		auto b = vec(5.0f, 6.0f, 7.0f);

		float dot = a | b;	///<скалярное произведение
		c = a ^ b;			///<векторное произведение		

		d = 5.0f * (c ^ a) + (a | c) - a * b;

		///все операции кроме векторного произведения определены
		///для векторов произвольной размерности
		vec10f e, f, g;

		e << 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f;
		f = e;

		dot = e | f;
		g = e + f - dot;

		///сумма квадратов элементов вектора
		float s = e.sqr_sum();
		///длина вектора
		float l = f.length();
		///нормализация
		e = normalize(g);///<нормализует g
	}

	///3. Матричные операции
	///Соответствие размерностей всех матриц в выражениях проверяются в compile time
	{
		mat3x3f m3x3, a, b;
		mat2x3f m2x3{ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f };
		mat3x2f m3x2{ 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f };

		///умножение
		///!!!!!Не иcпользуйте больше одной операции матричного 
		///умножения в выражении, очень большой overhead! 
		///Подробности http://www.gamedev.ru/forum/?group=2&topic=1661&page=3 пост 42
		m3x3 = 3.0f * m3x2 * m2x3;

		a = { 9.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f };
		b = { 9.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f };

		///обращение матриц
		m3x3 = b * inverse(a);///<обращает a и возвращает ее по ссылке. Проверка на сингулярность assert-ом
		std::cout << "Inversion\n" << m3x3 << std::endl;

		if (inverse_if(a))///<обращает a если определитель не ноль и возвращает false, в противном случае возвращает true.
			std::cout << "Zero determinant!" << std::endl;

		///для матриц 4х4, кот. являются комбинацией поворотов и сдвигов можно использовать
		///более быстрый алгоритм обращения, использующий  эту особенность
		mat4x4f ma, mb, res;

		ma << 0.57487148f, 0.20136166f, 0.79308027f, 0.0f,
			-0.78165251f, 0.42174986f, 0.45950660f, 0.0f,
			-0.24195445f, -0.88407046f, 0.39984691f, 0.0f,
			0.0f, 0.0f, -316.22778f, 1.0f;

		mb = ma;

		res = inverse_affine(ma) * mb;///< res - единичная матрица

		std::cout << "Affine inversion\n" << res << std::endl;

		///транспонирование		
		mat3x2f m1{ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f };

		mat2x3f m2 = m1.transpose();

		std::cout << m1 << m2 << std::endl;

		///составные выражения
		m3x3 = 3.0f * m3x2 * m2x3 - inverse(a) + b.transpose();
	}

	///Можно определить матрицу или вектор произвольной размерности 
	typedef mtr::Vector<20, long double> vec20ld;
	typedef mtr::Matrix<15, 4, unsigned int> mat15x4ui;
}

//---------------------------------------------------------------------------
int _tmain(int argc, _TCHAR* argv[]) {	
	
	usage();

	mul_test<mtr::mat4x4f>();
	mul_test<tst::mat4x4<float>>();

#ifdef GLM_TEST
	mul_test<glm::mat4x4>();
#endif

	auto gf_mtr = performance_test<mtr::mat4x4f>();
	auto gf_tst = performance_test<tst::mat4x4<float>>();

	std::cout << "\n\nlazy performance: " << gf_mtr << " GFLOPS\n";	
	std::cout << "standart performance: " << gf_tst << " GFLOPS\n\n";
	std::cout << "lazy faster then standart: " << gf_mtr / gf_tst << " times\n\n";

#ifdef GLM_TEST	
	auto gf_glm = performance_test<glm::mat4x4>();
	std::cout << "glm performance: " << gf_glm << " GFLOPS\n\n";
	std::cout << "lazy faster then glm: " << gf_mtr / gf_glm << " times";
#endif
	
	_getch();
	return 0;
}