//---------------------------------------------------------------------------
/// Written by Pavel Kazmerchuk 2005 - 2006.
/// Refactored to C++14 by Pavel Kazmerchuk 2016
/// Algorithm of matrices inversion - Федор 
/// http://www.gamedev.ru/community/cfd/articles/RecursiveMethod_
//---------------------------------------------------------------------------
#pragma once
//---------------------------------------------------------------------------
#include <assert.h>
#include <math.h>
#include <memory.h>
//---------------------------------------------------------------------------
#include <ostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <tuple>
//---------------------------------------------------------------------------
#ifdef _MSC_VER
#	pragma inline_depth( 255 )
#	pragma inline_recursion( on )
#	pragma auto_inline( on )
#	define inline __forceinline
#endif
//---------------------------------------------------------------------------
/// Lazy Matrix
//---------------------------------------------------------------------------
namespace mtr{
//---------------------------------------------------------------------------
// Declaration
//---------------------------------------------------------------------------
template<typename TLeftArg, typename TRightArg, typename TOp> class Expression;
template<size_t row, typename T> class VectorImpl;
template<size_t row, size_t col, typename T> class MatrixImpl;
template<typename M> class MatrixBase;
//---------------------------------------------------------------------------
template<size_t row, size_t col, typename T>
using Matrix = MatrixBase<MatrixImpl<row, col, T>>;
//---------------------------------------------------------------------------
template<size_t row, typename T>
using Vector = MatrixBase<VectorImpl<row, T>>;
//---------------------------------------------------------------------------
/// Helper classes	
//---------------------------------------------------------------------------
template<size_t f1, size_t f2, typename T1, typename T2>
struct Select2 {

	typedef T1 Result;
};
//---------------------------------------------------------------------------
template<size_t f, typename T1, typename T2>
struct Select2<0, f, T1, T2> {

	typedef T2 Result;
};
//---------------------------------------------------------------------------
template<size_t f, typename T1, typename T2>
struct Select2<f, 0, T1, T2> {

	typedef T2 Result;
};
//---------------------------------------------------------------------------
// unroll
//---------------------------------------------------------------------------
///< s - number of iterations, n - current iteration
template<size_t s, size_t n> struct Recurse;
//---------------------------------------------------------------------------
///< Inverse unroll
template<int j, int k>
struct Pivot_jk_If {

	template< int i, typename T >
	inline static void pivot_sub( T& a ) {

		a(i,j) -= a(i,k)*a(k,j);
	}

	template< int i, typename U, typename T >
	inline static void pivot_sub( U const& v, T& a ) {

		a(k,i) -= v(j)*a(j,i);
	}

	template< int i, typename U, typename T >
	inline static void pivot_sub2( U const& v, T& a ) {

		a(i,k) -= v(j)*a(i,j);
	}
};
//---------------------------------------------------------------------------
template<int k>
struct Pivot_jk_If<k,k> {

	template< int i, typename T >
	inline static void pivot_sub( T& a ) {
	}	

	template< int i, typename U, typename T >
	inline static void pivot_sub( U const& v, T& a ) {
	}

	template< int i, typename U, typename T >
	inline static void pivot_sub2( U const& v, T& a ) {
	}
};
//---------------------------------------------------------------------------
template<int i, int k>
struct Pivot_ik_If {

	template<int nDim, typename T>
	inline static void pivot_sub( T& a ) {

		return Recurse<nDim, 0>::template pivot2_sub<k,i>(a);
	}

	template<int nDim, typename U, typename T>
	inline static void pivot_sub( U const& v, T& a ) {

		return Recurse<nDim, 0>::template pivot2_sub<k,i>(v, a);
	}

	template<int nDim, typename U, typename T>
	inline static void pivot_sub2( U const& v, T& a ) {

		return Recurse<nDim, 0>::template pivot2_sub2<k,i>(v, a);
	}
};
//---------------------------------------------------------------------------
template<int k>
struct Pivot_ik_If<k,k> {

	template<int nDim, typename T>
	inline static void pivot_sub( T& a ) {		
	}

	template<int nDim, typename U, typename T>
	inline static void pivot_sub( U const& v, T& a ) {
	}

	template<int nDim, typename U, typename T>
	inline static void pivot_sub2( U const& v, T& a ) {
	}
};
//---------------------------------------------------------------------------
template<size_t s, size_t n>///< s - number of iterations, n - current iteration
struct Recurse {

	enum { count = n+1 };

	///mathematical matrix multiplication unroll
	template< typename TLeft, typename TRight >
	inline static auto mul( const TLeft& A, const TRight& B, size_t i, size_t j ) {

		enum {
			rows = TLeft::cols,
			cols = TRight::cols
		};

		return A( i*(size_t)rows + n ) * B( n*(size_t)cols + j) + Recurse<s, count>::mul( A, B, i, j );		
	}

	///dot product unroll
	template< typename TLeft, typename TRight >
	inline static auto dot( const TLeft& A, const TRight& B ) {		

		return A(n) * B(n) + Recurse<s, count>::dot( A, B );
	}

	///normalize unroll
	template< typename T >
	inline static void normalize( T& a, const typename T::Element& b ) {

		a(n) *= b;
		Recurse<s, count>::normalize( a, b );
	}

	///get row, get column unroll
	template< size_t shift, typename TLeft, typename TRight >
	inline static void getrowcol(size_t i, TLeft& A, const TRight& B) {	

		A(n) = B(i);
		Recurse<s, count>::template getrowcol<shift>(shift + i, A, B);		
	}

	///+=
	template< typename T >
	inline static void plus_eq( T& a, const typename T::Element& b ) {		

		a(n) += b;
		Recurse<s, count>::plus_eq( a, b );	
	}

	template< typename T >
	inline static void plus_eq( T& a, T const& b ) {	

		a(n) += b(n);
		Recurse<s, count>::plus_eq( a, b );	
	}

	template< typename T >
	inline static void minus_eq( T& a, T const& b ) {		

		a(n) -= b(n);
		Recurse<s, count>::minus_eq( a, b );	
	}

	template< typename T >
	inline static void mul_const( T& a, const typename T::Element& b ) {

		a(n) *= b;
		Recurse<s, count>::mul_const( a, b );	
	}

	template< typename T >
	inline static void incrFill( T& a ) {	

		a(n) = n;
		Recurse<s, count>::incrFill( a );	
	}

	///< Inverse unroll
	template<int nDim, int k, typename T>
	inline static void pivot1_sub( T& a ) {

		Pivot_ik_If<n,k>::template pivot_sub<nDim>(a),
			Recurse<s, count>::template pivot1_sub<nDim, k>(a);
	}

	template< int k, int i, typename T>
	inline static void pivot2_sub( T& a ) {

		Pivot_jk_If<n,k>::template pivot_sub<i>(a),
			Recurse<s, count>::template pivot2_sub<k, i>(a);
	}

	template<int nDim, int k, typename U, typename T>
	inline static void pivot1_sub( U const& v, T& a ) {

		a(k,n)=0;
		Pivot_ik_If<n,k>::template pivot_sub<nDim>(v, a),
			Recurse<s, count>::template pivot1_sub<nDim, k>(v, a);
	}

	template< int k, int i, typename U, typename T>
	inline static void pivot2_sub( U const& v, T& a ) {

		Pivot_jk_If<n,k>::template pivot_sub<i>(v, a),
			Recurse<s, count>::template pivot2_sub<k, i>(v, a);
	}

	template<int nDim, int k, typename U, typename T>
	inline static void pivot1_sub2( U const& v, T& a ) {

		a(n,k)=0;
		Pivot_ik_If<n,k>::template pivot_sub2<nDim>(v, a),
			Recurse<s, count>::template pivot1_sub2<nDim, k>(v, a);
	}

	template< int k, int i, typename U, typename T>
	inline static void pivot2_sub2( U const& v, T& a ) {

		Pivot_jk_If<n,k>::template pivot_sub2<i>(v, a),
			Recurse<s, count>::template pivot2_sub2<k, i>(v, a);
	}

	template<typename T>
	inline static void add1(int noZero, T& a ) {

		a(n, noZero) += a(n, s);
		Recurse<s, count>::add1(noZero, a);
	}

	template<typename T>
	inline static void add2(int noZero, T& a ) {

		a(s, n) += a(noZero, n);
		Recurse<s, count>::add2(noZero, a);
	}

	template<typename T>
	inline static void mul1(T& a ) {

		a(s, n) *= a(s, s);
		Recurse<s, count>::mul1(a);
	}

	template<typename U, typename T>
	inline static void mul2(U& v, T& a ) {

		v(n) = a(s, n);
		a(s, s) -= v(n)*a(n, s);
		Recurse<s, count>::mul2(v, a);
	}

	template<typename U, typename T>
	inline static void mul3(U& v, T& a ) {

		v(n) = a(n,s) * a(s,s);
		Recurse<s, count>::mul3(v, a);
	}

	template<typename T>
	inline static void findPivot(int& noZero, typename T::Element& pivot, T const& a ) {

		const typename T::Element tmp = abs(a(n,s));
		if (pivot < tmp)
			pivot = tmp, noZero = n;
		Recurse<s, count>::findPivot(noZero, pivot, a);
	}
	///< Inverse unroll end
};
//---------------------------------------------------------------------------
template<size_t s>
struct Recurse<s, s> {

	template< typename TLeft, typename TRight >
	inline static auto mul( const TLeft&, const TRight&,  size_t, size_t) {

		return 0;
	}

	template< typename TLeft, typename TRight >
	inline static auto dot( const TLeft&, const TRight&) {		

		return 0;
	}

	template< typename T >
	inline static void normalize( const T&, const typename T::Element&) {		
	}

	template< size_t shift, typename TLeft, typename TRight >
	inline static void getrowcol(size_t, TLeft&, const TRight&) {		
	}

	template< typename T >
	inline static void plus_eq(T&, const typename T::Element&) {		
	}

	template< typename T >
	inline static void plus_eq(T&, T const&) {
	}

	template< typename T >
	inline static void minus_eq( T&, T const&) {		
	}	
	
	template< typename T >
	inline static void mul_const(T&, const typename T::Element&) {		
	}

	template< typename T >
	inline static void incrFill(T&) {		
	}

	///< Inverse unroll
	template< int nDim, int k, typename T>
	inline static void pivot1_sub(T&) {
	}
	template< int nDim, int k, typename T>
	inline static void pivot2_sub(T&) {
	}
	template< int nDim, int k, typename U, typename T>
	inline static void pivot1_sub(U&, T&) {
	}
	template< int nDim, int k, typename U, typename T>
	inline static void pivot2_sub(U&, T&) {
	}
	template<int nDim, int k, typename U, typename T>
	inline static void pivot1_sub2(U&, T&)  {
	}	
	template< int k, int i, typename U, typename T>
	inline static void pivot2_sub2(U&, T&) {
	}
	template<typename T>
	inline static void add1(int, T&) {
	}
	template<typename T>
	inline static void add2(int, T&) {
	}
	template<typename T>
	inline static void mul1(T&) {
	}
	template<typename U, typename T>
	inline static void mul2(U&, T&) {
	}
	template<typename U, typename T>
	inline static void mul3(U&, T&) {
	}
	template<typename T>
	inline static void findPivot(int&, typename T::Element&, T const&) {
	}
	///< Inverse unroll end
};
//---------------------------------------------------------------------------
/// Assigment unroll
//---------------------------------------------------------------------------
template<size_t I, size_t dim>
struct Assigment
{
	enum { count = I+1 };

	template<typename T, typename Exp>
	inline static void assign( T& V, const Exp& A) {

		V(I) = A(I);
		Assigment<count, dim>::assign(V, A);
	}
};
//---------------------------------------------------------------------------
template<size_t dim>
struct Assigment<dim, dim>
{
	template<typename T, typename Exp>     
	inline static void assign( T& , const Exp&  ){}
};

//---------------------------------------------------------------------------
/// Initialization.  Example usage: Matrix2x2f a;  a << 1,2,3,4;
//---------------------------------------------------------------------------
template<typename TLeft, typename TRight, size_t I, size_t N>
class Init
{
	TLeft& left_;
	Init& operator=(Init const&) = delete;

public:

	Init(TLeft& left, const TRight& right)
		: left_(left) {	

		static_assert( I < N, "Wrong number of initializers." );
		left_(I) = right;
	}

	Init<TLeft, TRight, I+1, N> operator ,(const TRight v) {
		return Init<TLeft, TRight, I+1, N>(left_, v);
	}
};
//---------------------------------------------------------------------------
/// Arguments
//---------------------------------------------------------------------------
template<typename T>
class Argument
{
	typedef T ArgType;
	const ArgType& arg_; 
	
	Argument& operator=(const Argument&) = delete;

public:

	enum {
		rows = ArgType::rows, 
		cols = ArgType::cols 
	};

	typedef typename ArgType::Element Element;

	Argument(const ArgType& arg) : arg_(arg) {}

	inline const Element operator()(size_t i) const {

		return arg_(i);
	}

	inline bool check_down(const void* p) const { 

		return arg_.check_down(p);
	}
};
//---------------------------------------------------------------------------
#define SCALAR(type) \
template<>\
class Argument<const type> {\
const type arg_;\
Argument& operator=(const Argument&) = delete;\
public:\
enum {\
rows = 0,\
cols = 0 \
};\
typedef type Element;\
Argument(const type& arg) : arg_(arg) {}\
inline Element operator()(size_t) const {\
return arg_;\
}\
inline bool check_down(const void*) const {\
return true; \
}\
};
//---------------------------------------------------------------------------
SCALAR(double)
SCALAR(float)
SCALAR(int)
SCALAR(unsigned char)
SCALAR(unsigned short)
SCALAR(unsigned int)
//---------------------------------------------------------------------------
// supported operations
//---------------------------------------------------------------------------
namespace op_type {

	enum : int {
		sum,
		add,
		mul,
		mul_scalar,
		diff,
		unary_minus,
		dot,
		cross,
		div,
		transpose,
		none
	};
};
//---------------------------------------------------------------------------
namespace detail {

	template<bool, size_t n, size_t m>
	struct static_if {
		enum { value = n };
	};

	template<size_t n, size_t m>
	struct static_if<false, n, m> {
		enum { value = m };
	};

	constexpr bool eq(size_t n, size_t m) {
		return n == m;
	}

	template<class T>
	constexpr bool is_scalar() {
		return eq(T::rows, 0) && eq(T::cols, 0);
	}

	template<class T, class U>
	constexpr bool is_matrix_mul() {
		return eq(T::cols, U::rows);
	}

	/// left argumrnt is the same vector as a right argument
	template<class T, class U>
	constexpr bool is_same_vec() {
		return
			(eq(T::cols, U::cols) && eq(T::rows, U::rows) && eq(T::cols, 1)) ||
			(eq(T::cols, U::cols) && eq(T::rows, U::rows) && eq(T::rows, 1));
	}

	template<class T, size_t ot>
	constexpr bool is_op() {
		return eq(T::type, ot);
	}

	template<class T>
	constexpr bool is_transposing() {
		return is_op<T, op_type::transpose>();
	}

	template<class T, class U>
	struct common_type {
		typedef typename std::common_type<typename T::Element, typename U::Element>::type type;
	};

	template<class T, class U, class O>
	class define_size {
		enum {
			tmp_rows = static_if<is_scalar<T>(), U::rows, T::rows>::value,
			tmp_cols = static_if<is_scalar<U>(), T::cols, U::cols>::value
		};

	public:
		enum {			
			rows = static_if<is_transposing<O>(), tmp_cols, tmp_rows>::value,
			cols = static_if<is_transposing<O>(), tmp_rows, tmp_cols>::value
		};
	};
}
//---------------------------------------------------------------------------
// Expressions
//---------------------------------------------------------------------------
template< typename TLeftArg, typename TRightArg, typename TOp>
class Expression
{
	typedef Argument<TLeftArg>  TLeft;
	typedef Argument<TRightArg> TRight;

	const TLeft  left_;
	const TRight right_;	

public:

	enum {
		rows = detail::define_size<TLeft, TRight, TOp>::rows,
		cols = detail::define_size<TLeft, TRight, TOp>::cols,
		dim = rows * cols
	};
	
	typedef typename detail::common_type<TLeft, TRight>::type Element;
	typedef TOp Operation;

	Expression( const TLeft& left, const TRight& right )
		: left_(left)
		, right_(right) {

		// check dimensions
		static_assert(
			detail::is_matrix_mul<TLeft, TRight>() || 
			detail::is_scalar<TLeft>() || 
			detail::is_scalar<TRight>() || 
			detail::is_same_vec<TLeft, TRight>() ||
			detail::is_transposing<Operation>(),
			"Wrong dimensions of arguments."
		);
	}

	inline Element operator()(size_t i) const { 

		return TOp::evaluate(i, left_, right_); 
	}

	inline bool check_down(const void* p) const { 

		return left_.check_down(p) && right_.check_down(p) ;
	}
};
//---------------------------------------------------------------------------
///Operations implementation
//---------------------------------------------------------------------------
struct Void {

	enum {
		type = op_type::none
	};

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t i, const TLeft& A, TRight) {

		return A(i);
	}
};
//---------------------------------------------------------------------------
///Addition
//---------------------------------------------------------------------------
struct Sum {

	enum {
		type = op_type::sum
	}; 

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t i, const TLeft& A, const TRight& B ) { 

		return A(i) + B(i); 
	}
};
//---------------------------------------------------------------------------
///Difference
//---------------------------------------------------------------------------
struct Diff {

	enum {
		type = op_type::diff
	};

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t i, const TLeft& A, const TRight& B ) { 

		return A(i) - B(i); 
	}
};
//---------------------------------------------------------------------------
///Difference
//---------------------------------------------------------------------------
struct UnaryMinus {

	enum {
		type = op_type::unary_minus
	};

	template< typename T>
	inline static typename T::Element evaluate(size_t i, const T& A, const T& B ) { 

		return -A(i); 
	}
};
//---------------------------------------------------------------------------
///Dot product
//---------------------------------------------------------------------------
struct Dot {

	enum {
		type = op_type::dot
	};

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t, const TLeft& A, const TRight& B ) {

		enum{
			s = TLeft::rows * TLeft::cols
		};

		return Recurse<s, 0>::dot(A, B); 
	}
};
//---------------------------------------------------------------------------
///Cross product
//---------------------------------------------------------------------------
struct Cross {

	enum {
		type = op_type::cross
	};

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t i, const TLeft& A, const TRight& B ) { 

		const size_t l = (i + 1)%3;
		const size_t m = (i + 2)%3;

		return A(l)*B(m) - A(m)*B(l); 
	}
};
//---------------------------------------------------------------------------
///Multiplication by scalar
//---------------------------------------------------------------------------
struct MulConst {

	enum {
		type = op_type::mul_scalar
	};

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t i, const TLeft& A, const TRight& B ) { 

		return A(i) * B(i);  
	}
};
//---------------------------------------------------------------------------
///Mathematical matrix multiplication
//---------------------------------------------------------------------------
struct Mul {

	enum {
		type = op_type::mul
	};

	template< typename TLeft, typename TRight> 
	inline static auto evaluate(size_t k, const TLeft& A, const TRight& B ) { 

		enum {
			//result matrix rows number
			rows = TLeft::rows,
			//result matrix columns number			
			cols = TRight::cols,
			//s = TLeft::cols = TRight::row - number of iteration			
			s = TLeft::cols
		};

		return Recurse<s, 0>::mul(A, B, k/(size_t)cols, k - k/(size_t)cols * (size_t)cols);
	}
};
//---------------------------------------------------------------------------
///Transpose matrix
//---------------------------------------------------------------------------
struct Trans { 

	enum {
		type = op_type::transpose
	};

	template< typename TLeft, typename TRight> 
	inline static auto evaluate(size_t k, const TLeft& A, const TRight& ) { 			

		return A( (k - k/(size_t)TLeft::rows * (size_t)TLeft::rows) * (size_t)TLeft::cols + k/(size_t)TLeft::rows );
	}
};
//---------------------------------------------------------------------------
///Division
//---------------------------------------------------------------------------
struct Div {

	enum {
		type = op_type::div
	};

	template< typename TLeft, typename TRight>
	inline static auto evaluate(size_t i, const TLeft& A, const TRight& B ) { 

		return A(i) / B(i); 
	}
};
//---------------------------------------------------------------------------
///Matrix Inversion
//---------------------------------------------------------------------------
template <int iDim>
struct CRecurseInv {

	template <typename U, typename V>
	inline static bool InvertMatrixRecursionInUnroll(MatrixBase<U>& v, MatrixBase<V>& a) {

		typedef typename V::Element T;

		const T eps = std::numeric_limits<T>::epsilon(); 

		int iNoZero = iDim;		T Pivot(0);		

		if (eps > abs(a(iDim,iDim)))
			for (int i = 0; i < iDim; ++i)
			{
				const T tmp = abs(a(i,iDim));
				if (Pivot < tmp)
					Pivot = tmp, iNoZero = i;
			}

			if (eps > abs(a(iNoZero,iDim)))  return true; 

			if (iNoZero != iDim)   ///<a22==0
				Recurse<iDim,0>::add2(iNoZero, a);

			a(iDim,iDim) = static_cast<T>(1)/a(iDim,iDim);

			Recurse<iDim,0>::mul1(a); ///< a21 = a21/a22

			Recurse<iDim,0>::template pivot1_sub<iDim, iDim>(a);	///< a11 = a11 - a12*a21

			if (CRecurseInv<iDim-1>::InvertMatrixRecursionInUnroll (v, a))  //b11 = (a11 - a12*a21)^(-1)
				return true;

			Recurse<iDim,0>::mul3(v, a);///<a12 = a12/a22 		

			Recurse<iDim,0>::template pivot1_sub2<iDim, iDim>(v, a);///< b12 = -b11*a12

			Recurse<iDim,0>::mul2(v, a);///< b22 =1/a22 - a21*b12

			Recurse<iDim,0>::template pivot1_sub<iDim, iDim>(v, a); ///< b21 = -a21 * b11


			if (iNoZero != iDim)  ///<a22==0
				Recurse<iDim,0>::add1(iNoZero, a);

			return false;
	}
};
//---------------------------------------------------------------------------
template <>
struct CRecurseInv<0> {

	template <typename U, typename V>
	inline static bool InvertMatrixRecursionInUnroll(MatrixBase<U>& , MatrixBase<V>& a) {

		typedef typename V::Element T;
		//if a matrix has a zero determinant
		if (std::numeric_limits<T>::epsilon() > abs( a(0,0) ))  
			return true;

		a(0,0) = static_cast<T>(1)/a(0,0);
		return false;
	}
};
//---------------------------------------------------------------------------
///Vector
//---------------------------------------------------------------------------
template<size_t row, typename T>
class VectorImpl {

public:

	typedef T Element;

	enum { 
		rows = row, 
		cols = 1,
		dim  = rows*cols  
	};

	VectorImpl(){}

	template<class ...U>
	constexpr VectorImpl(U... next)
		: data_{ next... } 
	{
	}

	//-----------------------------------------------------------------------
	void zero() {

		memset(data_, 0, sizeof(data_));
	}
	//-----------------------------------------------------------------------
	///indexed access
	inline Element& operator()(size_t i) {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	inline Element const& operator()(size_t i) const {

		return data_[i];
	}
	inline Element& operator[](size_t i) {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	inline Element const& operator[](size_t i) const {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	/// dot product
	inline Element dot(VectorImpl<dim, T> const& m) const {

		return Dot::evaluate(0, *this, m);
	}
	//-----------------------------------------------------------------------
	/// dot product sugar
	inline Element operator | (VectorImpl<dim, T> const& m) const {

		return dot(m);
	}
	//-----------------------------------------------------------------------
	/// sum of squares of vector elements
	inline Element sqr_sum() const {

		return dot(*this);
	}
	//-----------------------------------------------------------------------
	/// length of vector
	inline Element length() const {

		return sqrt(sqr_sum());
	}
	//-----------------------------------------------------------------------
	/// normalization
	inline void normalize() {

		assert(length() != (Element)(0));
		Recurse<dim, 0>::normalize(*this, Element(1) / length());
	}

protected:	

	Element data_[dim];	
};
//---------------------------------------------------------------------------
template<size_t row, size_t col, typename T>
class MatrixImpl {

public:

	typedef T Element;

	enum { 
		rows = row, 
		cols = col,
		dim = row*col
	};

	MatrixImpl(){}

	template<class ...U>
	constexpr MatrixImpl(U... next)
		: data_{ next... }
	{
	}

	//-----------------------------------------------------------------------
	void zero() {

		memset(data_, 0, sizeof(data_));
	}
	//-----------------------------------------------------------------------
	void identity() {

		static_assert(row == col, "Matrix is not square.");

		zero(); 

		for(size_t i = 0; i < col; ++i)
			data_[i*col + i] = Element(1);
	}
	//-----------------------------------------------------------------------
	inline Element& operator()(size_t i) {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	inline Element const& operator()(size_t i) const {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	inline Element& operator[](size_t i) {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	inline Element const& operator[](size_t i) const {

		return data_[i];
	}
	//-----------------------------------------------------------------------
	inline Element& operator()(size_t i, size_t j) {

		return data_[i*col + j];
	}
	//-----------------------------------------------------------------------
	inline Element const& operator()(size_t i, size_t j) const {

		return data_[i*col + j];
	}
	//-----------------------------------------------------------------------
	///getting i-th row
	inline Matrix<1, cols, T> get_row(size_t i) const {

		assert(i < rows);

		Matrix<1,cols,T> res;

		Recurse<cols, 0>::template getrowcol<1>((size_t)cols*i, res, *this);
		return res;
	}
	//-----------------------------------------------------------------------
	///getting i-th column
	inline Vector<rows,T> get_col(size_t i) const {

		assert(i < cols);

		Vector<rows, T> res;

		Recurse<rows, 0>::template getrowcol<cols>(i, res, *this);
		return res;
	}

protected:	

	Element data_[dim];	
};
//---------------------------------------------------------------------------
template<std::size_t i>
struct Int2Type {
	enum {val = i};
};
//---------------------------------------------------------------------------
template<typename M>
class MatrixBase : public M {

	using M::data_;	

public:
	
	using M::operator();

	typedef MatrixBase<M> Self;    
	typedef typename M::Element Element;
	typedef typename M::Element value_type;

	enum { 
		rows = M::rows, 
		cols = M::cols,
		dim  = M::dim  
	};

	MatrixBase(){}	

	template<class ...T>
	constexpr MatrixBase(Element first, T... next)
		: M{ first, next... } {

		static_assert(sizeof...(next)+1 == dim, "Wrong number of initializers.");
	}

	// initialization from column vectors
	template<class ...T> 
	constexpr MatrixBase(const Matrix<1, cols, Element>& first, T... next)
		: MatrixBase{std::make_index_sequence<cols>{}, first, next... } {}

	// initialization from row vectors
	template<class ...T>
	constexpr MatrixBase(const Vector<rows, Element>& first, T... next)
		: MatrixBase{ Int2Type<0>(), std::make_tuple(first, next...), std::make_index_sequence<cols>{}} {}

private:
	
	template<std::size_t ...indices, class ...T>
	constexpr MatrixBase(std::integer_sequence<std::size_t, indices...> s, const Matrix<1, cols, Element>& first, T... next)
		: MatrixBase{ s, next... , first[indices]...} {}

	template<std::size_t ...indices, class ...T>
	constexpr MatrixBase(std::integer_sequence<std::size_t, indices...> s, Element first, T... next)
		: MatrixBase{ first, next... } {}

	template<std::size_t i, class U, std::size_t ...indices, class ...T>
	constexpr MatrixBase(Int2Type<i>, U const& t, std::integer_sequence<std::size_t, indices...> s, T... next)
		: MatrixBase{ Int2Type<i+1>(), t, s, next..., std::get<indices>(t)[i]... } {}

	template<class U, std::size_t ...indices, class ...T>
	constexpr MatrixBase(Int2Type<rows>, U const& t, std::integer_sequence<std::size_t, indices...> s, T... next)
		: MatrixBase{next... } {}

public: 

	//-----------------------------------------------------------------------
	template<class T, class U, class O>
	inline Self& operator = ( const Expression<T,U,O>& expression ) {

		typedef Expression<T,U,O> E;

		static_assert(detail::eq(cols, E::cols) && detail::eq(rows, E::rows), "Matrix assigment mismatch dimensions.");

		// check arguments
		assert(expression.check_down(begin()) && "Target of assigment must not be used in expression (on the left of = ).");

		Assigment<0,dim>::assign( *this, expression );
		return *this;
	}
	//-----------------------------------------------------------------------
	template<class T, class U, class O>
	inline MatrixBase(Expression<T, U, O> const& expression) {

		typedef Expression<T, U, O> E;

		static_assert(detail::eq(cols, E::cols) && detail::eq(rows, E::rows), "Matrix assigment mismatch dimensions.");

		Assigment<0, dim>::assign(*this, expression);
	}
	//-----------------------------------------------------------------------
	// return iterator for beginning
	inline const Element* begin() const {

		return &data_[0];
	}
	//-----------------------------------------------------------------------
	inline Element* begin() {

		return &data_[0];
	}
	//-----------------------------------------------------------------------
	// return iterator for end
	inline const Element* end() const {

		return &data_[dim];
	}
	//-----------------------------------------------------------------------
	inline Element* end() {

		return &data_[dim];
	}
	//-----------------------------------------------------------------------
	inline size_t size() const {

		return dim;
	}
	//-----------------------------------------------------------------------
	inline bool check_down(const void* p) const { 

		return begin() != p; 
	}
	//-----------------------------------------------------------------------
	// return proxy object for initialization through comma
	Init<Self, Element, 0, dim> operator << (Element v) {

		return Init<Self, Element, 0, dim>(*this, v);
	}
	//-----------------------------------------------------------------------
	inline Expression<const Self, const Self, Trans> transpose() const {

		return Expression<const Self, const Self, Trans>(*this, *this);
	}
	//-----------------------------------------------------------------------
	inline Self& operator += (Element e) {

		Recurse<dim, 0>::plus_eq(*this, e);
		return *this;
	}
	//-----------------------------------------------------------------------
	inline Self& operator *= (Element e) {

		Recurse<dim, 0>::mul_const(*this, e);
		return *this;
	}
	//-----------------------------------------------------------------------
	inline Self& operator /= (Element e) {

		Recurse<dim, 0>::mul_const(*this, static_cast<Element>(1)/e);
		return *this;
	}
	//-----------------------------------------------------------------------
	inline Self& operator -= (Element e) {

		Recurse<dim, 0>::plus_eq(*this, -e);
		return *this;
	}
	//-----------------------------------------------------------------------
	inline Self& operator += (Self const& m) {

		Recurse<dim, 0>::plus_eq(*this, m);
		return *this;
	}
	//-----------------------------------------------------------------------
	inline Self& operator -= (Self const& m) {

		Recurse<dim, 0>::minus_eq(*this, m);
		return *this;
	}
	//-----------------------------------------------------------------------
	inline Expression<Self, Self, UnaryMinus> operator - () const {

		return Expression<Self, Self, UnaryMinus>(*this, *this);
	}	
	//-----------------------------------------------------------------------
	bool operator == (Self const& mr) const {

		for( unsigned i=0; i<dim; ++i )
			if(abs(operator()(i) - mr(i)) > std::numeric_limits<Element>::epsilon())
				return false;

		return true;
	} 
	//-----------------------------------------------------------------------
	bool operator != (Self const& mr) const {

		return !(*this == mr);
	}
	//-----------------------------------------------------------------------
	bool operator < (const Element c) const {

		for( unsigned i=0; i<dim; ++i )
			if(operator()(i) > c)
				return false;

		return true;
	}
	//-----------------------------------------------------------------------
	bool operator > (const Element c) const {

		for( unsigned i=0; i<dim; ++i )
			if(operator()(i) < c)
				return false;

		return true;
	}
	//-----------------------------------------------------------------------
	bool abs_lt(const Element c) const {

		assert( c >= Element(0) && "abs_lt() error. negative constant.");

		for( unsigned i=0; i<dim; ++i )
			if(abs(operator()(i)) > c)
				return false;

		return true;
	}
	//-----------------------------------------------------------------------
	bool abs_gt(const Element c) const {

		assert(c >= Element(0) && "abs_gt() error. negative constant."); 
		
		for (unsigned i = 0; i<dim; ++i)
			if (abs(operator()(i)) < c)
				return false;

		return true;
	}
};
//---------------------------------------------------------------------------
//////////////////////////////// Lazy Operations ////////////////////////////
//---------------------------------------------------------------------------
//Addition
//---------------------------------------------------------------------------
template< typename T, typename U>
inline auto operator + ( const MatrixBase<T>& a, const MatrixBase<U>& b ) {

	return Expression<
		const MatrixBase<T>, 
		const MatrixBase<U>, 
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator + ( const MatrixBase<T>& a, const Expression<U, V, O>& b ) {

	return Expression<
		const MatrixBase<T>, 
		const Expression<U, V, O>, 
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator + ( const Expression<U, V, O>& a, const MatrixBase<T>& b ) {

	return Expression<
		const Expression<U, V, O>,
		const MatrixBase<T>,
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template<typename S, typename T, typename O, typename U, typename V, typename P>
inline auto operator + ( const Expression<S, T, O>& a, const Expression<U, V, P>& b ) {

	return Expression<
		const Expression<S, T, O>,
		const Expression<U, V, P>,
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T >
inline auto operator + ( const MatrixBase<T>& a, const typename T::Element b ) {

	return Expression<
		const MatrixBase<T>, 
		const typename T::Element, 
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T >
inline auto operator + ( const typename T::Element a, const MatrixBase<T>& b ) {

	return Expression<
		const typename T::Element, 
		const MatrixBase<T>, 
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator + ( const Expression<U, V, O>& a, const typename Expression<U, V, O>::Element b ) {

	return Expression<
		const Expression<U, V, O>, 
		const typename Expression<U, V, O>::Element, 
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator + ( const typename Expression<U, V, O>::Element a, const Expression<U, V, O>& b ) {

	return Expression<
		const typename Expression<U, V, O>::Element, 
		const Expression<U, V, O>, 
		Sum
	>(a, b);
}
//---------------------------------------------------------------------------
//Subtraction
//---------------------------------------------------------------------------
template< typename T, typename U>
inline auto operator - ( const MatrixBase<T>& a, const MatrixBase<U>& b ) {

	return Expression<
		const MatrixBase<T>, 
		const MatrixBase<U>, 
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator - ( const MatrixBase<T>& a, const Expression<U, V, O>& b ) {

	return Expression<
		const MatrixBase<T>, 
		const Expression<U, V, O>, 
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator - ( const Expression<U, V, O>& a, const MatrixBase<T>& b ) {

	return Expression<
		const Expression<U, V, O>,
		const MatrixBase<T>,
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template<typename S, typename T, typename O, typename U, typename V, typename P>
inline auto operator - ( const Expression<S, T, O>& a, const Expression<U, V, P>& b ) {

	return Expression<
		const Expression<S, T, O>,
		const Expression<U, V, P>,
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T >
inline auto operator - ( const MatrixBase<T>& a, const typename T::Element b ) {

	return Expression<
		const MatrixBase<T>, 
		const typename T::Element, 
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T >
inline auto operator - ( const typename T::Element a, const MatrixBase<T>& b ) {

	return Expression<
		const typename T::Element, 
		const MatrixBase<T>, 
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator - ( const Expression<U, V, O>& a, const typename Expression<U, V, O>::Element b ) {

	return Expression<
		const Expression<U, V, O>, 
		const typename Expression<U, V, O>::Element, 
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator - ( const typename Expression<U, V, O>::Element a, const Expression<U, V, O>& b ) {

	return Expression<
		const typename Expression<U, V, O>::Element, 
		const Expression<U, V, O>, 
		Diff
	>(a, b);
}
//---------------------------------------------------------------------------
//Multiplication
//---------------------------------------------------------------------------
template< int n, typename T >
inline auto operator * (const Vector<n, T>& a, const Vector<n, T>& b) {

	return Expression<
		const Vector<n, T>,
		const Vector<n, T>,
		MulConst
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U>
inline auto operator * ( const MatrixBase<T>& a, const MatrixBase<U>& b ) {

	return Expression<
		const MatrixBase<T>, 
		const MatrixBase<U>, 
		Mul
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator * ( const MatrixBase<T>& a, const Expression<U, V, O>& b ) {

	return Expression<
		const MatrixBase<T>, 
		const Expression<U, V, O>, 
		typename Select2<
			Argument<const MatrixBase<T> >::rows,
			Argument<const Expression<U, V, O> >::rows,
			Mul,
			MulConst
		>::Result
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator * ( const Expression<U, V, O>& a, const MatrixBase<T>& b ) {

	return Expression<
		const Expression<U, V, O>,
		const MatrixBase<T>,
		typename Select2<
			Argument<const Expression<U, V, O> >::rows,
			Argument<const MatrixBase<T> >::rows,
			Mul,
			MulConst
		>::Result
	>(a, b);
}
//---------------------------------------------------------------------------
template<typename S, typename T, typename O, typename U, typename V, typename P>
inline auto operator * ( const Expression<S, T, O>& a, const Expression<U, V, P>& b ) {

	return Expression<
		const Expression<S, T, O>,
		const Expression<U, V, P>,
		typename Select2<
			Argument<const Expression<S, T, O> >::rows,
			Argument<const Expression<U, V, P> >::rows,
			Mul,
			MulConst
		>::Result
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T >
inline auto operator * ( const MatrixBase<T>& a, const typename T::Element b ) {

	return Expression<
		const MatrixBase<T>, 
		const typename T::Element, 
		MulConst
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T >
inline auto operator * ( const typename T::Element a, const MatrixBase<T>& b ) {

	return Expression<
		const typename T::Element, 
		const MatrixBase<T>, 
		MulConst
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator * ( const Expression<U, V, O>& a, const typename Expression<U, V, O>::Element b ) {

	return Expression<
		const Expression<U, V, O>, 
		const typename Expression<U, V, O>::Element, 
		MulConst
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator * ( const typename Expression<U, V, O>::Element a, const Expression<U, V, O>& b ) {

	return Expression<
		const typename Expression<U, V, O>::Element, 
		const Expression<U, V, O>, 
		MulConst
	>(a, b);
}
//---------------------------------------------------------------------------
// Division
//---------------------------------------------------------------------------
template< typename T >
inline auto operator / ( const MatrixBase<T>& a, const typename T::Element b ) {

	assert(b != typename T::Element(0) && "Division dy zero.");

	return Expression<
		const MatrixBase<T>, 
		const typename T::Element, 
		MulConst
	>(a, typename T::Element(1)/b);
}
//---------------------------------------------------------------------------
template< typename U, typename V, typename O >
inline auto operator / ( const Expression<U, V, O>& a, const typename Expression<U, V, O>::Element b ) {

	return Expression<
		const Expression<U, V, O>, 
		const typename Expression<U, V, O>::Element, 
		MulConst
	>(a, static_cast<typename Expression<U, V, O>::Element>(1)/b);
}
//---------------------------------------------------------------------------
///Cross product
//---------------------------------------------------------------------------
template< typename T, typename U >
inline auto operator ^ ( const MatrixBase<T>& a, const MatrixBase<U>& b ) {

	static_assert(
		detail::eq(MatrixBase<T>::dim, MatrixBase<U>::dim) &&
		detail::eq(MatrixBase<T>::dim, 3),
		"Cross product mismatch dimensions."
		);

	return Expression<
		const MatrixBase<T>,
		const MatrixBase<U>,
		Cross
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator ^ ( const MatrixBase<T>& a, const Expression<U, V, O>& b ) {

	static_assert(
		detail::eq(MatrixBase<T>::dim, Expression<U, V, O>::dim) &&
		detail::eq(MatrixBase<T>::dim, 3),
		"Cross product mismatch dimensions."
		);

	return Expression<
		const MatrixBase<T>,
		const Expression<U, V, O>,
		Cross
	>(a, b);
}
//---------------------------------------------------------------------------
template< typename T, typename U, typename V, typename O>
inline auto operator ^ ( const Expression<U, V, O>& a, const MatrixBase<T>& b ) {

	static_assert(
		detail::eq(MatrixBase<T>::dim, Expression<U, V, O>::dim) &&
		detail::eq(MatrixBase<T>::dim, 3),
		"Cross product mismatch dimensions."
		);

	return Expression<
		const Expression<U, V, O>,
		const MatrixBase<T>,
		Cross
	>(a, b);
}
//---------------------------------------------------------------------------
template<typename S, typename T, typename O, typename U, typename V, typename P>
inline auto operator ^ ( const Expression<S, T, O>& a, const Expression<U, V, P>& b ) {

	static_assert(
		detail::eq(Expression<S, T, O>::dim, Expression<U, V, P>::dim) &&
		detail::eq(Expression<S, T, O>::dim, 3),
		"Cross product mismatch dimensions."
		);

	return Expression<
		const Expression<S, T, O>,
		const Expression<U, V, P>,
		Cross
	>(a, b);
}
//---------------------------------------------------------------------------
//Out to stream
//---------------------------------------------------------------------------
template<typename T>
std::ostream& operator << ( std::ostream& s,  const MatrixBase<T>& m) {

	s << '[' << static_cast<size_t>(T::rows)<< 'x' << static_cast<size_t>(T::cols) << ']' << '\n'; 

	for(size_t i = 0; i < T::rows; ++i) {
		for(size_t j = 0; j < T::cols; ++j)
			s << m(i*static_cast<size_t>(T::cols) + j) << " ";

		s << '\n';		
	}

	return s;
}
//---------------------------------------------------------------------------
//Vector construct helper
//---------------------------------------------------------------------------
template<typename... Args>
constexpr auto vec(Args const&... args) {

	return Vector<sizeof...(args), typename std::common_type<Args...>::type>(args...);
}
//---------------------------------------------------------------------------
//Assigment helper function
//---------------------------------------------------------------------------
template<typename T, int r1, int c1, int r2, int c2>
void fit(Matrix<r1,c1,T>& dest, Matrix<r2,c2,T> const& src) {

	memcpy(&dest(0), &src(0), r1*c1*sizeof(T));
}
//---------------------------------------------------------------------------
template<typename T, typename U>
void fit(MatrixBase<T> & dest, MatrixBase<U> const& src) {

	memcpy(&dest(0), &src(0), T::dim*sizeof(T::Element));
}
//---------------------------------------------------------------------------
template<typename T, int n1, int n2>
void fit(Vector<n1,T>& dest, Vector<n2,T> const& src) {

	memcpy(&dest(0), &src(0), n1*sizeof(T));
}
//---------------------------------------------------------------------------
template<typename T>
void fit(Vector<4,T>& dest, Vector<3,T> const& src) {

	dest(0) = src(0);	
	dest(1) = src(1);
	dest(2) = src(2);
	dest(3) = (T)1;
}
//---------------------------------------------------------------------------
// Some math stuff
//---------------------------------------------------------------------------
template<typename T>
MatrixBase<T> min(MatrixBase<T> const& a, MatrixBase<T> const& b) {

	MatrixBase<T> c = a;

	for(size_t i=0; i<a.size(); ++i)
		if( b(i) < c(i) ) c(i) = b(i);

	return c;
}
//---------------------------------------------------------------------------
template<typename T>
MatrixBase<T> max(MatrixBase<T> const& a, MatrixBase<T> const& b) {

	MatrixBase<T> c = a;

	for(size_t i=0; i<a.size(); ++i)
		if( b(i) > c(i) ) c(i) = b(i);

	return c;
}
//-----------------------------------------------------------------------
template<typename T>
inline MatrixBase<T> floor(MatrixBase<T> const& a) {

	MatrixBase<T> res;

	for (unsigned i = 0, s = a.size(); i<s; ++i)
		res(i) = ::floor(a(i));

	return res;
}
//---------------------------------------------------------------------------
template<typename T, size_t n>
Vector<n,T> normalize(Vector<n,T> const& v) {

	auto res = v;
	res.normalize();

	return res;
}
//---------------------------------------------------------------------------
//Matrix inversion
//---------------------------------------------------------------------------
/*
Быстрый метод вычисления обратной матрицы к матрице перемещения 4x4.
Написана и отлажена 12.05.06.
Вход исходная матрица, выход - обращенная матрица, в той же области памяти
http://www.gamedev.ru/faq/?id=74 автор Федор
*/
//---------------------------------------------------------------------------
// translation indices
namespace detail {

	const size_t X = 12;
	const size_t Y = 13;
	const size_t Z = 14;
}
//---------------------------------------------------------------------------
template<typename T>
auto& inverse_affine(Matrix<4, 4, T>& a) {	

	std::swap(a(1, 0), a(0, 1));
	std::swap(a(2, 0), a(0, 2));
	std::swap(a(2, 1), a(1, 2));

	const T v0 = a(0, 0) * a(detail::X) + a(1, 0) * a(detail::Y) + a(2, 0) * a(detail::Z);
	const T v1 = a(0, 1) * a(detail::X) + a(1, 1) * a(detail::Y) + a(2, 1) * a(detail::Z);
	const T v2 = a(0, 2) * a(detail::X) + a(1, 2) * a(detail::Y) + a(2, 2) * a(detail::Z);

	a(detail::X) = -v0;
	a(detail::Y) = -v1;
	a(detail::Z) = -v2;

	return a;
}
//---------------------------------------------------------------------------
template <typename U>
bool inverse_if(MatrixBase<U>& a) {

	///recursive method of matrix inversion by Федор
	typedef typename U::Element T;
	const int n = U::rows;

	static_assert(detail::eq(U::rows, U::cols), "Matrix is not square.");

	Vector<n-1, T> v;

	if (CRecurseInv<n-1>::InvertMatrixRecursionInUnroll(v, a))
		return true;

	return false;
}
//---------------------------------------------------------------------------
template <typename U>
MatrixBase<U>& inverse(MatrixBase<U>& a) {

	///recursive method of matrix inversion by Федор
	typedef typename U::Element T;
	const int n = U::rows;

	static_assert(detail::eq(U::rows, U::cols), "Matrix is not square.");
	Vector<n-1, T> v;

	bool flag = CRecurseInv<n-1>::InvertMatrixRecursionInUnroll(v, a);

	assert(!flag && "Zero determinant.");

	return a;
}
//---------------------------------------------------------------------------
typedef Vector<1,float> vec1f;
typedef Vector<2,float> vec2f;
typedef Vector<3,float> vec3f;
typedef Vector<4,float> vec4f;
typedef Vector<5,float> vec5f;
typedef Vector<6,float> vec6f;
typedef Vector<7,float> vec7f;
typedef Vector<8,float> vec8f;
typedef Vector<9,float> vec9f;
typedef Vector<10,float> vec10f;
typedef Vector<11,float> vec11f;
typedef Vector<12,float> vec12f;
typedef Vector<13,float> vec13f;
typedef Vector<14,float> vec14f;
typedef Vector<15,float> vec15f;
typedef Vector<16,float> vec16f;
typedef Vector<17,float> vec17f;
typedef Vector<18,float> vec18f;
typedef Vector<19,float> vec19f;
typedef Vector<20,float> vec20f;
typedef Vector<21,float> vec21f;
typedef Vector<22,float> vec22f;
typedef Vector<23,float> vec23f;

typedef Vector<1,int> vec1i;
typedef Vector<2,int> vec2i;
typedef Vector<3,int> vec3i;
typedef Vector<4,int> vec4i;

typedef Vector<2,double> vec2d;
typedef Vector<3,double> vec3d;

typedef Vector<3, unsigned char> vec3uc;

typedef Vector<3,unsigned short> vec3us;
typedef Vector<6,unsigned short> vec6us;
typedef Vector<24,unsigned short> vec24us;

typedef Vector<3,unsigned int> vec3ui;
typedef Vector<2, unsigned int> vec2ui;

typedef Matrix<1,2,float> mat1x2f;
typedef Matrix<1,3,float> mat1x3f;
typedef Matrix<1,4,float> mat1x4f;
typedef Matrix<4,1,float> mat4x1f;
typedef Matrix<2,2,float> mat2x2f;
typedef Matrix<3,2,float> mat3x2f;
typedef Matrix<2,3,float> mat2x3f;
typedef Matrix<3,3,float> mat3x3f;
typedef Matrix<3,4,float> mat3x4f;
typedef Matrix<4,3,float> mat4x3f;
typedef Matrix<4,4,float> mat4x4f;
//---------------------------------------------------------------------------
} //< end of mtr namespace
 //---------------------------------------------------------------------------
#ifdef _MSC_VER
#	undef inline
#endif
//---------------------------------------------------------------------------
