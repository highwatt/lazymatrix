//---------------------------------------------------------------------------
/// Written by Pavel Kazmerchuk 2005 - 2006.
/// Rewritten for C++ 11 by Pavel Kazmerchuk 2014.
/// Algorithm of matrices inversion by "Федор" 
/// http://www.gamedev.ru/forum/?group=2&topic=1661&page=3
//---------------------------------------------------------------------------
#pragma once
//---------------------------------------------------------------------------
#include <assert.h>
#include <math.h>
//---------------------------------------------------------------------------
#include <type_traits>
#include <algorithm>
#include <limits>
#include <iosfwd>
//---------------------------------------------------------------------------
#undef abs
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
namespace mtr {
//---------------------------------------------------------------------------
/// Helper classes 
//---------------------------------------------------------------------------
// general unroll
//---------------------------------------------------------------------------
//template <size_t n> struct final{ }; 
////---------------------------------------------------------------------------
//template <size_t i, size_t n, typename Lambda>
//inline void unroll(Lambda const& f,  final<n>) {
//	
//	f(i);
//	unroll<i + 1>(f, final<n - 1>());
//} 
////---------------------------------------------------------------------------
//template <size_t i, typename Lambda>
//inline void unroll(const Lambda& f, final<0>) {
//	f(i);
//}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/// Assigment unroll
//---------------------------------------------------------------------------
template<size_t i, size_t dim>
struct Assigment {

    enum { count = i+1 };

    template<typename T, typename Exp>
    inline static void assign(T& v, Exp const& a) {

        v(i) = a(i);
        Assigment<count, dim>::assign(v, a);
    }
};
//---------------------------------------------------------------------------
template<size_t dim>
struct Assigment<dim, dim> {

    template<typename T, typename Exp>     
    inline static void assign( T& , const Exp&  ){}
};
//---------------------------------------------------------------------------
/// Arguments
//---------------------------------------------------------------------------
template<typename T>
struct Argument {

    enum {
		rows = T::rows,
		cols = T::cols
	};

    typedef typename T::Element Element;

    Argument(T const& arg)
		: arg_(arg) {}

    inline Element operator()(size_t i) const {
        return arg_(i);
    }

	inline bool check_down(const void* p) const { 
		return arg_.check_down(p); 
	}

private:

	T const& arg_;
};
//---------------------------------------------------------------------------
#define SPECIALIZE_ARGUMENT(T)\
template<>\
struct Argument<T> {\
\
	enum {\
		rows = 0,\
		cols = 0\
	};\
\
	typedef T Element;\
\
	Argument(T const& arg)\
		: arg_(arg) {}\
\
	inline T operator()(size_t i) const {\
		return arg_;\
	}\
\
	inline bool check_down(const void* p) const {\
		return true;\
	}\
\
private:\
\
	T const& arg_;\
};
//---------------------------------------------------------------------------
SPECIALIZE_ARGUMENT(double)
SPECIALIZE_ARGUMENT(float)
SPECIALIZE_ARGUMENT(int)
//---------------------------------------------------------------------------
// Expressions
//---------------------------------------------------------------------------
template<class Res, class LeftArg, class RightArg, class Op>
struct Expression {

	typedef Argument<LeftArg>  Left;
	typedef Argument<RightArg> Right;
	typedef Op Operation;
	typedef Res Result;

	enum {
		// if Left::row (Left::col must be checked also) equal zero, then left operand is scalar
		rows = Left::rows ? Left::rows : Right::rows, 
		// if Right::col (Right::row must be checked also) equal zero, then right operand is scalar
		cols = Right::cols ? Right::cols : Left::cols,
		// num elements
		n = rows * cols
	};	

	operator Result() {

		//static_assert( std::is_same<Op, Dot>::value || (rows == cols == 1), "Expression error! Illegal conversion");
		return operator()(0);
	}

    Expression(Left const& left, Right const& right)
        : left_(left)
		, right_(right) 
	{
		static_assert(
			   Left::cols == Right::rows	// mathematical matrix multiplication
			|| Left::rows == 0				// left argument is scalar
			|| Right::rows == 0			// right argument is scalar
			|| (Left::cols == Right::cols && Left::rows == Right::rows && Left::cols == 1) // left argument is the same vector as a right argument
			//|| std::is_same<TOp, Trans>::value // transposition
			, "Expression construction error! Mismatch dimensions of operands.");
	}

	inline Result operator()(size_t i) const {
		return Op::evaluate(i, left_, right_); 
	}

	inline bool check_down(const void* p) const { 
		return left_.check_down(p) && right_.check_down(p) ; 
	}    
	
private:

	const Left  left_;
	const Right right_;
};
//---------------------------------------------------------------------------
///Matrix
//---------------------------------------------------------------------------
template<size_t row, typename T>
struct VectorImpl {

    typedef VectorImpl<row,T> Matrix;    
    typedef T Element;

	enum { 
		rows = row, 
		cols = 1,
		n = rows * cols,  
	};

	void zero() {
		memset(data_, 0, sizeof(data_));
	}
	///indexed access
	inline Element& operator()(size_t i) {
        return data_[i];
    }
	inline Element operator()(size_t i) const {
        return data_[i];
    }
	inline Element& operator[](size_t i) {
        return data_[i];
    }
	inline Element operator[](size_t i) const {
        return data_[i];
    }
	
protected:	

    Element data_[n];	
};
//---------------------------------------------------------------------------
template<size_t row, size_t col, typename T>
struct MatrixImpl {

    typedef MatrixImpl<row,col,T> Matrix;    
    typedef T Element;

	enum { 
		rows = row, 
		cols = col,
		n = row * col,  
	};

	void zero() {
		memset(data_, 0, sizeof(data_));
	}

	void identity() {

		static_assert(row == col, "MatrixImpl::identity() error! Matrix is not square.");
		
		zero(); 

		for(size_t i = 0; i < col; ++i)
			data_[i*col + i] = Element(1);
	}

	inline Element& operator()(size_t i) {
        return data_[i];
    }

	inline Element operator()(size_t i)const {
        return data_[i];
    }

	inline Element& operator[](size_t i) {
        return data_[i];
    }

	inline Element operator[](size_t i)const {
        return data_[i];
    }

	inline Element& operator()(size_t i, size_t j)  {
        return data_[i*col + j];
    }  

	inline Element operator()(size_t i, size_t j)const {
        return data_[i*col + j];
    }
	
protected:	

    Element data_[n];	
};
//---------------------------------------------------------------------------
template<typename M>
struct Matrix : public M {

	typedef typename M::Element Element;
	typedef typename M::Element value_type;

	enum { 
		rows = M::rows, 
		cols = M::cols,
		n    = M::n  
	};

	inline const Element* begin() const {
        return &data_[0];
    }
	inline const Element* end()const {
        return &data_[n];
    }
	inline size_t size()const {
		return n;
	}
    inline bool check_down(const void* p) const { 
		return begin() != p; 
	}
	//---------------------------------------------------------------------------	
	template<class R, class T, class U, class O>
    inline Matrix& operator = ( const Expression<R,T,U,O>& expression ) {

		typedef Expression<R,T,U,O> E;	

		static_assert( 
			(cols == E::cols && rows == E::rows) || 
			(std::is_same<O,Trans>::value && cols == E::rows && rows == E::cols), 
			"MatrixBase::operator=() error! Mismatch dimensions."
		);

		assert( expression.check_down(begin()) );
	
		Assigment<0,dim>::assign(*this, expression);

        return *this;
    }

    inline Matrix& operator = ( const Matrix& m ) {	

		memcpy(&data_[0], m.begin(), m.size()*sizeof(Element));
        return *this;
    }
};
//---------------------------------------------------------------------------
template< typename TLeft, typename TRight>
inline Expression< const TLeft, const TRight, SELECT_OP(Mul, MulConst) >
operator * ( const TLeft& a, const TRight& b )
{
	typedef typename Select
		< 
			Argument<const TLeft>::rows, 
			const TLeft, 
			typename Argument<const TLeft>::Element
		>::Result U; 

	typedef typename Select
		< 
			Argument<const TRight>::rows, 
			const TRight, 
			typename Argument<const TRight>::Element
		>::Result V;
	
	return Expression<const U, const V, SELECT_OP(Mul, MulConst)>(a, b);
}
//---------------------------------------------------------------------------
template< typename TLeft, typename TRight>
inline Expression< const TLeft, const TRight, Sum >
operator + ( const TLeft& a, const TRight& b )
{
    return Expression<const TLeft, const TRight, Sum>(a, b);
}
//---------------------------------------------------------------------------
template< typename TLeft, typename TRight>
inline Expression< const TLeft, const TRight, Diff >
operator - ( const TLeft& a, const TRight& b )
{	
	return Expression<const TLeft, const TRight, Diff>(a, b);
}
//---------------------------------------------------------------------------
///Cross product
template< typename TLeft, typename TRight>
inline Expression< const TLeft, const TRight, Cross >
operator ^ ( const TLeft& a, const TRight& b )
{
    STATIC_CHECK( TLeft::dim == TRight::dim && TLeft::dim == 3 , Mismatch_dimensions );
	
	return Expression<const TLeft, const TRight, Cross>(a, b);
}
//---------------------------------------------------------------------------
template<typename T>
std::ostream& operator << ( std::ostream& s,  const MatrixBase<T>& m)
{
	s << '[' << T::rows << 'x' << T::cols << ']'; s.put(s.widen('\n'));
	
	for(size_t i = 0; i < T::rows; ++i)
	{
		for(size_t j = 0; j < T::cols; ++j)
		 	s << m(i*(size_t)T::cols + j) << " ";
				
		s.put(s.widen('\n'));		
	}

	s.flush();

return s;
}
//---------------------------------------------------------------------------
template<typename T, int r1, int c1, int r2, int c2>
void make(MatrixBase<MatrixImpl<r1,c1,T> >& mtr1, MatrixBase<MatrixImpl<r2,c2,T> > const& mtr2)
{
	for(size_t i = 0;i<r1*c1; ++i)
		mtr1(i) = mtr2(i);
}
//---------------------------------------------------------------------------
template<typename T, typename U>
void make(MatrixBase<T> & mtr1, MatrixBase<U> const& mtr2)
{
	for(size_t i = 0;i<T::dim; ++i)
		mtr1(i) = mtr2(i);
}
//---------------------------------------------------------------------------
template<typename T, int n1, int n2>
void make(MatrixBase<VectorImpl<n1,T> >& mtr1, MatrixBase<VectorImpl<n2,T> > const& mtr2)
{
	for(size_t i = 0;i<n1; ++i)
		mtr1(i) = mtr2(i);
}
//---------------------------------------------------------------------------
template<typename T>
MatrixBase<T>& inverseAffine4x4(MatrixBase<T>& a)
{
	//Быстрый метод вычисление обратной матрицы к матрице перемещения 4x4
	//Написана и отлажена 12.05.06
	//вход исходная матрица, выход - обращенная матрица, в той же области памяти
	
	typedef MatrixBase<T> M;

	STATIC_CHECK( (M::rows == M::cols) && (M::cols == 4), Mismatch_dimensions);	

	//Транспозиция левой верхней 3x3 матрицы
	std::swap (a(1,0), a(0,1));
	std::swap (a(2,0), a(0,2)); 
	std::swap (a(2,1), a(1,2));

	//Вычисления нижнего вектора-строки
	const M::Element v0 = a(0,0) * a(3,0) + a(1,0) * a(3,1) + a(2,0) * a(3,2);
	const M::Element v1 = a(0,1) * a(3,0) + a(1,1) * a(3,1) + a(2,1) * a(3,2);
	const M::Element v2 = a(0,2) * a(3,0) + a(1,2) * a(3,1) + a(2,2) * a(3,2);

	a(3,0) = -v0;
	a(3,1) = -v1;
	a(3,2) = -v2;

return a;
}
//---------------------------------------------------------------------------
template <typename U>
bool inverseif(MatrixBase<U>& a)
{
	///recursive method of matrix inversion by Федор
	typedef U::Element T;
	const int nDim = U::rows;

	STATIC_CHECK(U::rows == U::cols, Matrix_is_not_square);
	
	mtr::MatrixBase< mtr::VectorImpl<nDim - 1,T> > v; 
	  
	if (CRecurseInv<nDim-1>::InvertMatrixRecursionInUnroll(v, a))    
		return true;	
	
	return false;
}
//---------------------------------------------------------------------------
template <typename U>
MatrixBase<U>& inverse(MatrixBase<U>& a)
{
	///recursive method of matrix inversion by Федор
	typedef U::Element T;
	const int nDim = U::rows;

	STATIC_CHECK(U::rows == U::cols, Matrix_is_not_square);
	mtr::MatrixBase< mtr::VectorImpl<nDim - 1,T> > v; 

	const bool flag = CRecurseInv<nDim-1>::InvertMatrixRecursionInUnroll(v, a);
	  
	assert(!flag && "Zero determinant");
	
	return a;
}
//---------------------------------------------------------------------------
} // mtr namespace
//---------------------------------------------------------------------------
#define TYPEDEF_MATRIX(r,c,type) typedef mtr::MatrixBase< mtr::MatrixImpl<r,c,type> >
#define TYPEDEF_VECTOR(n,type) typedef mtr::MatrixBase< mtr::VectorImpl<n,type> >

TYPEDEF_VECTOR(1,float) Vector1f;
TYPEDEF_VECTOR(2,float) Vector2f;
TYPEDEF_VECTOR(3,float) Vector3f;
TYPEDEF_VECTOR(4,float) Vector4f;
TYPEDEF_VECTOR(5,float) Vector5f;
TYPEDEF_VECTOR(6,float) Vector6f;
TYPEDEF_VECTOR(7,float) Vector7f;
TYPEDEF_VECTOR(8,float) Vector8f;
TYPEDEF_VECTOR(9,float) Vector9f;
TYPEDEF_VECTOR(10,float) Vector10f;
TYPEDEF_VECTOR(11,float) Vector11f;
TYPEDEF_VECTOR(12,float) Vector12f;
TYPEDEF_VECTOR(13,float) Vector13f;
TYPEDEF_VECTOR(14,float) Vector14f;
TYPEDEF_VECTOR(15,float) Vector15f;
TYPEDEF_VECTOR(16,float) Vector16f;
TYPEDEF_VECTOR(17,float) Vector17f;
TYPEDEF_VECTOR(18,float) Vector18f;
TYPEDEF_VECTOR(19,float) Vector19f;
TYPEDEF_VECTOR(20,float) Vector20f;
TYPEDEF_VECTOR(21,float) Vector21f;
TYPEDEF_VECTOR(22,float) Vector22f;
TYPEDEF_VECTOR(23,float) Vector23f;

TYPEDEF_VECTOR(1,int) Vector1i;
TYPEDEF_VECTOR(2,int) Vector2i;
TYPEDEF_VECTOR(3,int) Vector3i;
TYPEDEF_VECTOR(4,int) Vector4i;

TYPEDEF_VECTOR(3,unsigned short) Vector3us;
TYPEDEF_VECTOR(3,unsigned int) Vector3ui;

TYPEDEF_MATRIX(1,2,float) Matrix1x2f;
TYPEDEF_MATRIX(1,3,float) Matrix1x3f;
TYPEDEF_MATRIX(1,4,float) Matrix1x4f;
TYPEDEF_MATRIX(4,1,float) Matrix4x1f;
TYPEDEF_MATRIX(2,2,float) Matrix2x2f;
TYPEDEF_MATRIX(3,2,float) Matrix3x2f;
TYPEDEF_MATRIX(2,3,float) Matrix2x3f;
TYPEDEF_MATRIX(3,3,float) Matrix3x3f;
TYPEDEF_MATRIX(3,4,float) Matrix3x4f;
TYPEDEF_MATRIX(4,3,float) Matrix4x3f;
TYPEDEF_MATRIX(4,4,float) Matrix4x4f;
//---------------------------------------------------------------------------
#ifdef _MSC_VER
#	undef inline
#endif