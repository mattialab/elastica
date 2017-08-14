#ifndef MATRIX3_H
#define MATRIX3_H

#include "UsualHeaders.h"
#include "Vector3.h"
#include "Tolerance.h"

using namespace std;

class Matrix3
{
protected:
	typedef std::vector<Vector3> VV3;
	typedef std::vector<REAL> VREAL;
	typedef std::vector<bool> VBOOL;
	typedef std::vector<int> VINT;

	static const unsigned int PRINTPRECISION = 6;

public:
	// This shouldnt be here, super dangerous, wtf!!
	REAL	r1c1, r1c2, r1c3,
			r2c1, r2c2, r2c3,
			r3c1, r3c2, r3c3;

    // constructors
    Matrix3();
    Matrix3(REAL r1c1, REAL r1c2, REAL r1c3, REAL r2c1, REAL r2c2, REAL r2c3, REAL r3c1, REAL r3c2, REAL r3c3);
    Matrix3(const Vector3& S); // skew matrix <-> vector mapping
    Matrix3(const Vector3& r1, const Vector3& r2, const Vector3& r3);
    Matrix3(const Matrix3& M);
    
    // matrix operators
    inline Matrix3 operator+(const Matrix3& M) const
    {
    	return Matrix3(	r1c1 + M.r1c1, r1c2 + M.r1c2, r1c3 + M.r1c3,
    					r2c1 + M.r2c1, r2c2 + M.r2c2, r2c3 + M.r2c3,
    					r3c1 + M.r3c1, r3c2 + M.r3c2, r3c3 + M.r3c3);
    }

    inline Matrix3 operator-(const Matrix3& M) const
    {
    	return Matrix3(	r1c1 - M.r1c1, r1c2 - M.r1c2, r1c3 - M.r1c3,
    					r2c1 - M.r2c1, r2c2 - M.r2c2, r2c3 - M.r2c3,
    					r3c1 - M.r3c1, r3c2 - M.r3c2, r3c3 - M.r3c3);
    }

    inline Matrix3 operator*(const Matrix3& M) const // matrix multiplication
    {
    	return Matrix3(	r1c1 * M.r1c1 + r1c2 * M.r2c1 + r1c3 * M.r3c1,
    					r1c1 * M.r1c2 + r1c2 * M.r2c2 + r1c3 * M.r3c2,
    					r1c1 * M.r1c3 + r1c2 * M.r2c3 + r1c3 * M.r3c3,
    					r2c1 * M.r1c1 + r2c2 * M.r2c1 + r2c3 * M.r3c1,
    					r2c1 * M.r1c2 + r2c2 * M.r2c2 + r2c3 * M.r3c2,
    					r2c1 * M.r1c3 + r2c2 * M.r2c3 + r2c3 * M.r3c3,
    					r3c1 * M.r1c1 + r3c2 * M.r2c1 + r3c3 * M.r3c1,
    					r3c1 * M.r1c2 + r3c2 * M.r2c2 + r3c3 * M.r3c2,
    					r3c1 * M.r1c3 + r3c2 * M.r2c3 + r3c3 * M.r3c3);
    }

    // vector operators
    inline Vector3 operator*(const Vector3& v) const
    {
    	return Vector3(	r1c1 * v.x + r1c2 * v.y + r1c3 * v.z,
    					r2c1 * v.x + r2c2 * v.y + r2c3 * v.z,
    					r3c1 * v.x + r3c2 * v.y + r3c3 * v.z);
    }

    // scalar operators
    inline Matrix3 operator*(const REAL a) const
    {
    	return Matrix3(	a * r1c1, a * r1c2, a * r1c3,
    					a * r2c1, a * r2c2, a * r2c3,
    					a * r3c1, a * r3c2, a * r3c3);
    }

    inline Matrix3 operator/(const REAL a) const
    {
    	return Matrix3(	r1c1 / a, r1c2 / a, r1c3 / a,
    					r2c1 / a, r2c2 / a, r2c3 / a,
    					r3c1 / a, r3c2 / a, r3c3 / a);
    }

    inline Matrix3 operator-() const
    {
    	return Matrix3(	-r1c1, -r1c2, -r1c3,
    					-r2c1, -r2c2, -r2c3,
    					-r3c1, -r3c2, -r3c3);
    }

    // in-place modifiers
    inline Matrix3& operator+=(const Matrix3& M)
    {
    	r1c1 += M.r1c1; r1c2 += M.r1c2; r1c3 += M.r1c3;
    	r2c1 += M.r2c1; r2c2 += M.r2c2; r2c3 += M.r2c3;
    	r3c1 += M.r3c1; r3c2 += M.r3c2; r3c3 += M.r3c3;
    	return *this;
    }

    inline Matrix3& operator-=(const Matrix3& M)
    {
    	r1c1 -= M.r1c1; r1c2 -= M.r1c2; r1c3 -= M.r1c3;
    	r2c1 -= M.r2c1; r2c2 -= M.r2c2; r2c3 -= M.r2c3;
    	r3c1 -= M.r3c1; r3c2 -= M.r3c2; r3c3 -= M.r3c3;
    	return *this;
    }

    inline Matrix3& operator*=(const Matrix3& M)
    {
    	*this = *this * M;
    	return *this;
    }

    inline Matrix3& operator*=(const REAL a)
    {
    	r1c1 *= a; r1c2 *= a; r1c3 *= a;
    	r2c1 *= a; r2c2 *= a; r2c3 *= a;
    	r3c1 *= a; r3c2 *= a; r3c3 *= a;
    	return *this;
    }

    inline Matrix3& operator/=(const REAL a)
    {
    	r1c1 /= a; r1c2 /= a; r1c3 /= a;
    	r2c1 /= a; r2c2 /= a; r2c3 /= a;
    	r3c1 /= a; r3c2 /= a; r3c3 /= a;
    	return *this;
    }

    // general functions
    inline Matrix3 T() const // transpose
    {
    	return Matrix3(	r1c1, r2c1, r3c1,
    					r1c2, r2c2, r3c2,
    					r1c3, r2c3, r3c3);
    }

    inline void T(Matrix3& A) const // transpose
    {
    	A.r1c1 = r1c1;
		A.r1c2 = r2c1;
		A.r1c3 = r3c1;
		A.r2c1 = r1c2;
		A.r2c2 = r2c2;
		A.r2c3 = r3c2;
		A.r3c1 = r1c3;
		A.r3c2 = r2c3;
		A.r3c3 = r3c3;
    }

    inline Matrix3 I() const // inverse
    {
    	const REAL D = det();
    	assert(D>=Tolerance::tol());

    	const REAL _r1c1 = (r2c2 * r3c3 - r3c2 * r2c3) / D;
    	const REAL _r1c2 = (r3c2 * r1c3 - r1c2 * r3c3) / D;
    	const REAL _r1c3 = (r1c2 * r2c3 - r2c2 * r1c3) / D;

    	const REAL _r2c1 = (r3c1 * r2c3 - r2c1 * r3c3) / D;
    	const REAL _r2c2 = (r1c1 * r3c3 - r3c1 * r1c3) / D;
    	const REAL _r2c3 = (r2c1 * r1c3 - r1c1 * r2c3) / D;

    	const REAL _r3c1 = (r2c1 * r3c2 - r3c1 * r2c2) / D;
    	const REAL _r3c2 = (r3c1 * r1c2 - r1c1 * r3c2) / D;
    	const REAL _r3c3 = (r1c1 * r2c2 - r2c1 * r1c2) / D;

    	return Matrix3(	_r1c1, _r1c2 , _r1c3,
    					_r2c1, _r2c2 , _r2c3,
    					_r3c1, _r3c2 , _r3c3);
    }

    inline void diagI(Matrix3& b) const // inverse of a diagonal matrix in-place
    {
    	assert(r2c1==0.0);
    	assert(r3c1==0.0);
    	assert(r1c2==0.0);
    	assert(r3c2==0.0);
    	assert(r1c3==0.0);
    	assert(r2c3==0.0);
    	assert(r1c1>=Tolerance::tol());
    	assert(r2c2>=Tolerance::tol());
    	assert(r3c3>=Tolerance::tol());

    	assert(b.r2c1==0.0);
    	assert(b.r3c1==0.0);
    	assert(b.r1c2==0.0);
    	assert(b.r3c2==0.0);
    	assert(b.r1c3==0.0);
    	assert(b.r2c3==0.0);

    	b.r1c1 = 1.0/r1c1;
    	b.r2c2 = 1.0/r2c2;
    	b.r3c3 = 1.0/r3c3;
    }

    inline REAL tr() const // trace
    {
    	return r1c1 + r2c2 + r3c3;
    }

    inline REAL det() const // determinant
    {
    	const REAL D = (+r1c1  *  (r2c2 * r3c3 - r2c3 * r3c2)
    					-r1c2  *  (r2c1 * r3c3 - r2c3 * r3c1)
    					+r1c3  *  (r2c1 * r3c2 - r2c2 * r3c1));

    	return D;
    }

    inline Vector3 log(const REAL tol = Tolerance::tol()) const // matrix logarithm (\in so(3)) represented as vector
    {
    	const REAL theta = arcCos((tr() - 1.0) / 2.0);
    	return (theta < tol)?Vector3():( theta/(2.0*sin(theta))*Vector3(r2c3-r3c2, r3c1-r1c3, r1c2-r2c1) ); // if theta ~= 0 return zero vector
    }

    inline void log(Vector3& A, const REAL tol = Tolerance::tol()) const // matrix logarithm (\in so(3)) represented as vector
    {
    	// Clean vector
    	A.x = 0.0;
    	A.y = 0.0;
    	A.z = 0.0;

    	const REAL theta = arcCos((tr() - 1.0) / 2.0);

        const REAL factor = (theta<=tol)?0.0:(theta/(2.0*sin(theta)));
        A.x = factor*(r2c3-r3c2);
        A.y = factor*(r3c1-r1c3);
        A.z = factor*(r1c2-r2c1);
    }

    // norms
    inline REAL F() const // frobenius
    {
    	return sqrt(r1c1 * r1c1 + r1c2 * r1c2 + r1c3 * r1c3 +
    				r2c1 * r2c1 + r2c2 * r2c2 + r2c3 * r2c3 +
    				r3c1 * r3c1 + r3c2 * r3c2 + r3c3 * r3c3);
    }

    // utility functions
    inline Vector3 S() const // vector associated with skew matrix
    {
    	return Vector3(r2c3, r3c1, r1c2);
    }

    inline Vector3 operator[](int rowNum) const // returns the ith row
    {
    	switch (rowNum)
    	{
    	case 0 : return Vector3(r1c1, r1c2, r1c3);
    	case 1 : return Vector3(r2c1, r2c2, r2c3);
    	case 2 : return Vector3(r3c1, r3c2, r3c3);
    	default: throw "Bad [] access attemped on Matrix3.\n";
    	}
    }

    inline bool isOrthogonal(const REAL tol = Tolerance::tol()) const // returns true if orthogonal within tol
    {
    	Matrix3 QQT = (*this) * (*this).T();
    	if ((QQT - Matrix3()).F() < tol)
    		return true;
    	else
    		return false;
    }

    inline bool goodNumerics() const
    {
    	if ( (r1c1!=r1c1) || (r1c2!=r1c2) || (r1c3!=r1c3) || (r2c1!=r2c1) || (r2c2!=r2c2) || (r2c3!=r2c3) || (r3c1!=r3c1) || (r3c2!=r3c2) || (r3c3!=r3c3) )
    		return false;
    	else
    		return true;
    }

    // speed functions
    static inline void matrix_times_vector(const Matrix3& a, const Vector3& b, Vector3& c)
    {
    	c.x = a.r1c1 * b.x + a.r1c2 * b.y + a.r1c3 * b.z;
    	c.y = a.r2c1 * b.x + a.r2c2 * b.y + a.r2c3 * b.z;
    	c.z = a.r3c1 * b.x + a.r3c2 * b.y + a.r3c3 * b.z;
    }

    static inline void matrix_times_transposed(const Matrix3& a, const Matrix3& b, Matrix3& c)
    {
    	c.r1c1 = a.r1c1 * b.r1c1 + a.r1c2 * b.r1c2 + a.r1c3 * b.r1c3;
    	c.r1c2 = a.r1c1 * b.r2c1 + a.r1c2 * b.r2c2 + a.r1c3 * b.r2c3;
    	c.r1c3 = a.r1c1 * b.r3c1 + a.r1c2 * b.r3c2 + a.r1c3 * b.r3c3;

    	c.r2c1 = a.r2c1 * b.r1c1 + a.r2c2 * b.r1c2 + a.r2c3 * b.r1c3;
    	c.r2c2 = a.r2c1 * b.r2c1 + a.r2c2 * b.r2c2 + a.r2c3 * b.r2c3;
    	c.r2c3 = a.r2c1 * b.r3c1 + a.r2c2 * b.r3c2 + a.r2c3 * b.r3c3;

    	c.r3c1 = a.r3c1 * b.r1c1 + a.r3c2 * b.r1c2 + a.r3c3 * b.r1c3;
    	c.r3c2 = a.r3c1 * b.r2c1 + a.r3c2 * b.r2c2 + a.r3c3 * b.r2c3;
    	c.r3c3 = a.r3c1 * b.r3c1 + a.r3c2 * b.r3c2 + a.r3c3 * b.r3c3;
    }

    // I/O
    friend std::ostream& operator<< (std::ostream& out, const Matrix3& A);
};


////////////////// operators not declarable inside Matrix3 //////////

Matrix3 operator*(REAL a, const Matrix3& M);

// treating Vector3's as row vectors here: we need to be careful in using it!
Vector3 operator*(const Vector3& v, const Matrix3& M);

// exponential of a vector (\in so(3)) maps to an orthogonal matrix (\in SO(3)) using rodiguez formula.
// The size (angle) of the rotation is just the magnitude of the vector.
Matrix3 exp(const Vector3& S, const REAL tol = Tolerance::tol());
void exp(const Vector3& S, Matrix3& B);

// generates a uniform gaussian random variable in each element
Matrix3 randMatrix3();

#endif


