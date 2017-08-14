#ifndef VECTORFUNCTIONS_H_
#define VECTORFUNCTIONS_H_

#include <typeinfo>
#include <algorithm>
#include "UsualHeaders.h"
#include "ArithmeticPrecision.h"
#include "Vector3.h"
#include "Matrix3.h"

// All templated functions cannot be split in .h/.cpp, so this file is properly done in this sense

using namespace std;

namespace VEF
{
typedef std::vector<Vector3> VV3;
typedef std::vector<Matrix3> VM3;
typedef std::vector<REAL> VREAL;
typedef std::vector<bool> VBOOL;
typedef std::vector<int> VINT;
}

/////////////////////////////// Binary Vector operators ////////////////////////////////

template <typename T, typename U>
inline vector<decltype(T() * U())> operator*(const vector<T>& A, const vector<U>& B)
{

#ifndef NDEBUG
	const unsigned int n = A.size();
	if (n  !=  B.size())
	{
		ostringstream oss (ostringstream::out);
		oss << "Size mismatch in operator*(" << typeid(A).name() << "," <<  typeid(B).name()
					<< ", sizes: (" << (int) A.size() << ',' << (int) B.size() << "), crashing!";
		throw(oss.str());
	}
#endif

	vector<decltype(T() * U())> C;
	C.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	typename vector<U>::const_iterator b = B.begin();
	while (a != A.end())
	{
		C.push_back((*a) * (*b));
		++a;
		++b;
	}
	return C;
}

template <typename T, typename U>
inline vector<decltype(T() + U())> operator+(const vector<T>& A, const vector<U>& B)
{

#ifndef NDEBUG
	const unsigned int n = A.size();
	if (n  !=  B.size())
	{
		ostringstream oss (ostringstream::out);
		oss << "Size mismatch in operator+(" << typeid(A).name() << "," <<  typeid(B).name()
					<< ", sizes: (" << (int) A.size() << ',' << (int) B.size() << "), crashing!";
		throw(oss.str());
	}
#endif

	vector<decltype(T() + U())> C;
	C.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	typename vector<U>::const_iterator b = B.begin();
	while (a != A.end())
	{
		C.push_back((*a) + (*b));
		++a;
		++b;
	}
	return C;
}

template <typename T, typename U>
inline vector<decltype(T() - U())> operator-(const vector<T>& A, const vector<U>& B)
{

#ifndef NDEBUG
	const unsigned int n = A.size();
	if (n  !=  B.size())
	{
		ostringstream oss (ostringstream::out);
		oss << "Size mismatch in operator-(" << typeid(A).name() << "," <<  typeid(B).name()
					<< ", sizes: (" << (int) A.size() << ',' << (int) B.size() << "), crashing!";
		throw(oss.str());
	}
#endif

	vector<decltype(T() - U())> C;
	C.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	typename vector<U>::const_iterator b = B.begin();
	while (a != A.end())
	{
		C.push_back((*a) - (*b));
		++a;
		++b;
	}
	return C;
}

template <typename T, typename U>
inline vector<decltype(T() / U())> operator/(const vector<T>& A, const vector<U>& B){
#ifndef NDEBUG
	unsigned int n = A.size();
	if (n  !=  B.size()){
		ostringstream oss (ostringstream::out);
		oss << "Size mismatch in operator/(" << typeid(A).name() << "," <<  typeid(B).name()
			<< ", sizes: (" << (int) A.size() << ',' << (int) B.size() << "), crashing!";
		throw(oss.str());
	}
#endif
	vector<decltype(T() / U())> C;
	C.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	typename vector<U>::const_iterator b = B.begin();
	while (a != A.end()){
		C.push_back((*a) / (*b));
		++a; ++b;
	}
	return C;
}

template <typename T, typename U>
inline vector<decltype(T() % U())> operator%(const vector<T>& A, const vector<U>& B){
#ifndef NDEBUG
	unsigned int n = A.size();
	if (n  !=  B.size()){
		ostringstream oss (ostringstream::out);
		oss << "Size mismatch in operator%(" << typeid(A).name() << "," <<  typeid(B).name()
			<< ", sizes: (" << (int) A.size() << ',' << (int) B.size() << "), crashing!";
		throw(oss.str());
	}
#endif
	vector<decltype(T() % U())> C;
	C.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	typename vector<U>::const_iterator b = B.begin();
	while (a != A.end()){
		C.push_back((*a) % (*b));
		++a; ++b;
	}
	return C;
}

///////////////////// Vector-Object operators //////////////////////////////

template <typename T, typename U>
inline vector<decltype(T()*U())> operator+(const vector<T>& A, const U t){
	vector<decltype(T()*U())> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a) + t);
		++a;}
	return B;
}

template <typename T, typename U>
inline vector<decltype(T()*U())> operator-(const vector<T>& A, const U t){
	vector<decltype(T()*U())> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a) - t);
		++a;}
	return B;
}

template <typename T, typename U>
inline vector<decltype(T()*U())> operator*(const vector<T>& A, const U t){
	vector<decltype(T()*U())> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a) * t);
		++a;}
	return B;
}

template <typename T, typename U>
inline vector<decltype(T()/U())> operator/(const vector<T>& A, const U t){
	vector<decltype(T()/U())> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a) / t);
		++a;}
	return B;
}

template <typename T, typename U>
inline void operator+=(vector<T>& A, const U t){
	typename vector<T>::iterator a = A.begin();
	while (a != A.end()){
		*a += t;
		++a;
	}
}

template <typename T, typename U>
inline void operator-=(vector<T>& A, const U t){
	typename vector<T>::iterator a = A.begin();
	while (a != A.end()){
		*a -= t;
		++a;
	}
}

template <typename T, typename U>
inline void operator*=(vector<T>& A, const U t){
	typename vector<T>::iterator a = A.begin();
	while (a != A.end()){
		*a *= t;
		++a;
	}
}

template <typename T, typename U>
inline void operator/=(vector<T>& A, const U t){
	typename vector<T>::iterator a = A.begin();
	while (a != A.end()){
		*a /= t;
		++a;
	}
}

///////////////////// Object-Vector operators ///////////////////////////////


template <typename T, typename U>
inline vector<decltype(T()*U())> operator*(const T& t, const vector<U>& A){
	vector<decltype(T()*U())> B;
	B.reserve(A.size());
	typename vector<U>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back(t * (*a));
		++a;}
	return B;
}

template <typename T, typename U>
inline vector<decltype(T()/U())> operator/(const T& t, const vector<U>& A){
	vector<decltype(T()/U())> B;
	B.reserve(A.size());
	typename vector<U>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back(t / (*a));
		++a;}
	return B;
}


///////////////////// Vector-Scalar operators ///////////////////////////////

template <typename T>
inline vector<T> operator*(const vector<T>& A, REAL r){
	vector<T> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a)*r);
		++a;}
	return B;
}

template <typename T>
inline vector<T> operator*(REAL r, const vector<T>& A){
	vector<T> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back(r*(*a));
		++a;}
	return B;
}

template <typename T>
inline vector<T> operator/(const vector<T>& A, REAL r){
	vector<T> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a)/r);
		++a;}
	return B;
}

///////////////////// Unary Vector operators //////////////////////////////

template <typename T>
inline vector<T> operator-(const vector<T>& A){
	vector<T> B;
	B.reserve(A.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back(-(*a));
		++a;}
	return B;
}

////////////////// In-Place operators ////////////////////////////////////

template <typename T>
inline void operator+=(vector<T>& A, const vector<T>& B){
	typename vector<T>::iterator a = A.begin();
	typename vector<T>::const_iterator b = B.begin();
	while (a != A.end()){
		(*a) += (*b);
		++a; ++b;}
}

template <typename T>
inline void operator-=(vector<T>& A, const vector<T>& B){
	typename vector<T>::iterator a = A.begin();
	typename vector<T>::const_iterator b = B.begin();
	while (a != A.end()){
		(*a) -= (*b);
		++a; ++b;}
}

template <typename T>
inline void operator*=(vector<T>& A, const vector<T>& B){
	typename vector<T>::iterator a = A.begin();
	typename vector<T>::const_iterator b = B.begin();
	while (a != A.end()){
		(*a) *= (*b);
		++a; ++b;}
}

template <typename T>
inline void operator/=(vector<T>& A, const vector<T>& B){
	typename vector<T>::iterator a = A.begin();
	typename vector<T>::const_iterator b = B.begin();
	while (a != A.end()){
		(*a) /= (*b);
		++a; ++b;}
}

template <typename T>
inline void operator*=(vector<T>& A, REAL r){
	typename vector<T>::iterator a = A.begin();
	while (a != A.end()){
		(*a) *= r;
		++a;}
}

template <typename T>
inline void operator/=(vector<T>& A, REAL r){
	typename vector<T>::iterator a = A.begin();
	while (a != A.end()){
		(*a) /= r;
		++a;}
}

///////////////////////  Input-Output ///////////////////////////////////

template <typename T>
inline ostream& operator<<(ostream& out, const vector<T>& A){
	typename vector<T>::const_iterator a = A.begin();
	if (A.size() > 0){
		while(a != A.end()-1){
			out << *a << endl;
			++a;
		}
		out << *a;
	}
	return out;
}

////////////////////// VM3 calls //////////////////////////


inline VEF::VM3 vT(const VEF::VM3& A){
	VEF::VM3 B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).T());
		++a;}
	return B;
}

inline VEF::VM3 vI(const VEF::VM3& A){
	VEF::VM3 B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).I());
		++a;}
	return B;
}

inline VEF::VREAL vTr(const VEF::VM3& A){
	VEF::VREAL B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).tr());
		++a;}
	return B;
}

inline VEF::VREAL vDet(const VEF::VM3& A){
	VEF::VREAL B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).det());
		++a;}
	return B;
}

inline VEF::VV3 vLog(const VEF::VM3& A){
	VEF::VV3 B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).log());
		++a;
	}
	return B;
}

inline VEF::VREAL vF(const VEF::VM3& A){
	VEF::VREAL B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).F());
		++a;}
	return B;
}

inline VEF::VV3 vS(const VEF::VM3& A){
	VEF::VV3 B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).S());
		++a;
	}
	return B;
}

inline VEF::VM3 vExp(const VEF::VV3& A, REAL tol = Tolerance::tol()){
	VEF::VM3 B;
	B.reserve(A.size());
	VEF::VV3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back(exp(*a));
		++a;
	}
	return B;
}

inline VEF::VM3 vRandMatrix3(int n){
	VEF::VM3 A(n);
	generate(A.begin(), A.end(), randMatrix3);
	return A;
}

inline VEF::VBOOL vIsOrthogonal(const VEF::VM3& A, REAL tol = Tolerance::tol()){
	VEF::VBOOL B;
	B.reserve(A.size());
	VEF::VM3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).isOrthogonal());
		++a;
	}
	return B;
}

inline VEF::VM3 vRotDiff(const VEF::VM3& A){
	unsigned int n = A.size();
	VEF::VM3 B;
	B.reserve(n-1);
	for (unsigned int i=1; i<n; i++)
		B.push_back(A[i]*(A[i-1].T()));
	return B;
}

//////////////////////////// VV3 calls //////////////////////////

inline VEF::VV3 vUnitize(const VEF::VV3& A){
	VEF::VV3 B;
	B.reserve(A.size());
	VEF::VV3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).unitize());
		++a;
	}
	return B;
}

inline VEF::VREAL vLength(const VEF::VV3& A){
	VEF::VREAL B;
	B.reserve(A.size());
	VEF::VV3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).length());
		++a;
	}
	return B;
}

inline VEF::VBOOL vIsUnit(const VEF::VV3& A, REAL tol = Tolerance::tol()){
	VEF::VBOOL B;
	B.reserve(A.size());
	VEF::VV3::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back((*a).isUnit());
		++a;
	}
	return B;
}

inline VEF::VV3 vRandVector3(int n, REAL mu, REAL sigma){
	VEF::VV3 A(n);
	for (int i=0; i<n; i++)
		A[i] = randVector3(mu, sigma);
	return A;
}

inline VEF::VV3 vRandVector3(int n){
	return vRandVector3(n, 0, 1);
	//VV3 A(n);
	//generate(A.begin(), A.end(), randVector3);
	//return A;
}

inline VEF::VV3 vRandXOnly(int n){
	VEF::VV3 A(n);
	generate(A.begin(), A.end(), randXOnly);
	return A;
}



/////////////////////// REAL calls /////////////////////////////////

inline VEF::VREAL vSqrt(const VEF::VREAL& A){
	VEF::VREAL B;
	B.reserve(A.size());
	VEF::VREAL::const_iterator a = A.begin();
	while (a != A.end()){
		B.push_back(sqrt(*a));
		++a;
	}
	return B;
}

/////////////////////// discrete operators /////////////////////////

template <class T>
inline vector<T> vDiff(const vector<T>& A){
	unsigned int n = A.size();
	vector<T> B;
	B.reserve(n-1);
	for (unsigned int i=1; i<n; i++)
		B.push_back(A[i] - A[i-1]);
	return B;
}

template <class T>
inline vector<T> vDelta(const vector<T>& A){
	unsigned int n = A.size();
	vector<T> B;
	B.reserve(n+1);
	B.push_back(A[0]);
	for (unsigned int i=1; i<n; i++)
		B.push_back(A[i] - A[i-1]);
	B.push_back(-A[n-1]);
	return B;
}

template <class T>
inline vector<T> vMidAvg(const vector<T>& A){
	vector<T> B;
	unsigned int n = A.size();
	B.reserve(n+1);
	B.push_back(A[0]/2);
	for (unsigned int i=1; i<n; i++)
		B.push_back((A[i]+A[i-1])/2);
	B.push_back(A[n-1]/2);
	return B;
}

template <class T>
inline vector<T> vMidAvgInterior(const vector<T>& A){
	vector<T> B;
	unsigned int n = A.size();
	B.reserve(n-1);
	for (unsigned int i=1; i<n; i++)
		B.push_back((A[i]+A[i-1])/2);
	return B;
}

template <class T>
inline T vSum(const vector<T>& A){
	typename vector<T>::const_iterator a = A.begin();
	T b = *a;
	++a;
	while(a != A.end()){
		b += *a;
		++a;
	}
	return b;
}

//  returns the sums of A[start] through A[end-1]
template <class T>
inline T vPartialSum(const vector<T>& A, int start, int end){
	T result = T();
	for (int i=start; i<end; i++)
		result += A[i];
	return result;
}

// returns the cumulative sum vector of a vector
template <class T>
inline vector<T> vCumSum(const vector<T>& A){
	vector<T> B;
	unsigned int n = A.size();
	B.reserve(n+1);
	B.push_back(A[0]);
	for (unsigned int i=1; i<n; i++)
		B.push_back(A[i]+B[i-1]);
	return B;
}

template <class T>
inline T vMean(const vector<T>& A){
	return vSum(A)/A.size();
}

template <class T>
inline T vMax(const vector<T>& A){
	typename vector<T>::const_iterator a = A.begin();
	T b = *a;
	++a;
	while(a != A.end()){
		if (*a > b)
			b = *a;
		++a;
	}
	return b;
}

template <class T>
inline T vMin(const vector<T>& A){
	typename vector<T>::const_iterator a = A.begin();
	T b = *a;
	++a;
	while(a != A.end()){
		if (*a < b)
			b = *a;
		++a;
	}
	return b;
}


inline VEF::VREAL vAbs(const VEF::VREAL& A){
	VEF::VREAL B;
	B.reserve(A.size());
	VEF::VREAL::const_iterator a = A.begin();
	while(a != A.end()){
		B.push_back(fabs(*a));
		++a;
	}
	return B;
}


///////////////////// discrete generators ////////////////////////////

template <class T>
inline vector<T> vRange(T start, T end, int steps){
	vector<T> A;
	A.reserve(steps);
	if (steps>1){
		T stepSize = (end-start)/(steps-1);
		for (int i=0; i<steps; i++)
			A.push_back(start + stepSize*i);
	}
	else if (steps == 1)
		A.push_back(start);
	else
		throw("vRange must be called with int steps > 0");
	return A;
}

//////////////////// recombination operators ////////////////////////////

template <class T>
inline vector<T> vCat(const vector<T>& A, const vector<T>& B){
	vector<T> C;
	C.reserve(A.size() + B.size());
	typename vector<T>::const_iterator a = A.begin();
	while (a != A.end()){
		C.push_back(*a);
		a++;
	}
	typename vector<T>::const_iterator b = B.begin();
	while (b != B.end()){
		C.push_back(*b);
		b++;
	}
	return C;
}

#endif
