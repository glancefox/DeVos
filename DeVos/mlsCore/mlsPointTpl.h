

#ifndef MLS_POINTTPL_H
#define MLS_POINTTPL_H


#include<math.h>
#include<limits>
#include<algorithm>

template<typename Type> class Point2DTpl
{
public:
	union 
	{
		struct
		{
			Type x, y;
		};
		Type u[2];
	};

	// Default inline constructor;
	// Initiate the point to (0,0);
	inline explicit Point2DTpl() : x(0), y(0) {}

	//Constructor from a 2D point;
	inline Point2DTpl(Type a, Type b) : x(a), y(b) {}

	//Calculate the Square Norm;
	inline Type norm2() const { return (x*x) + (y*y); }

	//Returns vector norm;
	inline Type norm() const { return std::sqrt(norm2()); }

	//Normalization the norm to unity;
	inline void normalize()
	{
		Type n = norm2();
		if (n > 0)
		{
			*this /= std::sqrt(n);
		}
	}


	//Dot product;
	inline Type dot(const Point2DTpl& pt) const
	{
		return this->x*pt.x + this->y*pt.y;
	}
	//Cross product;
	inline Type cross(const Point2DTpl& pt) const
	{
		return this->x*pt.y - this->y*pt.x;
	}
	//Inverse operator;
	inline Point2DTpl& operator - ()
	{
		this->x = -x;
		this->y = -y;
		return *this;
	}
	//Inplace addition operator;
	inline Point2DTpl& operator +=(const Point2DTpl& pt)
	{
		this->x += pt.x;
		this->y += pt.y;
		return *this;
	}
	//Inplace subtraction operator;
	inline Point2DTpl& operator -=(const Point2DTpl& pt)
	{
		this->x -= pt.x;
		this->y -= pt.y;
		return *this;
	}
	//Inplace multiplication operator;
	inline Point2DTpl& operator*= (Type scalar)
	{
		this->x *= scalar;
		this->y *= scalar;
		return *this;
	}
	//Inplace division by a scalar;
	inline Point2DTpl& operator/=(Type scalar)
	{
		this->x /= scalar;
		this->y /= scalar;
		return *this;
	}
	//Addition operator;
	inline Point2DTpl operator+(const Point2DTpl& pt)const
	{
		return Point2DTpl(this->x + pt.x, this->y + pt.y);
	}
	//Substraction operator£»
	inline Point2DTpl operator-(const Point2DTpl& pt) const
	{
		return Point2DTpl(this->x - pt.x, this->y - pt.y);
	}
	//Multiplication operator;
	inline Point2DTpl operator* (Type scalar) const
	{
		return Point2DTpl(this->x*scalar, this->y*scalar);
	}
	//Division operator;
	inline Point2DTpl operator/ (Type scalar) const
	{
		return Point2DTpl(this->x / scalar, this->y / scalar);
	}
	//Direct coordinates access;
	inline Type& operator[] (unsigned i)
	{
		return this->u[i];
	}
	//Director coordinate access (const)
	inline const Type& operator[](unsigned i) const
	{
		return this->u[i];
	}
};

//This is the code from CloudCompare library.
template <class Type> class Tuple3Tpl
{
public:

	// The 3 tuple values as a union (array/separate values)
	union
	{
		struct
		{
			Type x, y, z;
		};
		Type u[3];
	};

	//! Default constructor
	/** Inits tuple to (0, 0, 0).
	**/
	inline Tuple3Tpl() : x(0), y(0), z(0) {}

	//! Constructor from a triplet of values
	/** Inits typle to (a,b,c).
	**/
	inline Tuple3Tpl(Type a, Type b, Type c) : x(a), y(b), z(c) {}

	//! Constructor from an array of 3 elements
	inline explicit Tuple3Tpl(const Type p[]) : x(p[0]), y(p[1]), z(p[2]) {}

	//! Inverse operator
	inline Tuple3Tpl operator - () const { Tuple3Tpl V(-x, -y, -z); return V; }
	//! In-place addition operator
	inline Tuple3Tpl& operator += (const Tuple3Tpl& v) { x += v.x; y += v.y; z += v.z; return *this; }
	//! In-place subtraction operator
	inline Tuple3Tpl& operator -= (const Tuple3Tpl& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
	//! In-place multiplication (by a scalar) operator
	inline Tuple3Tpl& operator *= (Type v) { x *= v; y *= v; z *= v; return *this; }
	//! In-place division (by a scalar) operator
	inline Tuple3Tpl& operator /= (Type v) { x /= v; y /= v; z /= v; return *this; }
	//! Addition operator
	inline Tuple3Tpl operator + (const Tuple3Tpl& v) const { return Tuple3Tpl(x + v.x, y + v.y, z + v.z); }
	//! Subtraction operator
	inline Tuple3Tpl operator - (const Tuple3Tpl& v) const { return Tuple3Tpl(x - v.x, y - v.y, z - v.z); }
	//! Multiplication operator
	inline Tuple3Tpl operator * (Type s) const { return Tuple3Tpl(x*s, y*s, z*s); }
	//! Division operator
	inline Tuple3Tpl operator / (Type s) const { return Tuple3Tpl(x / s, y / s, z / s); }
};




template <typename Type> class Point3DTpl
{
public:

	//Also organize the 3 tuple values as a union;
	union 
	{
		struct 
		{
			Type x, y, z;
		};
		Type u[3];
	};


	//Default constructor;
	inline Point3DTpl() : x(0), y(0), z(0) {}

	//Constructor from a triplet of values;
	inline Point3DTpl(Type a, Type b, Type c)
		: x(a)
		, y(b)
		, z(c)
	{
	}

	//Constructor from an array of 3 elements;
	inline explicit Point3DTpl(const Type p[])
		: x(p[0])
		, y(p[1])
		, z(p[2])
	{
	}

	//Inverse operator
	inline Point3DTpl operator-() const
	{
		return Point3DTpl(-this->x, -this->y, -this->z);
	}
	//Inplace addition operator for two points;
	inline Point3DTpl& operator+=(const Point3DTpl& pt)
	{
		this->x += pt.x;
		this->y += pt.y;
		this->z += pt.z;
		return *this;
	}
	//Inplace subtraction operator for two points
	inline Point3DTpl& operator -=(const Point3DTpl& pt)
	{
		this->x -= pt.x;
		this->y -= pt.y;
		this->z -= pt.z;
		return *this;
	}
	//Inplace multiplication a scalar to a 3D point;
	inline Point3DTpl& operator *= (Type scalar)
	{
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;
		return *this;
	}
	//Inplace subdivision a scalar to a 3D point;
	inline Point3DTpl& operator /= (Type scalar)
	{
		this->x /= scalar;
		this->y /= scalar;
		this->z /= scalar;
		return *this;
	}
	//Dot product;
	inline Type dot(const Point3DTPl& pt)const
	{
		return this->x*pt.x + this->y*pt.y + this->z*pt.z;
	}
	//Cross product;
	inline Point3DTpl cross(const Point3DTpl& pt)const
	{
		return Point3DTpl(this->y*pt.z - this->z*pt.y, this->z*pt.x - this->x*pt.z, this->x*pt.y - this->y*pt.x);
	}
	//Vector square norm;
	inline  Type norm2() const
	{
		return x*x + y*y + z*z;
	}
	//! Returns vector square norm (forces double precision output)
	inline double norm2d() const 
	{ 
		return static_cast<double>(x)*x + static_cast<double>(y)*y + static_cast<double>(z)*z; 
	}
	//Vector norm;
	inline Type norm() const 
	{ 
		return static_cast<Type>(std::sqrt(norm2d())); 
	}
	//Vector norm(forces double precision output)
	inline double normd() const 
	{ 
		return std::sqrt(norm2d()); 
	}
	//Sets vector norm to unity
	inline void normalize() 
	{ 
		Type n = norm(); 
		if (n > std::numeric_limits<Type>::epsilon()) 
			*this /= n; 
	}
	//! Returns a normalized vector which is orthogonal to this one
	inline Point3DTpl orthogonal() const
	{
		Point3DTpl orth; 
		vorthogonal(u, orth.u); 
		return orth; 
	}


	//Addition operator for two 3D points;
	inline Point3DTpl operator +(const Point3DTpl& pt) const
	{
		return Point3DTpl(this->x + pt.x, this->y + pt.y, this->z + pt.z);
	}
	//Subtraction opereator for two 3D points;
	inline Point3DTpl operator-(const Point3DTpl& pt) const
	{
		return Point3DTpl(this->x - pt.x, this->y - pt.y, this->z - pt.z);
	}
	//Multiplication operator for a scalar and a 3D point;
	inline Point3DTpl operator* (Type scalar)const
	{
		return Point3DTpl(this->x*scalar, this->y*scalar, this->z*scalar);
	}
	
	//Divide 
	static inline void vdivide(Type p[], Type scalar) 
	{ 
		p[0] /= scalar;
		p[1] /= scalar;
		p[2] /= scalar;
	}
	static inline void vdivide(const Type p[], Type scalar, Type r[]) 
	{ 
		r[0] = p[0] / scalar;
		r[1] = p[1] / scalar;
		r[2] = p[2] / scalar;
	}
	static inline void vnormalize(Type p[]) 
	{
		Type n = vnorm2(p); 
		if (n > 0) 
			vdivide(p, std::sqrt(n)); 
	}
	static inline void vmultiply(const Type p[], Type scalar, Type r[]) 
	{
		r[0] = p[0] * scalar;
		r[1] = p[1] * scalar;
		r[2] = p[2] * scalar;
	}
	static inline void vmultiply(Type p[], Type scalar) 
	{ 
		p[0] *= scalar;
		p[1] *= scalar;
		p[2] *= scalar;
	}
	static inline Type vdot(const Type p[], const Type q[]) 
	{ 
		return (p[0] * q[0]) + (p[1] * q[1]) + (p[2] * q[2]); 
	}
	static inline void vcross(const Type p[], const Type q[], Type r[]) 
	{ 
		r[0] = (p[1] * q[2]) - (p[2] * q[1]); 
		r[1] = (p[2] * q[0]) - (p[0] * q[2]); 
		r[2] = (p[0] * q[1]) - (p[1] * q[0]); 
	}
	static inline void vcopy(const Type p[], Type q[]) 
	{ 
		q[0] = p[0]; 
		q[1] = p[1]; 
		q[2] = p[2]; 
	}
	static inline void vset(Type p[], Type scalar) 
	{
		p[0] = p[1] = p[2] = scalar;
	}
	static inline void vset(Type p[], Type x, Type y, Type z) 
	{
		p[0] = x; 
		p[1] = y; 
		p[2] = z; 
	}
	static inline void vadd(const Type p[], const Type q[], Type r[]) 
	{
		r[0] = p[0] + q[0]; 
		r[1] = p[1] + q[1]; 
		r[2] = p[2] + q[2]; 
	}
	// note misspelling: should be vsubtract
	static inline void vsubstract(const Type p[], const Type q[], Type r[]) 
	{
		r[0] = p[0] - q[0]; 
		r[1] = p[1] - q[1]; 
		r[2] = p[2] - q[2]; 
	}
	static inline void vcombination(Type a, const Type p[], Type b, const Type q[], Type r[]) 
	{ 
		r[0] = (a*p[0]) + (b*q[0]); 
		r[1] = (a*p[1]) + (b*q[1]); 
		r[2] = (a*p[2]) + (b*q[2]); 
	}
	static inline void vcombination(const Type p[], Type b, const Type q[], Type r[]) 
	{ 
		r[0] = p[0] + (b*q[0]); 
		r[1] = p[1] + (b*q[1]); 
		r[2] = p[2] + (b*q[2]); 
	}
	static inline Type vnorm2(const Type p[])
	{
		return (p[0] * p[0]) + (p[1] * p[1]) + (p[2] * p[2]);
	}
	static inline void vnormalize(Type& p[]) 
	{
		Type n = vnorm2(p); 
		if (n > 0)
		{
			vdivide(p, std::sqrt(n));
		}
	}
	static inline Type vdistance2(const Type p[], const Type q[]) 
	{
		return ((p[0] - q[0])*(p[0] - q[0])) + ((p[1] - q[1])*(p[1] - q[1])) + ((p[2] - q[2])*(p[2] - q[2])); 
	}
	static inline Type vnorm(const Type p[]) 
	{
		return std::sqrt(vnorm2(p)); 
	}
	static inline Type vdistance(const Type p[], const Type q[]) 
	{
		return std::sqrt(vdistance2(p, q)); 
	}
	static void vorthogonal(const Type p[], Type q[])
	{
		if (std::abs(p[0]) <= std::abs(p[1]) && std::abs(p[0]) <= std::abs(p[2]))
		{
			q[0] = 0; q[1] = p[2]; q[2] = -p[1];
		}
		else if (std::abs(p[1]) <= std::abs(p[0]) && std::abs(p[1]) <= std::abs(p[2]))
		{
			q[0] = -p[2]; q[1] = 0; q[2] = p[0];
		}
		else
		{
			q[0] = p[1]; q[1] = -p[0]; q[2] = 0;
		}
		vnormalize(q);
	}
	static Type vangle_rad(const Type p[], const Type q[])
	{
		Type productNorm = vnorm(p) * vnorm(q);
		if (productNorm < std::numeric_limits<Type>::epsilon())
		{
			return std::numeric_limits<Type>::quiet_NaN();
		}

		Type cosAngle = vdot(p, q) / productNorm;
		return acos(std::max(std::min(cosAngle, static_cast<Type>(1.0)), static_cast<Type>(-1.0)));
	}
};


#endif // !MLS_POINTTPL_H

 