

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
	//Dot production;

	
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
	



};




#endif // !MLS_POINTTPL_H

 