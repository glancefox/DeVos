

#ifndef MLS_BOUNDINGBOX_H
#define MLS_BOUNDINGBOX_H

#include"mlsSquareMatrix.h"

#include<algorithm>
#include<cstdint>


namespace mlsCore
{
	class BoundingBox
	{
	public:

		//! Default constructor
		BoundingBox(); 
		//! Constructor from two vectors (lower min. and upper max. corners)
		BoundingBox(const Point3Df& minCorner, const Point3Df& maxCorner);

		//! Returns the 'sum' of this bounding-box and another one
		BoundingBox operator + (const BoundingBox& aBBox) const; 
		//! In place 'sum' of this bounding-box with another one
		const BoundingBox& operator += (const BoundingBox& aBBox);
		//! Shifts the bounding box with a vector
		const BoundingBox& operator += (const Point3Df& aVector);
		//! Shifts the bounding box with a vector
		const BoundingBox& operator -= (const Point3Df& aVector);
		//! Scales the bounding box
		const BoundingBox& operator *= (float scaleFactor);
		//! Rotates the bounding box
		const BoundingBox& operator *= (const SquareMatrix& aMatrix);

		//! Resets the bounding box
		/** (0,0,0) --> (0,0,0)
		**/
		void clear();

		//! 'Enlarges' the bounding box with a point
		void add(const Point3Df& aPoint);

		//! Returns min corner (const)
		inline const Point3Df& minCorner() const { return m_bbMin; }
		//! Returns max corner (const)
		inline const Point3Df& maxCorner() const { return m_bbMax; }

		//! Returns min corner
		inline Point3Df& minCorner() { return m_bbMin; }
		//! Returns max corner
		inline Point3Df& maxCorner() { return m_bbMax; }

		//! Returns center
		Point3Df getCenter() const;
		//! Returns diagonal vector
		Point3Df getDiagVec() const;
		//! Returns diagonal length
		inline PointTypeF getDiagNorm() const { return getDiagVec().norm(); }
		//! Returns diagonal length (double precision)
		double getDiagNormd() const { return getDiagVec().normd(); }
		//! Returns minimal box dimension
		PointTypeF getMinBoxDim() const;
		//! Returns maximal box dimension
		PointTypeF getMaxBoxDim() const;
		//! Returns the bounding-box volume
		double computeVolume() const; 

		//! Sets bonding box validity
		inline void setValidity(bool state) { m_valid = state; }

		//! Returns whether bounding box is valid or not
		inline bool isValid() const { return m_valid; }

		//! Computes min gap (absolute distance) between this bounding-box and another one
		/** \return min gap (>=0) or -1 if at least one of the box is not valid
		**/
		PointTypeF minDistTo(const BoundingBox& box) const;

		//! Returns whether a points is inside the box or not
		/** Warning: box should be valid!
		**/
		inline bool contains(const Point3Df& P) const
		{
			return (P.x >= m_bbMin.x && P.x <= m_bbMax.x &&
				P.y >= m_bbMin.y && P.y <= m_bbMax.y &&
				P.z >= m_bbMin.z && P.z <= m_bbMax.z);
		}

	protected:
		//! Lower min. corner
		Point3Df m_bbMin;
		//! Upper max. corner
		Point3Df m_bbMax;
		//! Validity
		bool m_valid;
	};

}//namespace mlsCore

#endif // !MLS_BOUNDINGBOX_H
