

#include"mlsBoundingBox.h"


using namespace mlsCore;


BoundingBox::BoundingBox()
	: m_bbMin(0, 0, 0)
	, m_bbMax(0, 0, 0)
	, m_valid(false)
{
	PointTypeF max = (std::numeric_limits<PointTypeF>::max)();
	//m_bbMax(-max, -max, -max);

}


BoundingBox::BoundingBox(const Point3Df& bbMinCorner, const Point3Df& bbMaxCorner)
	: m_bbMin(bbMinCorner)
	, m_bbMax(bbMaxCorner)
	, m_valid(true)
{}

void BoundingBox::clear()
{
	m_bbMin = m_bbMax = Point3Df(0, 0, 0);
	m_valid = false;
}

Point3Df BoundingBox::getCenter() const
{
	return (m_bbMax + m_bbMin) * static_cast<PointTypeF>(0.5);
}

Point3Df BoundingBox::getDiagVec() const
{
	return (m_bbMax - m_bbMin);
}

PointTypeF BoundingBox::getMinBoxDim() const
{
	Point3Df V = getDiagVec();

	return std::min(V.x, std::min(V.y, V.z));
}

PointTypeF BoundingBox::getMaxBoxDim() const
{
	Point3Df V = getDiagVec();

	return std::max(V.x, std::max(V.y, V.z));
}

double BoundingBox::computeVolume() const
{
	Point3Df V = getDiagVec();

	return static_cast<double>(V.x) * static_cast<double>(V.y) * static_cast<double>(V.z);
}

BoundingBox BoundingBox::operator + (const BoundingBox& aBBox) const
{
	if (!m_valid)
		return aBBox;
	if (!aBBox.isValid())
		return *this;


}

const BoundingBox& BoundingBox::operator += (const BoundingBox& aBBox)
{
	if (aBBox.isValid())
	{
		add(aBBox.minCorner());
		add(aBBox.maxCorner());
	}

	return *this;
}

const BoundingBox& BoundingBox::operator += (const Point3Df& aVector)
{
	if (m_valid)
	{
		m_bbMin += aVector;
		m_bbMax += aVector;
	}

	return *this;
}

const BoundingBox& BoundingBox::operator -= (const Point3Df& aVector)
{
	if (m_valid)
	{
		m_bbMin -= aVector;
		m_bbMax -= aVector;
	}

	return *this;
}

const BoundingBox& BoundingBox::operator *= (float scaleFactor)
{
	if (m_valid)
	{
		m_bbMin *= scaleFactor;
		m_bbMax *= scaleFactor;
	}

	return *this;
}

void BoundingBox::add(const Point3Df& aPoint)
{
	if (m_valid)
	{
		if (aPoint.x < m_bbMin.x)
			m_bbMin.x = aPoint.x;
		else if (aPoint.x > m_bbMax.x)
			m_bbMax.x = aPoint.x;

		if (aPoint.y < m_bbMin.y)
			m_bbMin.y = aPoint.y;
		else if (aPoint.y > m_bbMax.y)
			m_bbMax.y = aPoint.y;

		if (aPoint.z < m_bbMin.z)
			m_bbMin.z = aPoint.z;
		else if (aPoint.z > m_bbMax.z)
			m_bbMax.z = aPoint.z;
	}
	else
	{
		m_bbMax = m_bbMin = aPoint;
		m_valid = true;
	}
}

PointTypeF BoundingBox::minDistTo(const BoundingBox& box) const
{
	if (m_valid && box.isValid())
	{
		Point3Df d(0, 0, 0);

		for (uint8_t dim = 0; dim < 3; ++dim)
		{
			//if the boxes overlap in one dimension, the distance is zero (in this dimension)
			if (box.m_bbMin.u[dim] > m_bbMax.u[dim])
				d.u[dim] = box.m_bbMin.u[dim] - m_bbMax.u[dim];
			else if (box.m_bbMax.u[dim] < m_bbMin.u[dim])
				d.u[dim] = m_bbMin.u[dim] - box.m_bbMax.u[dim];
		}

		return d.norm();
	}
	else
	{
		return std::numeric_limits<PointTypeF>::quiet_NaN();
	}
}

const BoundingBox& BoundingBox::operator *= (const SquareMatrix& mat)
{
	if (m_valid)
	{
		Point3Df boxCorners[8];

		boxCorners[0] = m_bbMin;
		boxCorners[1] = Point3Df(m_bbMin.x, m_bbMin.y, m_bbMax.z);
		boxCorners[2] = Point3Df(m_bbMin.x, m_bbMax.y, m_bbMin.z);
		boxCorners[3] = Point3Df(m_bbMax.x, m_bbMin.y, m_bbMin.z);
		boxCorners[4] = m_bbMax;
		boxCorners[5] = Point3Df(m_bbMin.x, m_bbMax.y, m_bbMax.z);
		boxCorners[6] = Point3Df(m_bbMax.x, m_bbMax.y, m_bbMin.z);
		boxCorners[7] = Point3Df(m_bbMax.x, m_bbMin.y, m_bbMax.z);

		clear();

		for (int i = 0; i < 8; ++i)
		{
			add(mat*boxCorners[i]);
		}
	}

	return *this;
}