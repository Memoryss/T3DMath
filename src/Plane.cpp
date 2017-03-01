#include "Plane.h"

namespace T3D {



	Plane::Plane() : m_normal(0, 0, 0), m_distance(0.0f)
	{

	}

	Plane::Plane(const Plane &rhs) : m_normal(rhs.m_normal) ,m_distance(rhs.m_distance)
	{
		
	}

	Plane::Plane(const Vec3 &normal, const float distance) : m_normal(normal), m_distance(distance)
	{
		ReCalculate(normal, distance);
	}

	Plane::Plane(const Vec3 &normal, const Vec3 &point)
	{
		ReCalculate(normal, point);
	}

	Plane::Plane(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3)
	{
		ReCalculate(p1, p2, p3);
	}

	PlaneSide Plane::GetSide(const Vec3 &point)
	{
		float distance = GetDistance(point);
		if (distance > 0.0f)
		{
			return PS_Positive;
		}
		else if (distance < 0.0f)
		{
			return PS_Negative;
		}
		else
		{
			return PS_Self;
		}
	}

	float Plane::GetDistance(const Vec3 &point)
	{
		return m_normal.Dot(point) + m_distance;
	}

	void Plane::ReCalculate(const Vec3 &normal, const float distance)
	{
		m_normal = normal;
		m_distance = distance;
	}

	void Plane::ReCalculate(const Vec3 &normal, const Vec3 &point)
	{
		m_normal = normal;
		m_distance = -(m_normal.Dot(point));
	}

	void Plane::ReCalculate(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3)
	{
		Vec3 edge1 = p2 - p1;
		Vec3 edge2 = p3 - p1;
		m_normal = edge1.Cross(edge2);
		m_normal.Normalize();
		m_distance = -(m_normal.Dot(p1));
	}

	void Plane::Normalize()
	{
		float len = m_normal.Normalize();
		m_distance /= len;
	}

	bool Plane::operator==(const Plane &rhs) const
	{
		return (rhs.m_distance == m_distance && rhs.m_normal == m_normal);
	}

}