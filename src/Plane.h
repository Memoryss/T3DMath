#ifndef __PLANE_H__
#define __PLANE_H__

#include "Vector.h"

namespace T3D {

	/************************************************************************/
	/* 平面方程                                                                     */
	/************************************************************************/

	//顶点和平面的关系
	enum PlaneSide
	{
		PS_Self, //在平面上
		PS_Positive, //在平面正向
		PS_Negative, //在平面负向
	};

	class Plane
	{
	public:
		Plane();
		Plane(const Plane &rhs);
		Plane(const Vec3 &normal, const float distance);   //法线和法线上到原点的距离
		Plane(const Vec3 &normal, const Vec3 &point); //法线和平面上的点
		Plane(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3); //平面上的三个点

		//判断点在平面的哪一侧
		PlaneSide GetSide(const Vec3 &point);

		//点到平面的距离
		float GetDistance(const Vec3 &point);

		//重新设置面法线和距离
		void ReCalculate(const Vec3 &normal, const float distance);
		void ReCalculate(const Vec3 &normal, const Vec3 &point);
		void ReCalculate(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3);

		void Normalize();

		Vec3 m_normal; //法线
		float m_distance; //到原点的距离（从垂点到原点的向量，如果与法线方向相反则为负，反之为正）

		bool operator==(const Plane &rhs) const;
	};
}

#endif