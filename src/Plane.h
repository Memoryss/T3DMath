#ifndef __PLANE_H__
#define __PLANE_H__

#include "Vector.h"

namespace T3D {

	/************************************************************************/
	/* ƽ�淽��                                                                     */
	/************************************************************************/

	//�����ƽ��Ĺ�ϵ
	enum PlaneSide
	{
		PS_Self, //��ƽ����
		PS_Positive, //��ƽ������
		PS_Negative, //��ƽ�渺��
	};

	class Plane
	{
	public:
		Plane();
		Plane(const Plane &rhs);
		Plane(const Vec3 &normal, const float distance);   //���ߺͷ����ϵ�ԭ��ľ���
		Plane(const Vec3 &normal, const Vec3 &point); //���ߺ�ƽ���ϵĵ�
		Plane(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3); //ƽ���ϵ�������

		//�жϵ���ƽ�����һ��
		PlaneSide GetSide(const Vec3 &point);

		//�㵽ƽ��ľ���
		float GetDistance(const Vec3 &point);

		//���������淨�ߺ;���
		void ReCalculate(const Vec3 &normal, const float distance);
		void ReCalculate(const Vec3 &normal, const Vec3 &point);
		void ReCalculate(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3);

		void Normalize();

		Vec3 m_normal; //����
		float m_distance; //��ԭ��ľ��루�Ӵ��㵽ԭ�������������뷨�߷����෴��Ϊ������֮Ϊ����

		bool operator==(const Plane &rhs) const;
	};
}

#endif