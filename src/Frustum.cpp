#include "Frustum.h"

#include <math.h>

namespace T3D {

	Frustum::Frustum() : m_near(0.1f), m_far(1000.0f), m_fov(90), m_aspect(1.3333f), m_projType(PT_PERSPECTIVE)
	{

	}

	Frustum::Frustum(const Frustum &rhs)
	{
		m_near = rhs.m_near;
		m_far = rhs.m_far;
		m_fov = rhs.m_fov;
		m_aspect = rhs.m_aspect;
		m_projType = rhs.m_projType;
	}

	Frustum::Frustum(float near, float far, float aspect, float yfov, ProjectionType type/* = PT_PERSPECTIVE*/) : m_near(near), m_far(far), m_aspect(aspect), m_fov(yfov), m_projType(type)
	{

	}

	Frustum::~Frustum()
	{

	}

	void Frustum::SetProjType(ProjectionType type)
	{
		m_projType = type;
	}

	void Frustum::SetFov(float degree)
	{
		m_fov = degree;
	}

	void Frustum::SetNearDist(float dist)
	{
		m_near = dist;
	}

	void Frustum::SetFarDist(float dist)
	{
		m_far = dist;
	}

	void Frustum::SetAspect(float aspect)
	{
		m_aspect = aspect;
	}

	const T3D::Plane & Frustum::GetPlane(FrustumPlaneType type)
	{
		return m_planes[type];
	}

	const T3D::Matrix44 & Frustum::GetViewMatrix() const
	{
		return m_viewMatrix;
	}

	const T3D::Matrix44 & Frustum::GetProjectionMatrix() const
	{
		return m_projectionMatrix;
	}

	void Frustum::GetViewProjectionMatrix(Matrix44 &mat) const
	{
		mat = m_projectionMatrix * m_viewMatrix;
	}

	const T3D::Vec3 Frustum::GetCorners(size_t index)
	{
		return m_corners[index];
	}

	void Frustum::UpdateProj()
	{
		//Ӧ���ȸ�����׶������Ͷ���Ȼ���ڸ���ͶӰ����
		//ͨ���仯����ʹ��׶���ǶԳƵ�
		float dist_inv = 1 / (m_far - m_near);

		m_projectionMatrix = Matrix44::IDENTITY;
		m_projectionMatrix[0][0] = m_near / m_corners[0].x;
		m_projectionMatrix[1][1] = m_near / m_corners[0].y;
		m_projectionMatrix[2][2] = -(m_far + m_near) * dist_inv;
		m_projectionMatrix[2][3] = -2 * m_far * m_near * dist_inv;
		m_projectionMatrix[3][2] = 1;
	}

	void Frustum::UpdateOrtho()
	{
		float dist_inv = 1 / (m_far - m_near);

		m_projectionMatrix = Matrix44::IDENTITY;
		m_projectionMatrix[0][0] = 1 / m_corners[0].x;
		m_projectionMatrix[1][1] = 1 / m_corners[0].y;
		m_projectionMatrix[2][2] = -2 * dist_inv;
		m_projectionMatrix[2][3] = -(m_far + m_near) * dist_inv;
	}

	void Frustum::UpdateView(const Vec3 &pos, const Quaternion &orientation)
	{
		// TODO
		//view������Կ��ɣ�ͨ���仯������Ƶ�����ռ�������ԭ�㲢�ҷ������غ�
		Matrix44 trans;
		trans.makeTrans(-pos);  //�Ƚ�����ƶ�����������ϵ��ԭ��
		
		Matrix44 rot(orientation);  //������Ԫ�������仯����
		Matrix44 rot_inv = rot.inverse();  //�任��������

		m_viewMatrix = rot_inv * trans;   //�Ƚ������ƶ�����������ϵԭ�㣬�ڽ�����ת
	}

	void Frustum::Update(const Vec3 &pos, const Quaternion &orientation)
	{
		//�����������
		UpdateView(pos, orientation);

		//���¶�������
		float yAngle = m_fov * 0.5f;
		float tanY = tanf(yAngle);

		float xAngle = atanf(m_aspect * tanY);
		float tanX = tanf(xAngle);

		float nearX = m_near * tanX;
		float farX = m_far * tanX;

		float nearY = m_near * tanY;
		float farY = m_far * tanY;

		float nearZ = m_near;
		float farZ = m_far;

		

		Vec3 corners[8];
		corners[0].x = nearX;
		corners[0].y = nearY;
		corners[0].z = nearZ;

		corners[1].x = -nearX;
		corners[1].y = nearY;
		corners[1].z = nearZ;

		corners[2].x = -nearX;
		corners[2].y = -nearY;
		corners[2].z = nearZ;

		corners[3].x = nearX;
		corners[3].y = -nearY;
		corners[3].z = nearZ;

		corners[4].x = farX;
		corners[4].y = farY;
		corners[4].z = farZ;

		corners[5].x = -farX;
		corners[5].y = farY;
		corners[5].z = farZ;

		corners[6].x = -farX;
		corners[6].y = -farY;
		corners[6].z = farZ;

		corners[7].x = farX;
		corners[7].y = -farY;
		corners[7].z = farZ;

		//������׶�嶥������
		for (size_t index = 0; index < 8; ++index)
		{
			corners[index] = m_viewMatrix * corners[index];  //��updateview֮�����
		}

		//����ͶӰ����
		if (m_projType == PT_PERSPECTIVE)
		{
			UpdateProj();
		}
		else
		{
			UpdateOrtho();
		}
	}

	bool Frustum::operator==(const Frustum &rhs)
	{
		return (m_projType == rhs.m_projType && m_aspect == rhs.m_aspect && m_far == rhs.m_far && m_near == rhs.m_near && m_fov == rhs.m_fov);
	}

	T3D::Frustum & Frustum::operator=(const Frustum &rhs)
	{
		m_near = rhs.m_near;
		m_far = rhs.m_far;
		m_fov = rhs.m_fov;
		m_aspect = rhs.m_aspect;
		m_projType = rhs.m_projType;
		
		return *this;
	}

}