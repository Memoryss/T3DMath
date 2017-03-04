#include "Frustum.h"

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

	void Frustum::UpdateFrustum(const Vec3 &pos, const Quaternion &orientation)
	{
		// TODO
		//view矩阵可以看成，通过变化将相机移到世界空间的坐标的原点并且方向轴重合
		Matrix44 trans;
		trans.makeTrans(-pos);  //先将相机移动到世界坐标系的原点
		
		Matrix44 rot(orientation);  //根据四元数创建变化矩阵
		Matrix44 rot_inv = rot.inverse();  //变换矩阵求逆
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