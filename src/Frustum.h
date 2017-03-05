#ifndef __FRUSTUM_H__
#define __FRUSTUM_H__

#include "Plane.h"
#include "Matrix.h"

namespace T3D {

	enum ProjectionType
	{
		PT_ORTHOGRAPHIC, //正交投影
		PT_PERSPECTIVE, //透视投影
	};

	enum FrustumPlaneType
	{
		FR_PLANE_NEAR = 0,
		FR_PLANE_FAR,
		FR_PLANE_RIGHT,
		FR_PLANE_LEFT,
		FR_PLANE_TOP,
		FR_PLANE_BOTTOM,
		FRUSTUM_PLANES
	};

	/************************************************************************/
	/* 视锥体                                                                     */
	/************************************************************************/
	class Frustum
	{
	public:
		Frustum();
		Frustum(const Frustum &rhs);
		Frustum(float near, float far, float aspect, float yfov, ProjectionType type = PT_PERSPECTIVE);

		virtual ~Frustum();

		Frustum & operator=(const Frustum &rhs);

		bool operator==(const Frustum &rhs);

		void SetProjType(ProjectionType type);

		void SetFov(float degree);

		void SetNearDist(float dist);

		void SetFarDist(float dist);

		void SetAspect(float aspect);

		const Plane & GetPlane(FrustumPlaneType type);

		const Matrix44 & GetViewMatrix() const;
		const Matrix44 & GetProjectionMatrix() const;
		void GetViewProjectionMatrix(Matrix44 &mat) const;

	protected:
		void UpdateView(const Vec3 &pos, const Quaternion &orientation);

		void UpdateProj();

	protected:
		ProjectionType m_projType; //投影类型

		float m_fov; //y方向的fov
		float m_far; //远裁剪面的距离
		float m_near; //近裁剪面的距离
		float m_aspect; //宽高比

		Matrix44 m_viewMatrix;
		Matrix44 m_projectionMatrix;

		Plane m_planes[FRUSTUM_PLANES];
	};
}
#endif // !

