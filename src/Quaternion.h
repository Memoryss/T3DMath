#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include "Vector.h"

namespace T3D {

	/************************************************************************/
	/* 四元数                                                                     */
	/************************************************************************/
	
	class Matrix33;

	class Quaternion
	{
	public:
		float w; //实部
		float x; // i
		float y; // j
		float z; // k

		Quaternion();

		Quaternion(float fw, float fx, float fy, float fz);

		//从旋转矩阵中，构建四元数
		Quaternion(const Matrix33 &rot);

		//从角度和坐标轴中构建四元数
		Quaternion(const float degree, const Vec3 &v);

		//从坐标轴中构建四元数
		Quaternion(const Vec3 &xaxis, const Vec3 &yaxis, const Vec3 &zaxis);

		void swap(Quaternion &quat);

		void Identity();

		float operator[] (const size_t index) const;
		float & operator[] (const size_t index);

		void FromRotationMatrix(const Matrix33 &mat);
		void ToRotationMatrix(Matrix33 &mat) const;

		//统一使用角度表示 如果需要弧度，自行转化
		void FromAngleAxis(const float degree, const Vec3 &rAxis);
		void ToAngleAxis(float &degree, Vec3 &rAxis);

		void FromAxes(const Vec3 &xAxis, const Vec3 &yAxis, const Vec3 &zAxis);
		void ToAxes(Vec3 &xAxis, Vec3 &yAxis, Vec3 &zAxis);

		Quaternion & operator=(const Quaternion &quat);
		Quaternion operator+(const Quaternion &quat) const;
		Quaternion operator-(const Quaternion &quat) const;
		Quaternion operator*(const Quaternion &quat) const;
		Quaternion operator*(const float scale) const;
		Quaternion operator-() const;
		bool operator==(const Quaternion &quat) const;

		float Dot(const Quaternion &quat) const;
		float Length() const;
		float Normalize();
		Quaternion Inverse() const;

		/// Rotation of a vector by a quaternion
		Vec3 operator*(const Vec3 &v) const;

		float getRoll() const;
		float getPitch() const;
		float getYaw() const;

		static const Quaternion IDENTITY;
	};

}

#endif
