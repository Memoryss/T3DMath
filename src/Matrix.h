#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "Vector.h"
#include "Quaternion.h"

namespace T3D {

	class Matrix33
	{
	public:
		Matrix33();
		explicit Matrix33(float arr[3][3]);
		Matrix33(const Matrix33 &mat);
		Matrix33(float m00, float m01, float m02,
			float m10, float m11, float m12,
			float m20, float m21, float m22);

		Matrix33& operator=(const Matrix33 &mat);

		void swap(Matrix33 &mat);

		float * operator[] (int row);
		const float * operator[] (int row) const;

		void SetColumn(size_t iCol, const Vec3 &v);
		Vec3 GetColumn(size_t iCol);

		bool operator==(const Matrix33 &mat) const;

		Matrix33 operator+(const Matrix33 &mat) const;

		Matrix33 operator-(const Matrix33 &mat) const;

		Matrix33 operator*(const Matrix33 &mat) const;

		Matrix33 operator-() const;

		Vec3 operator*(const Vec3 &v) const;

		Matrix33 operator*(const float scale) const;

		//矩阵转置
		Matrix33 Transpose() const;

		//逆矩阵
		bool Inverse(Matrix33 &mat, float tolerance = 1e-06) const;

		//行列式
		float Determinant() const;

	protected:
		float m[3][3] = { 0 };
	};

	class Matrix44
	{
	public:
		Matrix44();
		Matrix44(float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23,
			float m30, float m31, float m32, float m33);

		Matrix44(const Matrix44 &mat);

		Matrix44(const Matrix33 &mat);

		Matrix44(const Quaternion &rot);

		void swap(Matrix44 &mat);

		float * operator[](size_t row);
		const float * operator[](size_t row) const;

		Matrix44 operator*(const Matrix44 &mat) const;
		Vec3 operator*(const Vec3 &v) const;
		Vec4 operator*(const Vec4 &v) const;
		Matrix44 operator*(float scale) const;

		Matrix44 operator+(const Matrix44 &mat) const;
		Matrix44 operator-(const Matrix44 &mat) const;
		bool operator==(const Matrix44 &mat) const;
		void operator=(const Matrix33 &mat3);

		Matrix44 transpose() const;
		float determinant() const;
		Matrix44 inverse() const;

		//平移
		void setTrans(const Vec3 &v);
		Vec3 getTrans();

		void makeTrans(const Vec3 &v);
		void makeTrans(float tx, float ty, float tz);

		static Matrix44 getTrans(const Vec3 &v);
		static Matrix44 getTrans(float tx, float ty, float tz);

		//缩放
		void setScale(const Vec3 &v);
		static Matrix44 getScale(const Vec3 &v);
		static Matrix44 getScale(float tx, float ty, float tz);

		static const Matrix44 ZERO;
		static const Matrix44 ZEROAFFINE;
		static const Matrix44 IDENTITY;

	protected:
		float m[4][4] = {0.0};
		
	};

}

#endif
