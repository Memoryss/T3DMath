#include "Matrix.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <utility>

namespace T3D {

	Matrix33::Matrix33()
	{

	}

	Matrix33::Matrix33(float arr[3][3])
	{
		memcpy(m, arr, 9 * sizeof(float));
	}

	Matrix33::Matrix33(const Matrix33 &mat)
	{
		memcpy(m, mat.m, 9 * sizeof(float));
	}

	Matrix33::Matrix33(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22)
	{
		m[0][0] = m00;
		m[0][1] = m01;
		m[0][2] = m02;
		m[1][0] = m10;
		m[1][1] = m11;
		m[1][2] = m12;
		m[2][0] = m20;
		m[2][1] = m21;
		m[2][2] = m22;
	}

	void Matrix33::swap(Matrix33 &mat)
	{
		float tempM[3][3];
		memcpy(tempM, mat.m, 9 * sizeof(float));
		memcpy(mat.m, m, 9 * sizeof(float));
		memcpy(m, tempM, 9 * sizeof(float));
	}

	float * Matrix33::operator[](int row)
	{
		assert(row < 3);
		return m[row];
	}

	const float * Matrix33::operator[](int row) const
	{
		assert(row < 3);
		return m[row];
	}

	void Matrix33::SetColumn(size_t iCol, const Vec3 &v)
	{
		m[0][iCol] = v.x;
		m[1][iCol] = v.y;
		m[2][iCol] = v.z;
	}

	Vec3 Matrix33::GetColumn(size_t iCol)
	{
		return Vec3(m[0][iCol], m[1][iCol], m[2][iCol]);
	}

	Matrix33 Matrix33::operator+(const Matrix33 &mat) const
	{
		Matrix33 m33;

		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				m33.m[i][j] = m[i][j] + mat.m[i][j];
			}
		}

		return m33;
	}

	Matrix33 Matrix33::operator-(const Matrix33 &mat) const
	{
		Matrix33 m33;

		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				m33.m[i][j] = m[i][j] - mat.m[i][j];
			}
		}

		return m33;
	}

	Matrix33 Matrix33::operator-() const
	{
		Matrix33 m33;

		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				m33.m[i][j] = -m[i][j];
			}
		}

		return m33;
	}

	Matrix33 Matrix33::operator*(const Matrix33 &mat) const
	{
		Matrix33 m33;

		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				m33.m[i][j] = m[i][0] * mat.m[0][j] + m[i][1] * mat.m[1][j] + m[i][2] * m[2][j];
			}
		}

		return m33;
	}

	Vec3 Matrix33::operator*(const Vec3 &v) const
	{
		Vec3 v3;
		v3.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
		v3.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
		v3.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;

		return v3;
	}

	Matrix33 Matrix33::operator*(const float scale) const
	{
		Matrix33 m33;

		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				m33.m[i][j] = m[i][0] * scale;
			}
		}

		return m33;
	}

	Matrix33 Matrix33::Transpose() const
	{
		Matrix33 m33;

		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				m33.m[i][j] = m[j][i];
			}
		}

		return m33;
	}

	bool Matrix33::Inverse(Matrix33 &rkInverse, float fTolerance /*= 1e-06*/) const
	{
		rkInverse[0][0] = m[1][1] * m[2][2] -
			m[1][2] * m[2][1];
		rkInverse[0][1] = m[0][2] * m[2][1] -
			m[0][1] * m[2][2];
		rkInverse[0][2] = m[0][1] * m[1][2] -
			m[0][2] * m[1][1];
		rkInverse[1][0] = m[1][2] * m[2][0] -
			m[1][0] * m[2][2];
		rkInverse[1][1] = m[0][0] * m[2][2] -
			m[0][2] * m[2][0];
		rkInverse[1][2] = m[0][2] * m[1][0] -
			m[0][0] * m[1][2];
		rkInverse[2][0] = m[1][0] * m[2][1] -
			m[1][1] * m[2][0];
		rkInverse[2][1] = m[0][1] * m[2][0] -
			m[0][0] * m[2][1];
		rkInverse[2][2] = m[0][0] * m[1][1] -
			m[0][1] * m[1][0];

		float fDet =
			m[0][0] * rkInverse[0][0] +
			m[0][1] * rkInverse[1][0] +
			m[0][2] * rkInverse[2][0];

		if (fabs(fDet) <= fTolerance)
			return false;

		float fInvDet = 1.0f / fDet;
		for (size_t iRow = 0; iRow < 3; iRow++)
		{
			for (size_t iCol = 0; iCol < 3; iCol++)
				rkInverse[iRow][iCol] *= fInvDet;
		}

		return true;
	}

	float Matrix33::Determinant() const
	{
		float fCofactor00 = m[1][1] * m[2][2] -
			m[1][2] * m[2][1];
		float fCofactor10 = m[1][2] * m[2][0] -
			m[1][0] * m[2][2];
		float fCofactor20 = m[1][0] * m[2][1] -
			m[1][1] * m[2][0];

		float fDet =
			m[0][0] * fCofactor00 +
			m[0][1] * fCofactor10 +
			m[0][2] * fCofactor20;

		return fDet;
	}

	bool Matrix33::operator==(const Matrix33 &mat) const
	{
		for (size_t row = 0; row < 3; ++row)
		{
			for (size_t col = 0; col < 3; ++col)
			{
				if (m[row][col] != mat.m[row][col]) return false;
			}
		}

		return true;
	}

	Matrix33& Matrix33::operator=(const Matrix33 &mat)
	{
		memcpy(m, mat.m, 9 * sizeof(float));
		return *this;
	}

	Matrix44::Matrix44()
	{

	}

	Matrix44::Matrix44(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33)
	{
		m[0][0] = m00;
		m[0][1] = m01;
		m[0][2] = m02;
		m[0][3] = m03;
		m[1][0] = m10;
		m[1][1] = m11;
		m[1][2] = m12;
		m[1][3] = m13;
		m[2][0] = m20;
		m[2][1] = m21;
		m[2][2] = m22;
		m[2][3] = m23;
		m[3][0] = m30;
		m[3][1] = m31;
		m[3][2] = m32;
		m[3][3] = m33;
	}

	static float MINOR(const Matrix44& m, const size_t r0, const size_t r1, const size_t r2,
			const size_t c0, const size_t c1, const size_t c2)
	{
		return m[r0][c0] * (m[r1][c1] * m[r2][c2] - m[r2][c1] * m[r1][c2]) -
			m[r0][c1] * (m[r1][c0] * m[r2][c2] - m[r2][c0] * m[r1][c2]) +
			m[r0][c2] * (m[r1][c0] * m[r2][c1] - m[r2][c0] * m[r1][c1]);
	}

	const Matrix44 Matrix44::ZERO(
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0);

	const Matrix44 Matrix44::ZEROAFFINE(
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 1);

	const Matrix44 Matrix44::IDENTITY(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);

	Matrix44::Matrix44(const Matrix44 &mat)
	{
		memcpy(m, mat.m, 16 * sizeof(mat));
	}

	Matrix44::Matrix44(const Quaternion &rot)
	{
		Matrix33 m3x3;
		rot.ToRotationMatrix(m3x3);
		operator=(IDENTITY);
		operator=(m3x3);
	}

	Matrix44::Matrix44(const Matrix33 &mat)
	{
		operator=(IDENTITY);
		operator=(mat);
	}

	void Matrix44::swap(Matrix44 &other)
	{
		std::swap(m[0][0], other.m[0][0]);
		std::swap(m[0][1], other.m[0][1]);
		std::swap(m[0][2], other.m[0][2]);
		std::swap(m[0][3], other.m[0][3]);
		std::swap(m[1][0], other.m[1][0]);
		std::swap(m[1][1], other.m[1][1]);
		std::swap(m[1][2], other.m[1][2]);
		std::swap(m[1][3], other.m[1][3]);
		std::swap(m[2][0], other.m[2][0]);
		std::swap(m[2][1], other.m[2][1]);
		std::swap(m[2][2], other.m[2][2]);
		std::swap(m[2][3], other.m[2][3]);
		std::swap(m[3][0], other.m[3][0]);
		std::swap(m[3][1], other.m[3][1]);
		std::swap(m[3][2], other.m[3][2]);
		std::swap(m[3][3], other.m[3][3]);
	}

	float * Matrix44::operator[](size_t row)
	{
		assert(row < 4);
		return m[row];
	}

	const float * Matrix44::operator[](size_t row) const
	{
		assert(row < 4);
		return m[row];
	}

	Matrix44 Matrix44::operator*(const Matrix44 &mat) const
	{
		Matrix44 m4;
		for (size_t row = 0; row < 4; ++row)
		{
			for (size_t col = 0; col < 4; ++col)
			{
				m4[row][col] = m[row][0] * mat.m[0][col] + m[row][1] * mat.m[1][col] + m[row][2] * mat.m[2][col] + m[row][3] * mat.m[3][col];
			}
		}

		return m4;
	}

	Vec3 Matrix44::operator*(const Vec3 &v) const
	{
		//¸½¼ÓÍ¸ÊÓ³ý·¨
		Vec3 r;

		float fInvW = 1.0f / (m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3]);

		r.x = (m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3]) * fInvW;
		r.y = (m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3]) * fInvW;
		r.z = (m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3]) * fInvW;

		return r;
	}

	Vec4 Matrix44::operator*(const Vec4 &v) const
	{
		return Vec4(
			m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w,
			m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w,
			m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w,
			m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w
		);
	}

	Matrix44 Matrix44::operator*(float scalar) const
	{
		return Matrix44(
			scalar*m[0][0], scalar*m[0][1], scalar*m[0][2], scalar*m[0][3],
			scalar*m[1][0], scalar*m[1][1], scalar*m[1][2], scalar*m[1][3],
			scalar*m[2][0], scalar*m[2][1], scalar*m[2][2], scalar*m[2][3],
			scalar*m[3][0], scalar*m[3][1], scalar*m[3][2], scalar*m[3][3]);
	}

	Matrix44 Matrix44::operator+(const Matrix44 &m2) const
	{
		Matrix44 r;

		r.m[0][0] = m[0][0] + m2.m[0][0];
		r.m[0][1] = m[0][1] + m2.m[0][1];
		r.m[0][2] = m[0][2] + m2.m[0][2];
		r.m[0][3] = m[0][3] + m2.m[0][3];

		r.m[1][0] = m[1][0] + m2.m[1][0];
		r.m[1][1] = m[1][1] + m2.m[1][1];
		r.m[1][2] = m[1][2] + m2.m[1][2];
		r.m[1][3] = m[1][3] + m2.m[1][3];

		r.m[2][0] = m[2][0] + m2.m[2][0];
		r.m[2][1] = m[2][1] + m2.m[2][1];
		r.m[2][2] = m[2][2] + m2.m[2][2];
		r.m[2][3] = m[2][3] + m2.m[2][3];

		r.m[3][0] = m[3][0] + m2.m[3][0];
		r.m[3][1] = m[3][1] + m2.m[3][1];
		r.m[3][2] = m[3][2] + m2.m[3][2];
		r.m[3][3] = m[3][3] + m2.m[3][3];

		return r;
	}

	Matrix44 Matrix44::operator-(const Matrix44 &m2) const
	{
		Matrix44 r;
		r.m[0][0] = m[0][0] - m2.m[0][0];
		r.m[0][1] = m[0][1] - m2.m[0][1];
		r.m[0][2] = m[0][2] - m2.m[0][2];
		r.m[0][3] = m[0][3] - m2.m[0][3];

		r.m[1][0] = m[1][0] - m2.m[1][0];
		r.m[1][1] = m[1][1] - m2.m[1][1];
		r.m[1][2] = m[1][2] - m2.m[1][2];
		r.m[1][3] = m[1][3] - m2.m[1][3];

		r.m[2][0] = m[2][0] - m2.m[2][0];
		r.m[2][1] = m[2][1] - m2.m[2][1];
		r.m[2][2] = m[2][2] - m2.m[2][2];
		r.m[2][3] = m[2][3] - m2.m[2][3];

		r.m[3][0] = m[3][0] - m2.m[3][0];
		r.m[3][1] = m[3][1] - m2.m[3][1];
		r.m[3][2] = m[3][2] - m2.m[3][2];
		r.m[3][3] = m[3][3] - m2.m[3][3];

		return r;
	}

	void Matrix44::operator=(const Matrix33 &mat3)
	{
			m[0][0] = mat3[0][0]; m[0][1] = mat3[0][1]; m[0][2] = mat3[0][2];
			m[1][0] = mat3[1][0]; m[1][1] = mat3[1][1]; m[1][2] = mat3[1][2];
			m[2][0] = mat3[2][0]; m[2][1] = mat3[2][1]; m[2][2] = mat3[2][2];
	}

	Matrix44 Matrix44::transpose() const
	{
		return Matrix44(m[0][0], m[1][0], m[2][0], m[3][0],
			m[0][1], m[1][1], m[2][1], m[3][1],
			m[0][2], m[1][2], m[2][2], m[3][2],
			m[0][3], m[1][3], m[2][3], m[3][3]);
	}

	float Matrix44::determinant() const
	{
		return m[0][0] * MINOR(*this, 1, 2, 3, 1, 2, 3) -
			m[0][1] * MINOR(*this, 1, 2, 3, 0, 2, 3) +
			m[0][2] * MINOR(*this, 1, 2, 3, 0, 1, 3) -
			m[0][3] * MINOR(*this, 1, 2, 3, 0, 1, 2);
	}

	Matrix44 Matrix44::inverse() const
	{
		float m00 = m[0][0], m01 = m[0][1], m02 = m[0][2], m03 = m[0][3];
		float m10 = m[1][0], m11 = m[1][1], m12 = m[1][2], m13 = m[1][3];
		float m20 = m[2][0], m21 = m[2][1], m22 = m[2][2], m23 = m[2][3];
		float m30 = m[3][0], m31 = m[3][1], m32 = m[3][2], m33 = m[3][3];

		float v0 = m20 * m31 - m21 * m30;
		float v1 = m20 * m32 - m22 * m30;
		float v2 = m20 * m33 - m23 * m30;
		float v3 = m21 * m32 - m22 * m31;
		float v4 = m21 * m33 - m23 * m31;
		float v5 = m22 * m33 - m23 * m32;

		float t00 = +(v5 * m11 - v4 * m12 + v3 * m13);
		float t10 = -(v5 * m10 - v2 * m12 + v1 * m13);
		float t20 = +(v4 * m10 - v2 * m11 + v0 * m13);
		float t30 = -(v3 * m10 - v1 * m11 + v0 * m12);

		float invDet = 1 / (t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03);

		float d00 = t00 * invDet;
		float d10 = t10 * invDet;
		float d20 = t20 * invDet;
		float d30 = t30 * invDet;

		float d01 = -(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
		float d11 = +(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
		float d21 = -(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
		float d31 = +(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

		v0 = m10 * m31 - m11 * m30;
		v1 = m10 * m32 - m12 * m30;
		v2 = m10 * m33 - m13 * m30;
		v3 = m11 * m32 - m12 * m31;
		v4 = m11 * m33 - m13 * m31;
		v5 = m12 * m33 - m13 * m32;

		float d02 = +(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
		float d12 = -(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
		float d22 = +(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
		float d32 = -(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

		v0 = m21 * m10 - m20 * m11;
		v1 = m22 * m10 - m20 * m12;
		v2 = m23 * m10 - m20 * m13;
		v3 = m22 * m11 - m21 * m12;
		v4 = m23 * m11 - m21 * m13;
		v5 = m23 * m12 - m22 * m13;

		float d03 = -(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
		float d13 = +(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
		float d23 = -(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
		float d33 = +(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

		return Matrix44(
			d00, d01, d02, d03,
			d10, d11, d12, d13,
			d20, d21, d22, d23,
			d30, d31, d32, d33);
	}

	void Matrix44::setTrans(const Vec3 &v)
	{
		m[0][3] = v.x;
		m[1][3] = v.y;
		m[2][3] = v.z;
	}

	void Matrix44::makeTrans(const Vec3 &v)
	{
		m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0; m[0][3] = v.x;
		m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0; m[1][3] = v.y;
		m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0; m[2][3] = v.z;
		m[3][0] = 0.0; m[3][1] = 0.0; m[3][2] = 0.0; m[3][3] = 1.0;
	}

	void Matrix44::makeTrans(float tx, float ty, float tz)
	{
		m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0; m[0][3] = tx;
		m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0; m[1][3] = ty;
		m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0; m[2][3] = tz;
		m[3][0] = 0.0; m[3][1] = 0.0; m[3][2] = 0.0; m[3][3] = 1.0;
	}

	Vec3 Matrix44::getTrans()
	{
		return Vec3(m[0][0], m[1][1], m[2][2]);
	}

	Matrix44 Matrix44::getTrans(const Vec3 &v)
	{
		Matrix44 t;
		t.m[0][0] = 1.0;
		t.m[1][1] = 1.0;
		t.m[2][2] = 1.0;
		t.m[3][3] = 1.0;
		t.m[0][3] = v.x;
		t.m[1][3] = v.y;
		t.m[2][3] = v.z;
		
		return t;
	}

	Matrix44 Matrix44::getTrans(float tx, float ty, float tz)
	{
		Matrix44 t;
		t.m[0][0] = 1.0;
		t.m[1][1] = 1.0;
		t.m[2][2] = 1.0;
		t.m[3][3] = 1.0;
		t.m[0][3] = tx;
		t.m[1][3] = ty;
		t.m[2][3] = tz;

		return t;
	}

	void Matrix44::setScale(const Vec3 &v)
	{
		m[0][0] = v.x;
		m[1][1] = v.y;
		m[2][2] = v.z;
	}

	Matrix44 Matrix44::getScale(const Vec3 &v)
	{
		Matrix44 r;
		r.m[0][0] = v.x;
		r.m[1][1] = v.y;
		r.m[2][2] = v.z;
		r.m[3][3] = 1.0;

		return r;
	}

	Matrix44 Matrix44::getScale(float tx, float ty, float tz)
	{
		Matrix44 r;
		r.m[0][0] = tx;
		r.m[1][1] = ty;
		r.m[2][2] = tz;
		r.m[3][3] = 1.0;

		return r;
	}

	bool Matrix44::operator==(const Matrix44 &mat) const
	{
		for (size_t row = 0; row < 4; ++row)
		{
			for (size_t col = 0; col < 4; ++col)
			{
				if (m[row][col] != mat.m[row][col]) return false;
			}
		}

		return true;
	}

}