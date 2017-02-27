#include "Matrix.h"

#include <string.h>
#include <math.h>

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
		return m[row];
	}

	const float * Matrix33::operator[](int row) const
	{
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
	}

}