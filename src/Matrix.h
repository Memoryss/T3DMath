#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "Vector.h"

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

		bool operator==(const Matrix33 &mat);

		Matrix33 operator+(const Matrix33 &mat);

		Matrix33 operator*(const Matrix33 &mat);

	protected:
		float m[3][3];


	};

}

#endif
