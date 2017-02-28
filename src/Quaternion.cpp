#include "Quaternion.h"

#include <math.h>
#include <algorithm>
#include <assert.h>

#include "TMath.h"

namespace T3D {



	Quaternion::Quaternion() : w(0), x(0), y(0), z(0)
	{

	}

	Quaternion::Quaternion(float fw, float fx, float fy, float fz) : w(fw), x(fx), y(fy), z(fz)
	{

	}

	Quaternion::Quaternion(const Matrix33 &rot)
	{
		this->FromRotationMatrix(rot);
	}

	Quaternion::Quaternion(const float degree, const Vec3 &v)
	{
		this->FromAngleAxis(degree, v);
	}

	Quaternion::Quaternion(const Vec3 &xaxis, const Vec3 &yaxis, const Vec3 &zaxis)
	{
		this->FromAxis(xaxis, yaxis, zaxis);
	}

	void Quaternion::swap(Quaternion &quat)
	{
		std::swap(quat.w, w);
		std::swap(quat.x, x);
		std::swap(quat.y, y);
		std::swap(quat.z, z);
	}

	float Quaternion::operator[](const size_t index) const
	{
		assert(index < 4);
		return (&w)[index];
	}

	float & Quaternion::operator[](const size_t index)
	{
		assert(index < 4);
		return (&w)[index];
	}

	void Quaternion::FromRotationMatrix(const Matrix33 &kRot)
	{
		// Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
		// article "Quaternion Calculus and Fast Animation".

		float fTrace = kRot[0][0] + kRot[1][1] + kRot[2][2];
		float fRoot;

		if (fTrace > 0.0)
		{
			// |w| > 1/2, may as well choose w > 1/2
			fRoot = sqrt(fTrace + 1.0f);  // 2w
			w = 0.5f*fRoot;
			fRoot = 0.5f / fRoot;  // 1/(4w)
			x = (kRot[2][1] - kRot[1][2])*fRoot;
			y = (kRot[0][2] - kRot[2][0])*fRoot;
			z = (kRot[1][0] - kRot[0][1])*fRoot;
		}
		else
		{
			// |w| <= 1/2
			static size_t s_iNext[3] = { 1, 2, 0 };
			size_t i = 0;
			if (kRot[1][1] > kRot[0][0])
				i = 1;
			if (kRot[2][2] > kRot[i][i])
				i = 2;
			size_t j = s_iNext[i];
			size_t k = s_iNext[j];

			fRoot = sqrt(kRot[i][i] - kRot[j][j] - kRot[k][k] + 1.0f);
			float* apkQuat[3] = { &x, &y, &z };
			*apkQuat[i] = 0.5f*fRoot;
			fRoot = 0.5f / fRoot;
			w = (kRot[k][j] - kRot[j][k])*fRoot;
			*apkQuat[j] = (kRot[j][i] + kRot[i][j])*fRoot;
			*apkQuat[k] = (kRot[k][i] + kRot[i][k])*fRoot;
		}
	}

	void Quaternion::ToRotationMatrix(Matrix33 &kRot)
	{
		float fTx = x + x;
		float fTy = y + y;
		float fTz = z + z;
		float fTwx = fTx*w;
		float fTwy = fTy*w;
		float fTwz = fTz*w;
		float fTxx = fTx*x;
		float fTxy = fTy*x;
		float fTxz = fTz*x;
		float fTyy = fTy*y;
		float fTyz = fTz*y;
		float fTzz = fTz*z;

		kRot[0][0] = 1.0f - (fTyy + fTzz);
		kRot[0][1] = fTxy - fTwz;
		kRot[0][2] = fTxz + fTwy;
		kRot[1][0] = fTxy + fTwz;
		kRot[1][1] = 1.0f - (fTxx + fTzz);
		kRot[1][2] = fTyz - fTwx;
		kRot[2][0] = fTxz - fTwy;
		kRot[2][1] = fTyz + fTwx;
		kRot[2][2] = 1.0f - (fTxx + fTyy);
	}

	void Quaternion::FromAngleAxis(const float degree, const Vec3 &rAxis)
	{
		w = Math::FastCos(degree * 0.5);
		float sinDg = Math::FastSin(degree * 0.5);
		x = sinDg * rAxis.x;
		y = sinDg * rAxis.y;
		z = sinDg * rAxis.z;
	}

	void Quaternion::ToAngleAxis(float &degree, Vec3 &rAxis)
	{
		float value = Math::Arccos(w);
		degree = value * 2;

		float fsin_inv = 1 / Math::FastSin(value);
		rAxis.x = x * fsin_inv;
		rAxis.y = y * fsin_inv;
		rAxis.z = z * fsin_inv
	}

	void Quaternion::FromAxes(const Vec3 &xaxis, const Vec3 &yaxis, const Vec3 &zaxis)
	{
		Matrix33 kRot;

		kRot[0][0] = xaxis.x;
		kRot[1][0] = xaxis.y;
		kRot[2][0] = xaxis.z;

		kRot[0][1] = yaxis.x;
		kRot[1][1] = yaxis.y;
		kRot[2][1] = yaxis.z;

		kRot[0][2] = zaxis.x;
		kRot[1][2] = zaxis.y;
		kRot[2][2] = zaxis.z;

		FromRotationMatrix(kRot);
	}

	void Quaternion::ToAxes(Vec3 &xAxis, Vec3 &yAxis, Vec3 &zAxis)
	{
		Matrix33 rot;
		ToRotationMatrix(rot);

		xAxis[0] = rot[0][0];
		xAxis[1] = rot[1][0];
		xAxis[2] = rot[2][0];

		yAxis[0] = rot[0][1];
		yAxis[1] = rot[1][1];
		yAxis[2] = rot[2][1];

		zAxis[0] = rot[0][2];
		zAxis[1] = rot[1][2];
		zAxis[2] = rot[2][2];
	}

	Quaternion Quaternion::operator+(const Quaternion &quat) const
	{
		return Quaternion(w + quat.w, x + quat.x, y + quat.y, z + quat.z);
	}

	Quaternion Quaternion::operator-(const Quaternion &quat) const
	{
		return Quaternion(w - quat.w, x - quat.x, y - quat.y, z - quat.z);
	}

	Quaternion Quaternion::operator-() const
	{
		return Quaternion(-w, -x, -y, -z);
	}

	float Quaternion::Dot(const Quaternion &quat) const
	{
		return w * quat.w + x * quat.x + y * quat.y + z * quat.z;
	}

	float Quaternion::Length() const
	{
		return w * w + x * x + y * y + z * z;
	}

	float Quaternion::Normalize()
	{
		float len = Length();
		float factor = 1 / sqrt(len);
		*this = *this * factor;
		return len;
	}

	Quaternion Quaternion::Inverse() const
	{
		float fNorm = w*w + x*x + y*y + z*z;
		if (fNorm > 0.0)
		{
			float fInvNorm = 1.0f / fNorm;
			return Quaternion(w*fInvNorm, -x*fInvNorm, -y*fInvNorm, -z*fInvNorm);
		}
		else
		{
			return Quaternion(0, 0, 0, 0);
		}
	}

	bool Quaternion::operator==(const Quaternion &quat) const
	{
		return w == quat.w && x == quat.x && y == quat.y && z == quat.z;
	}

	Quaternion Quaternion::operator*(const Quaternion &quat) const
	{
		return Quaternion(
			w * quat.w - x * quat.x - y * quat.y - z * quat.z,
			w * quat.x + x * quat.w + y * quat.z - z * quat.y,
			w * quat.y + y * quat.w + z * quat.x - x * quat.z,
			w * quat.z + z * quat.w + x * quat.y - y * quat.x
		);
	}

	Quaternion Quaternion::operator*(const float scale) const
	{
		return Quaternion(w * scale, x * scale, y * scale, z * scale);
	}

	Vec3 Quaternion::operator*(const Vec3 &v) const
	{
		Vec3 uv, uuv;
		Vec3 qvec(x, y, z);
		uv = qvec.Cross(v);
		uuv = qvec.Cross(uv);

		uv = uv * (2.0f * w);
		uuv = uuv * 2.0f;

		return v + uv + uuv;
	}

	float Quaternion::getRoll() const
	{
		//TODO
	}

	float Quaternion::getPitch() const
	{
		//TODO
	}

	float Quaternion::getYaw() const
	{
		//TODO
	}

	Quaternion & Quaternion::operator=(const Quaternion &quat)
	{
		w = quat.w;
		x = quat.x;
		y = quat.y;
		z = quat.z;
	}

}