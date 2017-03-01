#include "Vector.h"

#include <math.h>
#include <assert.h>

namespace T3D {

	/************************************************************************/
	/* vec2                                                                     */
	/************************************************************************/
	Vec2::Vec2() : x(0), y(0)
	{

	}

	Vec2::Vec2(float a) : x(a), y(a)
	{

	}

	Vec2::Vec2(float a, float b) : x(a), y(b)
	{

	}

	Vec2::Vec2(const Vec2 &v) : x(v.x), y(v.y)
	{

	}

	float & Vec2::operator[](int index)
	{
		assert(index < 2);
		return (&x)[index];
	}

	const float & Vec2::operator[](int index) const
	{
		assert(index < 2);
		return (&x)[index];
	}

	Vec2 & Vec2::operator=(const Vec2 &v)
	{
		x = v.x;
		y = v.y;
		return *this;
	} 
	
	void Vec2::setValue(const float a, const float b)
	{
		x = a;
		y = b;
	}

	void Vec2::setValue(const Vec2 &v)
	{
		x = v.x;
		y = v.y;
	}

	void Vec2::getValue(float &a, float &b)
	{
		a = x;
		b = y;
	}

	void Vec2::getValue(Vec2 &v)
	{
		v.x = x;
		v.y = y;
	}

	Vec2 * Vec2::getValue()
	{
		return this;
	}

	const Vec2 * Vec2::getValue() const
	{
		return this;
	}

	Vec2 Vec2::operator-() const
	{
		return Vec2(-x, -y);
	}

	Vec2 Vec2::operator-(const Vec2 &v) const
	{
		return Vec2(x - v.x, y - v.y);
	}

	Vec2 Vec2::operator+(const Vec2 &v) const
	{
		return Vec2(x + v.x, y + v.y);
	}

	Vec2 Vec2::operator*(const float a) const
	{
		return Vec2(x * a, y * a);
	}

	bool Vec2::operator==(const Vec2 &v) const
	{
		return (x == v.x && y == v.y);
	}

	float Vec2::Length() const
	{
		//TODO 可以使用优化算法开根号
		return (float)sqrt(x * x + y * y);
	}

	float Vec2::LengthSqr() const
	{
		return x * x + y * y;
	}

	float Vec2::Normalize()
	{
		float len = Length();
		float length_inv = 1 / len;
		x *= length_inv;
		y *= length_inv;

		return len;
	}

	float Vec2::Distance(const Vec2 &v)
	{
		return (*this - v).Length();
	}

	float Vec2::DistanceSqr(const Vec2 &v)
	{
		return (*this - v).LengthSqr();
	}

	float Vec2::Dot(const Vec2 &v)
	{
		return x * v.x + y * v.y;
	}


	/************************************************************************/
	/* vec3                                                                     */
	/************************************************************************/

	float & Vec3::operator[](int index)
	{
		assert(index < 3);
		return (&x)[index];
	}

	const float & Vec3::operator[](int index) const
	{
		assert(index < 3);
		return (&x)[index];
	}

	Vec3::Vec3() : x(0), y(0), z(0)
	{

	}

	Vec3::Vec3(float a) : x(a), y(a), z(a)
	{

	}

	Vec3::Vec3(float a, float b, float c) : x(a), y(b), z(c)
	{

	}

	Vec3::Vec3(const Vec3 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	void Vec3::setValue(const float a, const float b, const float c)
	{
		x = a;
		y = b;
		z = c;
	}

	void Vec3::setValue(const Vec3 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	void Vec3::getValue(float &a, float &b, float &c)
	{
		a = x;
		b = y;
		c = z;
	}

	void Vec3::getValue(Vec3 &v)
	{
		v.x = x;
		v.y = y;
		v.z = z;
	}

	Vec3 * Vec3::getValue()
	{
		return this;
	}

	const Vec3 * Vec3::getValue() const
	{
		return this;
	}

	Vec3 Vec3::operator-() const
	{
		return Vec3(-x, -y, -z);
	}

	Vec3 Vec3::operator-(const Vec3 &v) const
	{
		return Vec3(x - v.x, y - v.y, z - v.z);
	}

	Vec3 Vec3::operator+(const Vec3 &v) const
	{
		return Vec3(x + v.x, y + v.y, z + v.z);
	}

	Vec3 Vec3::operator*(const float a) const
	{
		return Vec3(x * a, y * a, z * a);
	}

	float Vec3::Length() const
	{
		return (float)sqrt(x * x + y * y + z * z);
	}

	float Vec3::LengthSqr() const
	{
		return x * x + y * y + z * z;
	}

	float Vec3::Distance(const Vec3 &v)
	{
		return (*this - v).Length();
	}

	float Vec3::DistanceSqr(const Vec3 &v)
	{
		return (*this - v).LengthSqr();
	}

	float Vec3::Normalize()
	{
		float len = Length();
		float length_inv = 1 / len;
		x *= length_inv;
		y *= length_inv;
		z *= length_inv;

		return len;
	}

	float Vec3::Dot(const Vec3 &v)
	{
		return x * v.x + y * v.y + z * v.z;
	}

	Vec3 Vec3::Cross(const Vec3 &v)
	{
		Vec3 v3;
		v3.x = y * v.z - z * v.y;
		v3.y = z * v.x - x * v.z;
		v3.z = x * v.y - y * v.x;
		return v3;
	}

	bool Vec3::operator==(const Vec3 &v) const
	{
		return (x == v.x && y == v.y && z == v.z);
	}

	Vec3 & Vec3::operator=(const Vec3 &v) 
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	float & Vec4::operator[](int index) 
	{
		assert(index < 4);
		return (&x)[index];
	}

	const float & Vec4::operator[](int index) const
	{
		assert(index < 4);
		return (&x)[index];
	}

	Vec4::Vec4() : x(0), y(0), z(0), w(0)
	{
		
	}

	Vec4::Vec4(float a, float b, float c, float d) : x(a), y(b), z(c), w(d)
	{

	}

	Vec4::Vec4(const Vec4 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = v.w;
	}

	void Vec4::setValue(const float a, const float b, const float c, const float d)
	{
		x = a;
		y = b;
		z = c;
		w = d;
	}

	void Vec4::setValue(const Vec4 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = v.w;
	}

	void Vec4::getValue(float &a, float &b, float &c, float &d)
	{
		a = x;
		b = y;
		c = z;
		d = w;
	}

	void Vec4::getValue(Vec4 &v)
	{
		v.x = x;
		v.y = y;
		v.z = z;
		v.w = w;
	}

	Vec4 * Vec4::getValue()
	{
		return this;
	}

	const Vec4 * Vec4::getValue() const
	{
		return this;
	}

	Vec4 Vec4::operator-() const
	{
		return Vec4(-x, -y, -z, -w);
	}

	Vec4 Vec4::operator-(const Vec4 &v) const
	{
		return Vec4(x - v.x, y - v.y, z - v.z, w - v.w);
	}

	Vec4 Vec4::operator+(const Vec4 &v) const
	{
		return Vec4(x + v.x, y + v.y, z + v.z, w + v.w);
	}

	Vec4 Vec4::operator*(const float a) const
	{
		return Vec4(x * a, y * a, z * a, w * a);
	}

	float Vec4::Dot(const Vec4 &v)
	{
		return x * v.x + y * v.y + z * v.z + w * v.w;
	}

	bool Vec4::operator==(const Vec4 &v) const
	{
		return (x == v.x && y == v.y && z == v.z && w == v.w);
	}

	Vec4 & Vec4::operator=(const Vec4 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = v.w;

		return *this;
	}

	Vec3 Vec4::XYZ() const 
	{
		return Vec3(x, y, z);
	}

	Vec2 Vec4::XY() const
	{
		return Vec2(x, y);
	}
}