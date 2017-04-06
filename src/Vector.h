#ifndef __VECTOR_H__
#define __VECTOR_H__

namespace T3D {

	class Vec2
	{
	public:
		float x;
		float y;

		float & operator[] (int index);
		const float & operator[] (int index) const;  //供常量对象使用

		Vec2(); //默认
		explicit Vec2(float a);  //(a, a)
		Vec2(float a, float b); //(a, b)

		Vec2(const Vec2 &v);
		Vec2 & operator=(const Vec2 &v);

		void setValue(const float a, const float b);

		void setValue(const Vec2 &v);

		void getValue(float &a, float &b);

		void getValue(Vec2 &v);

		Vec2 * getValue();
		const  Vec2 * getValue() const;

		Vec2 operator-() const;
		Vec2 operator-(const Vec2 &v) const;

		Vec2 operator+(const Vec2 &v) const;

		Vec2 operator*(const float a) const;

		bool operator==(const Vec2 &v) const;

		float Length() const;
		float LengthSqr() const;

		float Distance(const Vec2 &v);
		float DistanceSqr(const Vec2 &v);

		float Normalize();

		float Dot(const Vec2 &v);

		void Swap(Vec2 &v);

		static const Vec2 ZERO;
		static const Vec2 UNIT_X;
		static const Vec2 UNIT_Y;
		static const Vec2 NEGATIVE_UNIT_X;
		static const Vec2 NEGATIVE_UNIT_Y;
		static const Vec2 UNIT_SCALE;
	};

	class Vec3
	{
	public:
		float x;
		float y;
		float z;

		float & operator[] (int index);
		const float & operator[] (int index) const;  //供常量对象使用

		Vec3();
		explicit Vec3(float a);
		Vec3(float a, float b, float c);
		Vec3(const Vec3 &v);
		Vec3 & operator=(const Vec3 &v);

		void setValue(const float a, const float b, const float c);
		void setValue(const Vec3 &v);
		void getValue(float &a, float &b, float &c);

		void getValue(Vec3 &v);

		Vec3 * getValue();
		const  Vec3 * getValue() const;

		Vec3 operator-() const;
		Vec3 operator-(const Vec3 &v) const;

		Vec3 operator+(const Vec3 &v) const;

		Vec3 operator*(const float a) const;

		bool operator==(const Vec3 &v) const;

		float Length() const;
		float LengthSqr() const;

		float Distance(const Vec3 &v);
		float DistanceSqr(const Vec3 &v);

		float Normalize();

		float Dot(const Vec3 &v);
		
		Vec3 Cross(const Vec3 &v);

		void Swap(Vec3 &v);

		static const Vec3 ZERO;
		static const Vec3 UNIT_X;
		static const Vec3 UNIT_Y;
		static const Vec3 UNIT_Z;
		static const Vec3 NEGATIVE_UNIT_X;
		static const Vec3 NEGATIVE_UNIT_Y;
		static const Vec3 NEGATIVE_UNIT_Z;
		static const Vec3 UNIT_SCALE;
	};

	class Vec4
	{
	public:
		float x;
		float y;
		float z;
		float w;

		float & operator[] (int index);
		const float & operator[] (int index) const;  //供常量对象使用

		Vec4();
		Vec4(float a, float b, float c, float d);
		Vec4(const Vec4 &v);
		Vec4 & operator=(const Vec4 &v);

		void setValue(const float a, const float b, const float c, const float d);
		void setValue(const Vec4 &v);
		void getValue(float &a, float &b, float &c, float &d);

		void getValue(Vec4 &v);

		Vec4 * getValue();
		const  Vec4 * getValue() const;

		Vec4 operator-() const;
		Vec4 operator-(const Vec4 &v) const;

		Vec4 operator+(const Vec4 &v) const;

		Vec4 operator*(const float a) const;

		bool operator==(const Vec4 &v) const;

		float Dot(const Vec4 &v);

		Vec3 XYZ() const;
		Vec2 XY() const;

		void Swap(Vec4 &v);

		static const Vec4 ZERO;
	};

}

#endif