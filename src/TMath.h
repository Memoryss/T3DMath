#ifndef __TMATH_H__
#define __TMATH_H__

#include <math.h>

#define PI ((float)3.141592654f)
#define PI2 ((float)6.283185307f)
#define PI_DIV_2 ((float)1.570796327f)
#define PI_DIV_4 ((float)0.785398163f)
#define PI_INV ((float)0.318309886f)

#define DEG_TO_RAN(ang) ((ang) * PI / 180.0f)
#define RAD_TO_DEG(rads) ((rads) * 180.0f / PI)

typedef unsigned int uint32;
typedef unsigned char uint8;

namespace T3D {

	class Math
	{
	public:

		Math() { BuildSinCosTables(); }

		static void BuildSinCosTables();

		static float FastSin(float theta);
		static float FastCos(float theta);

		static float Arcsin(float value);
		static float Arccos(float value);

		//像素转化 将24位无alpha转为uint32 brga存储  (bmp brg)
		static uint32 uint8BGR_to_uint32BGRA(const uint8 * const ch);

		//像素比例相乘
		static uint32 ColorScaleFloat(uint32 color, float scale);

		//将计算好的sin cos保存起来
		static float cos_look[361];
		static float sin_look[361];
	};

	void Math::BuildSinCosTables()
	{
		for (size_t i = 0; i < 361; ++i)
		{
			sin_look[i] = sinf(DEG_TO_RAN(i));
			cos_look[i] = cosf(DEG_TO_RAN(i));
		}
	}

	float Math::FastCos(float theta)
	{
		theta = fmodf(theta, 360);
		if (theta < 0) theta += 360;

		int theta_int = (int)theta;
		float theta_frac = theta - theta_int;

		return cos_look[theta_int] + theta_frac * (cos_look[theta_int + 1] - cos_look[theta_int]);
	}

	float Math::FastSin(float theta)
	{
		theta = fmodf(theta, 360);
		if (theta < 0) theta += 360;

		int theta_int = (int)theta;
		float theta_frac = theta - theta_int;

		return sin_look[theta_int] + theta_frac * (sin_look[theta_int + 1] - sin_look[theta_int]);
	}

	float Math::Arccos(float value)
	{
		if (-1.0 < value)
		{
			if (value < 1.0)
				return RAD_TO_DEG(acos(value));
			else
				return 0.0;
		}
		else
		{
			return RAD_TO_DEG(PI);
		}
	}

	float Math::Arcsin(float fValue)
	{
		if (-1.0 < fValue)
		{
			if (fValue < 1.0)
				return RAD_TO_DEG(asin(fValue));
			else
				return 90;
		}
		else
		{
			return 270;
		}
	}

	uint32 Math::uint8BGR_to_uint32BGRA(const uint8 * const ch)
	{
		//BMP中 颜色排列是bgr
		uint32 texel = 0;
		texel += (*ch << 24 ) | (*(ch + 1) << 16) | (*(ch + 2) << 8) | 0xff;  //uint32 bgra
		return texel;
	}

	uint32 Math::ColorScaleFloat(uint32 color, float scale)
	{
		//TODO 可以做一下优化
		return (int)((color & 0xff000000 >> 24) * scale) << 24 |
			(int)((color & 0x00ff0000 >> 16) * scale) << 16 |
			(int)((color & 0x0000ff00 >> 8) * scale) << 8 |
			(int)((color & 0x000000ff >> 0) * scale) << 0;
	}
}

#endif
