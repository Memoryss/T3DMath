#ifndef __VECTOR_H__
#define __VECTOR_H__

namespace T3D {

	template<class T>
	class Vec2
	{
	public:
		T x, y;

		T & operator[] (int index);
		const T & operator[] (int index) const;  //供常量对象使用

		Vec2(); //默认
		explicit Vec2(T a);  //(a, a)
		Vec2(T a, T b); //(a, b)

		Vec2(const Vec2 &v);
		template<class S> Vec2(const Vec2<S> &v);
		const Vec2 & operator=(const Vec2 &v);

		template<class S>
		void setValue(S a, S b);

		template<class S>
		void setValue(const Vec2<S> &v);

		template<class S>
		void getValue(S &a, S &b);

		template<class S>
		void getValue(Vec2<S> &v);

		T * getValue();
		const  T * getValue() const;

		template<class>
	};

}

#endif