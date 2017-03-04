set (PUB_HDR
	Matrix.h
	Vector.h
	Quaternion.h
	TMath.h
	Plane.h
	Frustum.h
)

add_library(T3DMath
	Matrix.cpp
	Vector.cpp
	Quaternion.cpp
	Plane.cpp
	Frustum.cpp
	${PUB_HDR}
)

set_target_properties (T3DMath PROPERTIES
  DEBUG_POSTFIX "_d"
  PUBLIC_HEADER "${PUB_HDR}"
  )

install (TARGETS T3DMath
	ARCHIVE DESTINATION lib
	PUBLIC_HEADER DESTINATION include/T3DMath)