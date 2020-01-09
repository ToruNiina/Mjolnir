#ifndef MJOLNIR_CUDA_VECTOR_HPP
#define MJOLNIR_CUDA_VECTOR_HPP
#include <mjolnir/util/macro.hpp>
#include <mjolnir/math/vector_util.hpp>
#include <mjolnir/math/functions.hpp>
#include <vector_types.h>
#include <vector_functions.h>
#include <iosfwd>

// ---------------------------------------------------------------------------
// use double4 and float4 as a coordinate_type.
//
// Here, some operators and (meta-) functions are defined to use them in the
// same way as mjolnir::math::Vector.
//
// The following stuff is defined.
//
// - meta function to convert float/double and float4/double4
//   - mjolnir::cuda::real_type_of<T>
//   - mjolnir::cuda::vector_type_of<T>
// - bridges between mjolnir::math and float4/double4
//   - mjolnir::math::make_coordinate<float4|double4>
//   - mjolnir::math::X, Y, Z(float4|double4)
//   - mjolnir::math::operator<<(float4|double4)
// - math functions for float4|double4
//   - mjolnir::math::dot_product(float4|double4)
//   - mjolnir::math::cross_product(float4|double4)
//   - mjolnir::math::length_sq(float4|double4)
//   - mjolnir::math::length(float4|double4)
//   - mjolnir::math::rlength(float4|double4)
// - math operators for float4|double4
//   - mjolnir::cuda::operator+(float4|double4) [unary]
//   - mjolnir::cuda::operator-(float4|double4) [unary]
//   - mjolnir::cuda::operator+(float4|double4, float4|double4) [binary]
//   - mjolnir::cuda::operator-(float4|double4, float4|double4) [binary]
//   - mjolnir::cuda::operator*(float4|double4, float |double)
//   - mjolnir::cuda::operator*(float |double,  float4|double4)
//   - mjolnir::cuda::operator/(float4|double4, float |double)

namespace mjolnir
{

namespace cuda
{
template<typename T>
struct vector_type_of;
template<>
struct vector_type_of<double> {using type = double4;};
template<>
struct vector_type_of<float>  {using type = float4;};

template<typename T>
using vector_type_of_t = typename vector_type_of<T>::type;
} // cuda

// ---------------------------------------------------------------------------
// mjolnjir::math bridges

namespace math
{
// ---------------------------------------------------------------------------
// utility meta-functions to convert float4 -> float

template<>
struct real_type_of<double4> {using type = double;};
template<>
struct real_type_of<float4>  {using type = float;};

// ---------------------------------------------------------------------------
// utility functions to make (T, T, T) -> vecT

template<>
struct make_coordinate_impl<double4>
{
    MJOLNIR_CUDA_HOST_DEVICE static double4
    invoke(const double x, const double y, const double z) noexcept
    {
        return ::make_double4(x, y, z, 0.0);
    }
};
template<>
struct make_coordinate_impl<float4>
{
    MJOLNIR_CUDA_HOST_DEVICE static float4
    invoke(const float x, const float y, const float z) noexcept
    {
        return ::make_float4(x, y, z, 0.0f);
    }
};

template<typename realT>
MJOLNIR_CUDA_DEVICE_ONLY
cuda::vector_type_of_t<realT>
make_vec(const realT x, const realT y, const realT z) noexcept
{
    return make_coordinate_impl<cuda::vector_type_of_t<realT>>::invoke(x, y, z);
}

MJOLNIR_CUDA_HOST_DEVICE inline double  X(const double4& v) noexcept {return v.x;}
MJOLNIR_CUDA_HOST_DEVICE inline double& X(double4&       v) noexcept {return v.x;}
MJOLNIR_CUDA_HOST_DEVICE inline double  Y(const double4& v) noexcept {return v.y;}
MJOLNIR_CUDA_HOST_DEVICE inline double& Y(double4&       v) noexcept {return v.y;}
MJOLNIR_CUDA_HOST_DEVICE inline double  Z(const double4& v) noexcept {return v.z;}
MJOLNIR_CUDA_HOST_DEVICE inline double& Z(double4&       v) noexcept {return v.z;}

MJOLNIR_CUDA_HOST_DEVICE inline float   X(const float4& v)  noexcept {return v.x;}
MJOLNIR_CUDA_HOST_DEVICE inline float&  X(float4&       v)  noexcept {return v.x;}
MJOLNIR_CUDA_HOST_DEVICE inline float   Y(const float4& v)  noexcept {return v.y;}
MJOLNIR_CUDA_HOST_DEVICE inline float&  Y(float4&       v)  noexcept {return v.y;}
MJOLNIR_CUDA_HOST_DEVICE inline float   Z(const float4& v)  noexcept {return v.z;}
MJOLNIR_CUDA_HOST_DEVICE inline float&  Z(float4&       v)  noexcept {return v.z;}

// ---------------------------------------------------------------------------
// math functions
//
// XXX Note for the namespace:
// here the following cuda operators are deinfed under the mjolnir:: namepsace,
// not `mjolnir::math::`. The reason is the following. `float4` and `double4`
// are defined in the global namespace because they are builtin structs for
// cuda-C architecture. Thus argument dependent lookup (ADL) does not work and
// compiler cannot find the operators unless we explicitly write
// `using (mjolnir::)math::operator+` if we wrap the following functions by
// `mjolnir::math` namespace.

MJOLNIR_CUDA_HOST_DEVICE
inline double dot_product(const double4& lhs, const double4& rhs) noexcept
{
    return X(lhs) * X(rhs) + Y(lhs) * Y(rhs) + Z(lhs) * Z(rhs);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float dot_product(const float4& lhs, const float4& rhs) noexcept
{
    return X(lhs) * X(rhs) + Y(lhs) * Y(rhs) + Z(lhs) * Z(rhs);
}

MJOLNIR_CUDA_HOST_DEVICE
inline double4 cross_product(const double4& lhs, const double4& rhs) noexcept
{
    return ::make_double4(Y(lhs) * Z(rhs) - Z(lhs) * Y(rhs),
                          Z(lhs) * X(rhs) - X(lhs) * Z(rhs),
                          X(lhs) * Y(rhs) - Y(lhs) * X(rhs), 0.0);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4 cross_product(const float4& lhs, const float4& rhs) noexcept
{
    return ::make_float4(Y(lhs) * Z(rhs) - Z(lhs) * Y(rhs),
                         Z(lhs) * X(rhs) - X(lhs) * Z(rhs),
                         X(lhs) * Y(rhs) - Y(lhs) * X(rhs), 0.0f);
}

MJOLNIR_CUDA_HOST_DEVICE
inline double length_sq(const double4& lhs) noexcept
{
    return X(lhs) * X(lhs) + Y(lhs) * Y(lhs) + Z(lhs) * Z(lhs);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float  length_sq(const float4& lhs) noexcept
{
    return X(lhs) * X(lhs) + Y(lhs) * Y(lhs) + Z(lhs) * Z(lhs);
}

MJOLNIR_CUDA_HOST_DEVICE
inline double length(const double4& lhs) noexcept
{
#ifdef __CUDA_ARCH__
    return ::sqrt(length_sq(lhs));
#else
    return std::sqrt(length_sq(lhs));
#endif
}
MJOLNIR_CUDA_HOST_DEVICE
inline float  length(const float4& lhs) noexcept
{
#ifdef __CUDA_ARCH__
    return ::sqrtf(length_sq(lhs));
#else
    return std::sqrt(length_sq(lhs));
#endif
}

MJOLNIR_CUDA_HOST_DEVICE
inline double rlength(const double4& lhs) noexcept
{
#ifdef __CUDA_ARCH__
    return ::rsqrt(length_sq(lhs));
#else
    return ::mjolnir::math::rsqrt(length_sq(lhs));
#endif
}
MJOLNIR_CUDA_HOST_DEVICE
inline float rlength(const float4& lhs) noexcept
{
#ifdef __CUDA_ARCH__
    return ::rsqrtf(length_sq(lhs));
#else
    return ::mjolnir::math::rsqrt(length_sq(lhs));
#endif
}

} // math

// ===========================================================================
// output operators for logger

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const double4& vec)
{
    os << '(' << math::X(vec) << ',' << math::Y(vec) << ',' << math::Z(vec) << ')';
    return os;
}
template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const float4& vec)
{
    os << '(' << math::X(vec) << ',' << math::Y(vec) << ',' << math::Z(vec) << ')';
    return os;
}

// ===========================================================================
// math operators

// ---------------------------------------------------------------------------
// unary operator + and -

MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator+(const double4& v) noexcept
{
    return v;
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4 operator+(const float4& v) noexcept
{
    return v;
}
MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator-(const double4& v) noexcept
{
    return make_double4(-v.x, -v.y, -v.z, -v.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4 operator-(const float4& v) noexcept
{
    return make_float4(-v.x, -v.y, -v.z, -v.w);
}

// ---------------------------------------------------------------------------
// binary operator +, -, * and /.

MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator+(const double4& lhs, const double4& rhs) noexcept
{
    return make_double4(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z, lhs.w+rhs.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4  operator+(const float4& lhs, const float4& rhs) noexcept
{
    return make_float4(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z, lhs.w+rhs.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline double4& operator+=(double4& lhs, const double4& rhs) noexcept
{
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    lhs.w += rhs.w;
    return lhs;
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4&  operator+=(float4& lhs, const float4& rhs) noexcept
{
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    lhs.w += rhs.w;
    return lhs;
}

MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator-(const double4& lhs, const double4& rhs) noexcept
{
    return make_double4(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z, lhs.w-rhs.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4  operator-(const float4& lhs, const float4& rhs) noexcept
{
    return make_float4(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z, lhs.w-rhs.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline double4& operator-=(double4& lhs, const double4& rhs) noexcept
{
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    lhs.w -= rhs.w;
    return lhs;
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4&  operator-=(float4& lhs, const float4& rhs) noexcept
{
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    lhs.w -= rhs.w;
    return lhs;
}

MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator*(const double4& lhs, const double rhs) noexcept
{
    return make_double4(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs, lhs.w*rhs);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4  operator*(const float4& lhs, const float rhs) noexcept
{
    return make_float4(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs, lhs.w*rhs);
}
MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator*(const double lhs, const double4& rhs) noexcept
{
    return make_double4(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z, lhs*rhs.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4  operator*(const float lhs, const float4& rhs) noexcept
{
    return make_float4(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z, lhs*rhs.w);
}
MJOLNIR_CUDA_HOST_DEVICE
inline double4& operator*=(double4& lhs, const double rhs) noexcept
{
    lhs.x *= rhs;
    lhs.y *= rhs;
    lhs.z *= rhs;
    lhs.w *= rhs;
    return lhs;
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4&  operator*=(float4& lhs, const float rhs) noexcept
{
    lhs.x *= rhs;
    lhs.y *= rhs;
    lhs.z *= rhs;
    lhs.w *= rhs;
    return lhs;
}

MJOLNIR_CUDA_HOST_DEVICE
inline double4 operator/(const double4& lhs, const double rhs) noexcept
{
    return make_double4(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs, lhs.w/rhs);
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4  operator/(const float4& lhs, const float rhs) noexcept
{
    return make_float4(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs, lhs.w/rhs);
}
MJOLNIR_CUDA_HOST_DEVICE
inline double4& operator/=(double4& lhs, const double rhs) noexcept
{
    lhs.x /= rhs;
    lhs.y /= rhs;
    lhs.z /= rhs;
    lhs.w /= rhs;
    return lhs;
}
MJOLNIR_CUDA_HOST_DEVICE
inline float4&  operator/=(float4& lhs, const float rhs) noexcept
{
    lhs.x /= rhs;
    lhs.y /= rhs;
    lhs.z /= rhs;
    lhs.w /= rhs;
    return lhs;
}

} // mjolnir
#endif // MJOLNIR_CUDA_VECTOR_HPP
