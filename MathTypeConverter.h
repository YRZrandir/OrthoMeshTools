#ifndef MATH_TYPE_CONVERTER_H
#define MATH_TYPE_CONVERTER_H
#define MATH_TYPE_CONVERTER_EIGEN
#define MATH_TYPE_CONVERTER_CGAL
#define MATH_TYPE_CONVERTER_ASSIMP

#ifdef MATH_TYPE_CONVERTER_GLM
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#endif
#ifdef MATH_TYPE_CONVERTER_EIGEN
#include <Eigen/Eigen>
#endif
#ifdef MATH_TYPE_CONVERTER_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Direction_3.h>
#endif
#ifdef MATH_TYPE_CONVERTER_ASSIMP
#include <assimp/vector2.h>
#include <assimp/vector3.h>
#include <assimp/matrix3x3.h>
#include <assimp/matrix4x4.h>
#include <assimp/quaternion.h>
#endif

//////////////////////////////////////////////
//TO GLM
//////////////////////////////////////////////
#ifdef MATH_TYPE_CONVERTER_GLM
#ifdef MATH_TYPE_CONVERTER_EIGEN
/*
* Eigen => GLM
*/
template<typename T, int R, int C> inline
std::enable_if_t<(R > 1 && C > 1), glm::mat<C, R, T>>
ToGLM( const Eigen::Matrix<T, R, C>& m )
{
    glm::mat<C, R, T> result;
    for (int i = 0; i < R; i++)
        for (int j = 0; j < C; j++)
            result[j][i] = m( i, j );
    return result;
}

template<typename T, int L> inline
glm::vec<L, T>
ToGLM( Eigen::Matrix<T, L, 1> v )
{
    glm::vec<L, T> result;
    for (int i = 0; i < L; i++)
        result[i] = v( i );
    return result;
}

template <typename T, int L> inline
glm::vec<L, T>
ToGLM( Eigen::Matrix<T, 1, L> v )
{
    glm::vec<L, T> result;
    for (int i = 0; i < L; i++)
        result[i] = v( i );
    return result;
}

template <typename M> inline
std::enable_if_t<
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M>&& M::ColsAtCompileTime == 1 && M::RowsAtCompileTime != 1,
    glm::vec<M::RowsAtCompileTime, typename M::Scalar>
>
ToGLM( const M& m )
{
    glm::vec<M::RowsAtCompileTime, typename M::Scalar> result;
    for (int i = 0; i < M::RowsAtCompileTime; i++)
    {
        result[i] = m( i );
    }
    return result;
}

template <typename M> inline
std::enable_if_t<
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M>&& M::RowsAtCompileTime == 1,
    glm::vec<M::ColsAtCompileTime, typename M::Scalar>
>
ToGLM( const M& m )
{
    glm::vec<M::ColsAtCompileTime, typename M::Scalar> result;
    for (int i = 0; i < M::ColsAtCompileTime; i++)
    {
        result[i] = m( i );
    }
    return result;
}

template <typename T> inline
glm::qua<T>
ToGLM( Eigen::Quaternion<T> q )
{
    return glm::qua<T>( q.w(), q.x(), q.y(), q.z() );
}
#endif
#ifdef MATH_TYPE_CONVERTER_CGAL
/*
* CGAL => GLM
*/
template <typename Kernel> inline
glm::vec<2, typename Kernel::FT>
ToGLM( CGAL::Point_2<Kernel> p )
{
    return glm::vec<2, typename Kernel::FT>( p.x(), p.y() );
}

template <typename Kernel> inline
glm::vec<3, typename Kernel::FT>
ToGLM( CGAL::Point_3<Kernel> p )
{
    return glm::vec<3, typename Kernel::FT>( p.x(), p.y(), p.z() );
}

template <typename Kernel> inline
glm::vec<2, typename Kernel::FT>
ToGLM( CGAL::Vector_2<Kernel> v )
{
    return glm::vec<2, typename Kernel::FT>( v.x(), v.y() );
}

template <typename Kernel> inline
glm::vec<3, typename Kernel::FT>
ToGLM( CGAL::Vector_3<Kernel> v )
{
    return glm::vec<3, typename Kernel::FT>( v.x(), v.y(), v.z() );
}

template <typename Kernel> inline
glm::vec<2, typename Kernel::FT>
ToGLM( CGAL::Direction_2<Kernel> d )
{
    return ToGLM( d.vector() );
}

template <typename Kernel> inline
glm::vec<3, typename Kernel::FT>
ToGLM( CGAL::Direction_3<Kernel> d )
{
    return ToGLM( d.vector() );
}
#endif
#ifdef MATH_TYPE_CONVERTER_ASSIMP
template <typename T> inline
glm::vec<2, T>
ToGLM( const aiVector2t<T>& v )
{
    return glm::vec<2, T>( v.x, v.y );
}

template <typename T> inline
glm::vec<3, T>
ToGLM( const aiVector3t<T>& v )
{
    return glm::vec<3, T>( v.x, v.y, v.z );
}

template <typename T> inline
glm::mat<3, 3, T>
ToGLM( const aiMatrix3x3t<T>& m )
{
    return glm::mat<3, 3, T>(
        m.a1, m.b1, m.c1,
        m.a2, m.b2, m.c2,
        m.a3, m.b3, m.c3 );
}

template<typename T> inline
glm::mat<4, 4, T>
ToGLM( const aiMatrix4x4t<T>& m )
{
    return glm::mat<4, 4, T>(
        m.a1, m.b1, m.c1, m.d1,
        m.a2, m.b2, m.c2, m.d2,
        m.a3, m.b3, m.c3, m.d3,
        m.a4, m.b4, m.c4, m.d4 );
}

template<typename T> inline
glm::qua<T>
ToGLM( const aiQuaterniont<T>& q )
{
    return glm::qua<T>( q.w, q.x, q.y, q.z );
}
#endif
#endif

//////////////////////////////////////////////
//TO EIGEN
//////////////////////////////////////////////
#ifdef MATH_TYPE_CONVERTER_EIGEN
#ifdef MATH_TYPE_CONVERTER_GLM
template<typename T, int R, int C> inline
Eigen::Matrix<T, R, C> ToEigen( const glm::mat<C, R, T>& m )
{
    Eigen::Matrix<T, R, C> result;
    for (int i = 0; i < R; i++)
        for (int j = 0; j < C; j++)
            result( i, j ) = m[j][i];
    return result;
}

template<typename T, int L> inline
Eigen::Matrix<T, L, 1> ToEigen( glm::vec<L, T> v )
{
    Eigen::Matrix<T, L, 1> result;
    for (int i = 0; i < L; i++)
        result( i ) = v[i];
    return result;
}

template <typename T> inline
Eigen::Quaternion<T> ToEigen( glm::qua<T> q )
{
    return Eigen::Quaternion<T>( q.w, q.x, q.y, q.z );
}
#endif
#ifdef MATH_TYPE_CONVERTER_CGAL
template <typename Kernel> inline
Eigen::Matrix<typename Kernel::FT, 2, 1>
ToEigen( CGAL::Point_2<Kernel> p )
{
    return Eigen::Matrix<typename Kernel::FT, 2, 1>( p.x(), p.y() );
}

template <typename Kernel> inline
Eigen::Matrix<typename Kernel::FT, 3, 1>
ToEigen( CGAL::Point_3<Kernel> p )
{
    return Eigen::Matrix<typename Kernel::FT, 3, 1>( p.x(), p.y(), p.z() );
}

template <typename Kernel> inline
Eigen::Matrix<typename Kernel::FT, 2, 1>
ToEigen( CGAL::Vector_2<Kernel> v )
{
    return Eigen::Matrix<typename Kernel::FT, 2, 1>( v.x(), v.y() );
}

template<typename Kernel> inline
Eigen::Matrix<typename Kernel::FT, 3, 1>
ToEigen( CGAL::Vector_3<Kernel> v )
{
    return Eigen::Matrix<typename Kernel::FT, 3, 1>( v.x(), v.y(), v.z() );
}

template <typename Kernel> inline
Eigen::Matrix<typename Kernel::FT, 2, 1>
ToEigen( CGAL::Direction_2<Kernel> d )
{
    return ToEigen( d.vector() );
}

template <typename Kernel> inline
Eigen::Matrix<typename Kernel::FT, 3, 1>
ToEigen( CGAL::Direction_3<Kernel> d )
{
    return ToEigen( d.vector() );
}
#endif
#ifdef MATH_TYPE_CONVERTER_ASSIMP
template<typename T> inline
Eigen::Matrix<T, 2, 1>
ToEigen( const aiVector2t<T>& v )
{
    return Eigen::Matrix<T, 2, 1>( v.x, v.y );
}

template<typename T> inline
Eigen::Matrix<T, 3, 1>
ToEigen( const aiVector3t<T>& v )
{
    return Eigen::Matrix<T, 3, 1>( v.x, v.y, v.z );
}

template<typename T> inline
Eigen::Matrix<T, 3, 3>
ToEigen( const aiMatrix3x3t<T>& m )
{
    Eigen::Matrix<T, 3, 3> result;
    result( 0, 0 ) = m.a1;
    result( 0, 1 ) = m.a2;
    result( 0, 2 ) = m.a3;
    result( 1, 0 ) = m.b1;
    result( 1, 1 ) = m.b2;
    result( 1, 2 ) = m.b3;
    result( 2, 0 ) = m.c1;
    result( 2, 1 ) = m.c2;
    result( 2, 2 ) = m.c3;
    return result;
}

template<typename T> inline
Eigen::Matrix<T, 4, 4>
ToEigen( const aiMatrix4x4t<T>& m )
{
    Eigen::Matrix<T, 3, 3> result;
    result( 0, 0 ) = m.a1;
    result( 0, 1 ) = m.a2;
    result( 0, 2 ) = m.a3;
    result( 0, 3 ) = m.a4;
    result( 1, 0 ) = m.b1;
    result( 1, 1 ) = m.b2;
    result( 1, 2 ) = m.b3;
    result( 1, 3 ) = m.b4;
    result( 2, 0 ) = m.c1;
    result( 2, 1 ) = m.c2;
    result( 2, 2 ) = m.c3;
    result( 2, 3 ) = m.c4;
    result( 3, 0 ) = m.d1;
    result( 3, 1 ) = m.d2;
    result( 3, 2 ) = m.d3;
    result( 3, 3 ) = m.d4;
    return result;
}

template<typename T> inline
Eigen::Quaternion<T>
ToEigen( const aiQuaterniont<T>& q )
{
    return Eigen::Quaternion<T>( q.w, q.x, q.y, q.z );
}
#endif
#endif

//////////////////////////////////////////////
//TO CGAL
//////////////////////////////////////////////
#ifdef MATH_TYPE_CONVERTER_CGAL
#ifdef MATH_TYPE_CONVERTER_GLM
template <typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_2<Kernel>
ToCGAL( glm::vec<2, T> v )
{
    return CGAL::Vector_2<Kernel>( v.x, v.y );
}

template <typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_3<Kernel>
ToCGAL( glm::vec<3, T> v )
{
    return CGAL::Vector_3<Kernel>( v.x, v.y, v.z );
}

#endif
#ifdef MATH_TYPE_CONVERTER_EIGEN
template<typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_2<Kernel>
ToCGAL( const Eigen::Matrix<T, 2, 1>& v )
{
    return CGAL::Vector_2<Kernel>( v[0], v[1] );
}

template<typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_3<Kernel>
ToCGAL( const Eigen::Matrix<T, 3, 1>& v )
{
    return CGAL::Vector_3<Kernel>( v[0], v[1], v[2] );
}

template<typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_2<Kernel>
ToCGAL( const Eigen::Matrix<T, 1, 2>& v )
{
    return CGAL::Vector_2<Kernel>( v[0], v[1] );
}

template<typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_3<Kernel>
ToCGAL( const Eigen::Matrix<T, 1, 2>& v )
{
    return CGAL::Vector_3<Kernel>( v[0], v[1], v[2] );
}

template<typename M, typename Kernel = CGAL::Simple_cartesian<typename M::Scalar>> inline
std::enable_if_t <
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M> && (M::ColsAtCompileTime == 2 && M::RowsAtCompileTime == 1 || M::ColsAtCompileTime == 1 && M::RowsAtCompileTime == 2),
    CGAL::Vector_2<Kernel>
>
ToCGAL( const M& m )
{
    return CGAL::Vector_2<Kernel>( m[0], m[1] );
}

template<typename M, typename Kernel = CGAL::Simple_cartesian<typename M::Scalar>> inline
std::enable_if_t <
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M> && (M::ColsAtCompileTime == 3 && M::RowsAtCompileTime == 1 || M::ColsAtCompileTime == 1 && M::RowsAtCompileTime == 3),
    CGAL::Vector_3<Kernel>
>
ToCGAL( const M& m )
{
    return CGAL::Vector_3 <Kernel>( m[0], m[1], m[2] );
}

#endif
#ifdef MATH_TYPE_CONVERTER_ASSIMP
template<typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_2<Kernel>
ToCGAL( const aiVector2t<T>& v )
{
    return CGAL::Vector_2<Kernel>( v.x, v.y );
}

template<typename T, typename Kernel = CGAL::Simple_cartesian<T>> inline
CGAL::Vector_3<Kernel>
ToCGAL( const aiVector3t<T>& v )
{
    return CGAL::Vector_3<Kernel>( v.x, v.y, v.z );
}
#endif
#endif

//////////////////////////////////////////////
//TO ASSIMP
//////////////////////////////////////////////
#ifdef MATH_TYPE_CONVERTER_ASSIMP
#ifdef MATH_TYPE_CONVERTER_GLM
template<typename T> inline
aiVector2t<T>
ToAssimp( const glm::vec<2, T>& v )
{
    return aiVector2t<T>( v.x, v.y );
}

template<typename T> inline
aiVector3t<T>
ToAssimp( const glm::vec<3, T>& v )
{
    return aiVector3t<T>( v.x, v.y, v.z );
}

template<typename T> inline
aiMatrix3x3t<T>
ToAssimp( const glm::mat<3, 3, T>& m )
{
    return aiMatrix3x3t<T>( m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2] );
}

template<typename T> inline
aiMatrix4x4t<T>
ToAssimp( const glm::mat<4, 4, T>& m )
{
    return aiMatrix4x4t<T>(
        m[0][0], m[1][0], m[2][0], m[3][0],
        m[0][1], m[1][1], m[2][1], m[3][1],
        m[0][2], m[1][2], m[2][2], m[3][2],
        m[0][3], m[1][3], m[2][3], m[3][3] );
}

template<typename T> inline
aiQuaterniont<T>
ToAssimp( const glm::qua<T>& q )
{
    return aiQuaterniont<T>( q.w, q.x, q.y, q.z );
}
#endif
#ifdef MATH_TYPE_CONVERTER_EIGEN
template<typename T> inline
std::enable_if_t<std::is_arithmetic_v<T>, aiVector2t<T>>
ToAssimp( const Eigen::Matrix<T, 2, 1>& v )
{
    return aiVector2t<T>( v[0], v[1] );
}

template<typename T> inline
std::enable_if_t<std::is_arithmetic_v<T>, aiVector3t<T>>
ToAssimp( const Eigen::Matrix<T, 3, 1>& v )
{
    return aiVector3t<T>( v[0], v[1], v[2] );
}

template<typename T> inline
std::enable_if_t<std::is_arithmetic_v<T>, aiVector2t<T>>
ToAssimp( const Eigen::Matrix<T, 1, 2>& v )
{
    return aiVector2t<T>( v[0], v[1] );
}

template<typename T> inline
std::enable_if_t<std::is_arithmetic_v<T>, aiVector3t<T>>
ToAssimp( const Eigen::Matrix<T, 1, 3>& v )
{
    return aiVector3t<T>( v[0], v[1], v[2] );
}

template<typename T> inline
std::enable_if_t<std::is_arithmetic_v<T>, aiMatrix3x3t<T>>
ToAssimp( const Eigen::Matrix<T, 3, 3>& m )
{
    return aiMatrix3x3t<T>(
        m( 0, 0 ), m( 0, 1 ), m( 0, 2 ),
        m( 1, 0 ), m( 1, 1 ), m( 1, 2 ),
        m( 2, 0 ), m( 2, 1 ), m( 2, 2 ) );
}

template<typename T> inline
std::enable_if_t<std::is_arithmetic_v<T>, aiMatrix4x4t<T>>
ToAssimp( const Eigen::Matrix<T, 4, 4>& m )
{
    return aiMatrix4x4t<T>(
        m( 0, 0 ), m( 0, 1 ), m( 0, 2 ), m( 0, 3 ),
        m( 1, 0 ), m( 1, 1 ), m( 1, 2 ), m( 1, 3 ),
        m( 2, 0 ), m( 2, 1 ), m( 2, 2 ), m( 2, 3 ),
        m( 3, 0 ), m( 3, 1 ), m( 3, 2 ), m( 3, 3 ) );
}

template<typename M> inline
std::enable_if_t<
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M> && (M::ColsAtCompileTime == 2 && M::RowsAtCompileTime == 1 || M::ColsAtCompileTime == 1 && M::RowsAtCompileTime == 2),
    aiVector2t<typename M::Scalar>
>
ToAssimp( const M& m )
{
    return aiVector2t<typename M::Scalar>( m( 0 ), m( 1 ) );
}

template<typename M> inline
std::enable_if_t<
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M> && (M::ColsAtCompileTime == 3 && M::RowsAtCompileTime == 1 || M::ColsAtCompileTime == 1 && M::RowsAtCompileTime == 3),
    aiVector3t<typename M::Scalar>
>
ToAssimp( const M& m )
{
    return aiVector3t<typename M::Scalar>( m( 0 ), m( 1 ), m( 2 ) );
}

template<typename M> inline
std::enable_if_t<
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M>&& M::ColsAtCompileTime == 3 && M::RowsAtCompileTime == 3,
    aiMatrix3x3t<typename M::Scalar>
>
ToAssimp( const M& m )
{
    return aiMatrix3x3t<typename M::Scalar>(
        m( 0, 0 ), m( 0, 1 ), m( 0, 2 ),
        m( 1, 0 ), m( 1, 1 ), m( 1, 2 ),
        m( 2, 0 ), m( 2, 1 ), m( 2, 2 ) );
}

template<typename M> inline
std::enable_if_t<
    std::is_base_of_v<Eigen::DenseCoeffsBase<M, Eigen::ReadOnlyAccessors>, M>&& M::ColsAtCompileTime == 4 && M::RowsAtCompileTime == 4,
    aiMatrix4x4t<typename M::Scalar>
>
ToAssimp( const M& m )
{
    return aiMatrix4x4t<typename M::Scalar>(
        m( 0, 0 ), m( 0, 1 ), m( 0, 2 ), m( 0, 3 ),
        m( 1, 0 ), m( 1, 1 ), m( 1, 2 ), m( 1, 3 ),
        m( 2, 0 ), m( 2, 1 ), m( 2, 2 ), m( 2, 3 ),
        m( 3, 0 ), m( 3, 1 ), m( 3, 2 ), m( 3, 3 )
        );
}

template<typename T> inline
aiQuaterniont<T>
ToAssimp( const Eigen::Quaternion<T>& q )
{
    return aiQuaterniont<T>( q.w(), q.x(), q.y(), q.z() );
}


#endif
#ifdef MATH_TYPE_CONVERTER_CGAL
template<typename Kernel>
aiVector2t<typename Kernel::FT>
ToAssimp( const CGAL::Vector_2<Kernel>& v )
{
    return aiVector2t<typename Kernel::FT>( v.x(), v.y() );
}

template<typename Kernel>
aiVector3t<typename Kernel::FT>
ToAssimp( const CGAL::Vector_3<Kernel>& v )
{
    return aiVector3t<typename Kernel::FT>( v.x(), v.y(), v.z() );
}

template<typename Kernel>
aiVector2t<typename Kernel::FT>
ToAssimp( const CGAL::Point_2<Kernel>& p )
{
    return aiVector2t<typename Kernel::FT>( p.x(), p.y() );
}

template<typename Kernel>
aiVector3t<typename Kernel::FT>
ToAssimp( const CGAL::Point_3<Kernel>& p )
{
    return aiVector3t<typename Kernel::FT>( p.x(), p.y(), p.z() );
}

template<typename Kernel>
aiVector2t<typename Kernel::FT>
ToAssimp( const CGAL::Direction_2<Kernel>& d )
{
    return aiVector2t<typename Kernel::FT>( d.dx(), d.dy() );
}

template<typename Kernel>
aiVector3t<typename Kernel::FT>
ToAssimp( const CGAL::Direction_3<Kernel>& d )
{
    return aiVector3t<typename Kernel::FT>( d.dx(), d.dy(), d.dz() );
}
#endif
#endif

#endif