/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [TriVector] <|                                            *
 *                                                                             *
 * An efficient and user-friendly structure for 3D vectors.                    *
 *                                                                             *
 * (Maybe) updated version and other (related) cool stuff:                     *
 * https://github.com/emaballarin/distnbody-hm                                 *
 *                                                                             *
 *                                                                             *
 * Copyright (C) 2020-* Emanuele Ballarin <emanuele@ballarin.cc>.              *
 * All rights reserved. Distribution: Apache License 2.0.                      *
 * Full license text: https://ballarin.cc/legal/licenses/apacheii.txt          *
 *                                                                             *
 * ─────────────────────────────────────────────────────────────────────────── */


/******************
 * GLOBAL PRAGMAS *
 ******************/

#pragma once


/************
 * INCLUDES *
 ************/

#include <cmath>          // sqrt
#include <iostream>       // Overloading of `<<`
#include "precision.hpp"  // Handling of real numbers
#include "utilityfx.hpp"  // Utility functions


/******************************************************************************
 ******************************************************************************/

class TriVec {

    public:

    inline TriVec(real_t _x, real_t _y, real_t _z) : x(_x), y(_y), z(_z)
    {
    }

    inline TriVec() : x{}, y{}, z{}
    {
    }

    inline TriVec operator-() const
    {
        TriVec v{-x, -y, -z};
        return v;
    }

    inline TriVec& operator+=(const TriVec& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    inline TriVec& operator-=(const TriVec& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    template<typename Scalar>
    inline TriVec& operator*=(const Scalar& scalar)
    {
        x * scalar;
        y * scalar;
        z * scalar;
        return *this;
    }

    template<typename Scalar>
    inline TriVec& operator/=(const Scalar& scalar)
    {
        x * (static_cast<real_t>(1) / scalar);
        y * (static_cast<real_t>(1) / scalar);
        z * (static_cast<real_t>(1) / scalar);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const TriVec vec)
    {
        os << "[ "
           << vec.x << " " << vec.y << " " << vec.z
           << " ]";
        return os;
    }

    // Impose box boundary conditions, in-place
    inline void boxBC(real_t boxx, real_t boxy, real_t boxz)
    {
        assert(boxx >= 0);
        assert(boxy >= 0);
        assert(boxz >= 0);

        x = rangeloop<long int>(x, boxx);
        y = rangeloop<long int>(y, boxy);
        z = rangeloop<long int>(z, boxz);
    }

    // Impose CUBIC box boundary conditions, in-place
    inline void boxBC_cubic(real_t boxside)
    {
        assert(boxside >= 0);

        x = rangeloop<long int>(x, boxside);
        y = rangeloop<long int>(y, boxside);
        z = rangeloop<long int>(z, boxside);
    }

    real_t x;
    real_t y;
    real_t z;

};


inline TriVec operator+(const TriVec& lhs, const TriVec& rhs)
{
    return TriVec(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

inline TriVec operator-(const TriVec& lhs, const TriVec& rhs)
{
    return TriVec(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}


template<typename Scalar>
inline TriVec operator*(const Scalar& scalar, const TriVec& vec)
{
    return TriVec(scalar*(vec.x), scalar*(vec.y), scalar*(vec.z));
}

template<typename Scalar>
inline TriVec operator/(const TriVec& vec, const Scalar& scalar)
{
    return (static_cast<real_t>(1)/scalar)*vec;
}

template<typename Scalar>
inline TriVec operator*(const TriVec& vec, const Scalar& scalar)
{
    return scalar * vec;
}

inline real_t dot(const TriVec& lhs, const TriVec& rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

inline real_t mod(const TriVec& vec)
{
    return static_cast<real_t>(sqrt(dot(vec, vec)));
}
