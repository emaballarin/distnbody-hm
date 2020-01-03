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

#include <cmath>           // sqrt
#include <iostream>        // Overloading of `<<`
#include "precision.hpp"  // Handling of real numbers
#include "utilityfx.hpp"  // Utility functions


/******************************************************************************
 ******************************************************************************/

class TriVec {

    /* Here we want to be as intuitive as possible */
    public:

    inline TriVec(real_t _x, real_t _y, real_t _z) : storagearray{_x, _y, _z}
    {
    }

    inline TriVec() : storagearray{}
    {
    }

    inline TriVec operator-() const
    {
        TriVec v{-x, -y, -z};
        return v;
    }

    inline TriVec& operator+=(const TriVec& rhs)
    {
        storagearray[0] += rhs.storagearray[0];
        storagearray[1] += rhs.storagearray[1];
        storagearray[2] += rhs.storagearray[2];
        return *this;
    }

    inline TriVec& operator-=(const TriVec& rhs)
    {
        storagearray[0] -= rhs.storagearray[0];
        storagearray[1] -= rhs.storagearray[1];
        storagearray[2] -= rhs.storagearray[2];
        return *this;
    }

    template<typename Scalar>
    inline TriVec& operator*=(const Scalar& scalar)
    {
        storagearray[0] * scalar;
        storagearray[1] * scalar;
        storagearray[2] * scalar;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const TriVec vec)
    {
        os << "[ "
           << vec.x << " " << vec.y << " " << vec.z
           << " ]";
        return os;
    }


    real_t& x{storagearray[0]};
    real_t& y{storagearray[1]};
    real_t& z{storagearray[2]};


    private:
    /* But we also want to allow easy vectorization by the compiler */
    real_t storagearray[3];
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
