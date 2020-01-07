/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [Utility Functions] <|                                    *
 *                                                                             *
 * Utility functions for various usage cases.                                  *
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

#include "precision.hpp"  // Handling of real numbers
#include <cassert>        // Diagnostics
#include <cmath>          // Absolute value


/******************************************************************************
 ******************************************************************************/

/* Eye-candy to avoid many static_casts in simple computations with booleans */
template<typename integer>
inline integer oneif(const bool& condition)
{
    return static_cast<integer>(condition);
}

/* Templated implementation of Kronecker's Delta */
template<typename candidate, typename rettype>
inline rettype kdelta(const candidate& lhs, const candidate& rhs)
{
    return static_cast<rettype>(oneif<int>(lhs == rhs));
}

/* Loop over a modular arithmetic (loop over rings) */
template<typename integer>
inline integer ringloop(const integer& candidate, const integer& base)
{
    /* While waiting for C++20 Concepts... */
    assert(base > 0);  // It will work in Z\0; We just don't do that here ;)
    // assert(base != 0)    // Will it ever be needed...

    const integer remainder{candidate % base};
    return (remainder < 0) ? (base + remainder) : remainder;
}

/* Loop over a range (loop over ranges) */
template<typename integer>
inline real_t rangeloop(const real_t& candidate, const real_t& base)
{
    /* While waiting for C++20 Concepts... */
    assert(base > 0);

    const real_t remainder{candidate - base*static_cast<integer>(candidate/base)};
    return (remainder < 0) ? (base + remainder) : remainder;
}

/* Real two-ple (integer struct for two values) */
struct RealTwople
{
    real_t start;
    real_t stop;
};

/* Integer two-ple (integer struct for two values) */
template<typename integer>
struct IntegerTwople
{
    integer start;
    integer stop;
};

/* Start-Stop partitioner (workload slicer) */
template<typename integer>
inline IntegerTwople<integer> workslice(const integer& workload, const integer& slices, const integer& worker)
{
    /* While waiting for C++20 Concepts... */
    assert(workload >= 0);    // Cannot share a "negative workload"
    assert(slices > 0);       // Cannot divide by zero workers; Cannot have a negative number of workers
    assert(worker >= 0);      // A consequence of the above
    assert(worker < slices);  // A consequence of the above

    /*
     * We adopt here the 'minimal unbalance criterion' in order to assign to
     * each worker the right amount-of-work to be performed.
     * Correct use in for loops is: [start, stop)
     */

    const integer slicefloor{static_cast<integer>(workload / slices)};
    const integer remainder{workload - (slices * slicefloor)};
    const bool OneMoreNeeded{worker < remainder};

    const integer slicelen{slicefloor + oneif<integer>(OneMoreNeeded)};

    IntegerTwople<integer> ToBeRet{worker * slicelen + remainder * oneif<integer>(!OneMoreNeeded), 0};
    ToBeRet.stop = ToBeRet.start + slicelen;

    return ToBeRet;
}
