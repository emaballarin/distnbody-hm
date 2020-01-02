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


/************
 * INCLUDES *
 ************/

#include "utilityfx.hpp"


/******************************************************************************
 ******************************************************************************/


/* Eye-candy to avoid many static_casts in simple computations with booleans */
template<typename integer>
inline integer oneif(bool condition)
{
    return static_cast<integer>(condition);
}

/* Templated implementation of Kronecker's Delta */
template<typename candidate, typename rettype>
inline rettype kdelta(candidate lhs, candidate rhs)
{
    return static_cast<rettype>(oneif<int>(lhs == rhs));
}

/* Loop over a modular arithmetic (loop over rings) */
template<typename integer>
inline integer ringloop(integer candidate, integer base)
{
    return candidate % base;
}
