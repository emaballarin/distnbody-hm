/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [Precision] <|                                            *
 *                                                                             *
 * Uniform, project-wide management of real number precision.                  *
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


/***********
 * DEFINES *
 ***********/

#define realdefault double  // A preset

/********************
 * TYPE DEFINITIONS *
 ********************/

// Change precision of real numbers via a `-D...` at compile time...
#if defined(REALFLOAT)
typedef float real_t;
#elif defined(REALDOUBLE)
typedef double real_t;
#elif defined(REALLONG)
typedef long double real_t;
// ... or just rely on a sensible default.
#else
typedef realdefault real_t;
#endif
