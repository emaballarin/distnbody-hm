/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [Immutables] <|                                           *
 *                                                                             *
 * General constants, const-expressions and other stuff that is not expected   *
 * to change. Ever (maybe).                                                    *
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
#include <cmath>  // M_PI


/***********
 * DEFINES *
 ***********/

#define _USE_MATH_DEFINES  // Acceptable warning: ensures back-compatibility.


/******************************************************************************
 ******************************************************************************/

/*
 * `constexpr` is the way to go: less impact on memory, same optimizations as `const`.
 * The matter is tricky: just don't use `const`. Don't try to be "smarter" either!
 */

/* Pi */
constexpr real_t pi{static_cast<real_t>(M_PI)};  // Fix precision once; use everywhere!

/* 2Pi */
constexpr real_t twopi{static_cast<real_t>(2 * M_PI)};  // Useful if the precision of
                                                        // real_t is less than M_PI

/* Newtonian Gravitational constant */
constexpr real_t PhysicalGgrav{static_cast<real_t>(6.67430e-11)};  // m^3 kg^-1 s^-2
constexpr real_t FakeGgrav{static_cast<real_t>(10e-6)};  // m^3 kg^-1 s^-2
