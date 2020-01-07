/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [Physics] <|                                              *
 *                                                                             *
 * Allowing the physics to happen.                                             *
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

#include <cassert>        // Diagnostics
#include <cmath>          // Powers
#include "precision.hpp"  // Handling of real numbers
#include "immutables.hpp" // Physical constants
#include "trivector.hpp"  // 3D vectors


/***********
 * DEFINES *
 ***********/

//#define Ggrav PhysicalGgrav
#define Ggrav FakeGgrav

/******************************************************************************
 ******************************************************************************/


/* Simulator constants */

constexpr real_t veleps = 0.05; // Maximum allowed relative variation on velocity (modulus)

/*******************************************************************************/


/* Utility functions */

inline TriVec accel_from_Fm(TriVec _force, real_t _mass)
{
    assert(_mass >= 0);
    return _force/_mass;
}

/*******************************************************************************/


/* Dynamics Integrators */

inline TriVec DeltaVel(TriVec _accel, real_t _timestep)
{
    assert(_timestep > 0);
    return _accel*_timestep;
}

inline TriVec NewVel(TriVec _oldvel, TriVec _accel, real_t _timestep)
{
    assert(_timestep > 0);
    return _oldvel + DeltaVel(_accel, _timestep);
}

inline TriVec DeltaPos(TriVec _vel, real_t _timestep)
{
    assert(_timestep > 0);
    return _vel * _timestep;
}

inline TriVec NewPos(TriVec _oldpos, TriVec _vel, real_t _timestep)
{
    assert(_timestep > 0);
    return _oldpos + DeltaPos(_vel, _timestep);
}

inline TriVec DeltaPos_single_pass(TriVec _oldvel, TriVec _accel, real_t _timestep)
{
    assert(_timestep > 0);
    return NewVel(_oldvel, _accel, _timestep)*_timestep;
}

inline TriVec NewPos_single_pass(TriVec _oldpos, TriVec _oldvel, TriVec _accel, real_t _timestep)
{
    assert(_timestep > 0);
    return _oldpos + DeltaPos_single_pass(_oldvel, _accel, _timestep);
}

/*******************************************************************************/

/* Gravitation */

// Acceptable warnings: we just want to avoid +inf forces/energies

inline TriVec Fgrav_on1_from2(real_t _m1, real_t _m2, TriVec _pos1, TriVec _pos2)
{
    assert((_pos1.x != _pos2.x) || (_pos1.y != _pos2.y) || (_pos1.z != _pos2.z));
    TriVec _dist(_pos1 - _pos2);
    return (Ggrav*_m1*_m2)*_dist/(pow(mod(_dist), 3));
}

inline real_t Egrav_of1_from2(real_t _m1, real_t _m2, TriVec _pos1, TriVec _pos2)
{
    assert((_pos1.x != _pos2.x) || (_pos1.y != _pos2.y) || (_pos1.z != _pos2.z));
    TriVec _dist(_pos1 - _pos2);
    return -(Ggrav * _m1 * _m2) / mod(_dist);
}
