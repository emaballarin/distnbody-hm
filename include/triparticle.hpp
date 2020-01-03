/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [TriParticle] <|                                          *
 *                                                                             *
 * An efficient and user-friendly structure for 3D material points             *
 * (a.k.a. particles)                                                          *
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

#include <cassert>
#include "precision.hpp"  // Handling of real numbers
#include "utilityfx.hpp"  // Utility functions
#include "trivector.hpp"  // 3D vectors


/******************************************************************************
 ******************************************************************************/


class TriPart
{
    public:

    /********
    * CTORS *
    *********/

    inline TriPart(real_t _m) : position{}, velocity{}, energy{}, mass(_m)
    {
        assert(_m > 0); // Mass must be positive
    }

    inline TriPart(TriVec _pos, real_t _m) : position(_pos), velocity{}, energy{}, mass(_m)
    {
        assert(_m > 0);  // Mass must be positive
    }

    inline TriPart(TriVec _pos, TriVec _vel, real_t _m) : position(_pos), velocity(_vel), energy{}, mass(_m)
    {
        assert(_m > 0);  // Mass must be positive
    }

    // NOTE: Eneegy initialization MUST be done manually after allocation! (otw: E=0)

    /***********
     * GETTERS *
     ***********/

    inline TriVec pos()
    {
        return position;
    }

    inline real_t x()
    {
        return position.x;
    }

    inline real_t y()
    {
        return position.y;
    }

    inline real_t z()
    {
        return position.z;
    }

    inline TriVec vel()
    {
        return velocity;
    }

    inline real_t vx()
    {
        return velocity.x;
    }

    inline real_t vy()
    {
        return velocity.y;
    }

    inline real_t vz()
    {
        return velocity.z;
    }

    inline real_t E()
    {
        return energy;
    }

    inline real_t m()
    {
        return mass;
    }


    /***********
     * SETTERS *
     ***********/

    inline void pos(TriVec _pos)
    {
        position.x = _pos.x;
        position.y = _pos.y;
        position.z = _pos.z;
    }

    inline void pos(real_t _pos[3]) // Do not const-ify!
    {
        /* Do not assign null pointer dereferences */
        assert(_pos != nullptr);
        assert(_pos + 1 != nullptr);
        assert(_pos + 2 != nullptr);

        position.x = _pos[0];
        position.y = _pos[1];
        position.z = _pos[2];
    }

    inline void x(real_t _x)
    {
        position.x = _x;
    }

    inline void y(real_t _y)
    {
        position.y = _y;
    }

    inline void z(real_t _z)
    {
        position.z = _z;
    }

    inline void vel(TriVec _vel)
    {
        position.x = _vel.x;
        position.y = _vel.y;
        position.z = _vel.z;
    }

    inline void vel(real_t _vel[3])  // Do not const-ify!
    {
        /* Do not assign null pointer dereferences */
        assert(_vel != nullptr);
        assert(_vel + 1 != nullptr);
        assert(_vel + 2 != nullptr);

        position.x = _vel[0];
        position.y = _vel[1];
        position.z = _vel[2];
    }

    inline void vx(real_t _vx)
    {
        velocity.x = _vx;
    }

    inline void vy(real_t _vy)
    {
        velocity.y = _vy;
    }

    inline void vz(real_t _vz)
    {
        velocity.z = _vz;
    }

    inline void E(real_t _E)
    {
        assert(_E >= 0);
        energy = _E;
    }

   /**************
    * MEMBER FXS *
    **************/

    inline real_t Ekin()
   {
        return (mass*dot(velocity, velocity))/2;
   }

    private:
    TriVec position;
    TriVec velocity;
    real_t energy;
    const real_t mass;
};
