/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [Main] <|                                                 *
 *                                                                             *
 * MPI/OpenMP-parallel gravitational collisional N-body simulator via modern   *
 * C++ and hypermodularity.                                                    *
 *                                                                             *
 * (Maybe) updated version and other (related) cool stuff:                     *
 * https://github.com/emaballarin/distnbody-hm                                 *
 *                                                                             *
 *                                                                             *
 * Copyright (C) 2020-* Emanuele Ballarin <emanuele@ballarin.cc>.              *
 * All rights reserved. Distribution: Apache License 2.0.                      *
 * Full license text: https://ballarin.cc/legal/licenses/apacheii.txt          *
 *                                                                             *
 *                                                                             *
 * "Any fool can write code that a computer can understand. Good programmers   *
 *  write code that humans can understand."                                    *
 *                                               (M. Fowlers)                  *
 *                                                                             *
 * ─────────────────────────────────────────────────────────────────────────── */


/******************
 * GLOBAL PRAGMAS *
 ******************/

#pragma once


/*******************
 * GLOBAL INCLUDES *
 *******************/

#include <cassert>
#include <iostream>
#include <mpi.h>
//#include <openmpi/ompi/mpi/cxx/mpicxx.h>


/******************
 * LOCAL INCLUDES *
 ******************/

#include "immutables.hpp"
#include "physics.hpp"
#include "precision.hpp"
#include "triparticle.hpp"
#include "trivector.hpp"
#include "utilityfx.hpp"


/******************************************************************************
 ******************************************************************************/


/********
 * MAIN *
 ********/

int main(int argc, char* argv[])
{

    /* Define required MPI variables */
    int procId;
    int errCode;
    int P;
    MPI_Status status;
    int master{0};

    /* Initialize MPI */
    errCode = MPI_Init(&argc, &argv);
    errCode = MPI_Comm_size(MPI_COMM_WORLD, &P);
    errCode = MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    /*
     * ~~ BRIEFING: ~~
     *
     * From now on:
     * - P is nodes total;
     * - procId is the node index.
     *
     * This means that:
     * - Loops iterate in the range [0, P) or [0, P-1]
     * - procId is the procId+1-th node.
     */





    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    std::cout << procId << " of " << P << "\n\n";



    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    MPI_Finalize();
    return 0;
}
