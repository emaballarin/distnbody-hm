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
#include <vector>
#include <bits/stdc++.h>
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


/***********
 * DEFINES *
 ***********/
#define li long int
#define lli long long int
#define ui unsigned int
#define ul unsigned long int
#define ull unsigned long long int
#define SEED 35791246  // Seed for number generator


/******************************************************************************
 ******************************************************************************/

/*********************
 * UTILITY FUNCTIONS *
 *********************/

 inline real_t drand48_vel_adjust(real_t uniform0_1)
{
     return (uniform0_1*0.10) - 0.05;
}

//inline RealTwople rangesplit(real_t _rangestart, real_t _rangestop, ul _nSlices, ul _mySlice_zero_indexed)
//{
//     assert(_rangestop >= _rangestart);
//     const real_t _slicelen = (_rangestop -  _rangestart)/(_nSlices);
//     RealTwople _toberet{_rangestart + _mySlice_zero_indexed * _slicelen, 0};
//     _toberet.stop = _toberet.start + _slicelen;
//     return _toberet;
//}
//
//inline bool A_lies_in_BC (real_t A, real_t B, real_t C)
//{
//     assert(B > C);
//     return (A>=B && A<C);
//}


/********
 * MAIN *
 ********/

int main(int argc, char* argv[])
{

    /* Define global, problem-specific constants */
    constexpr real_t GlobalMass = 100;
    constexpr ull Nparticles = 100;  // Change this if required!

    /* Compute global, problem-specific constants */
    constexpr real_t singlemass = GlobalMass / Nparticles;

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
     * - Loops iterate in the range [0, P) or [0, P-1] --> USE THIS: for (auto i{0}; i < wrknum; ++i)
     * - procId is the procId+1-th node.
     */

    /* GENERATE THE INITIAL CONDITIONS (in a smart way!) */

    // How many we need? (per node)
    ull particles_to_gen = workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(procId)).stop
                           - workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(procId)).start;

    // Seed the RNG
    srand48(SEED * (procId + 1));

    // Generate position/velocity of particles in arrays

    auto genXpos = new real_t[particles_to_gen];
    auto genXvel = new real_t[particles_to_gen];
    auto genYpos = new real_t[particles_to_gen];
    auto genYvel = new real_t[particles_to_gen];
    auto genZpos = new real_t[particles_to_gen];
    auto genZvel = new real_t[particles_to_gen];

    for (ull i{0}; i < particles_to_gen; ++i)
    {
        genXpos[i] = drand48();  // Already in the [0, 1) range!
        genYpos[i] = drand48();  // Already in the [0, 1) range!
        genZpos[i] = drand48();  // Already in the [0, 1) range!
        genXvel[i] = drand48_vel_adjust(drand48());
        genYvel[i] = drand48_vel_adjust(drand48());
        genZvel[i] = drand48_vel_adjust(drand48());
    }

    /*
     * ~~ And now: ready for a trick? ~~
     */
    std::sort(genXpos, genXpos + particles_to_gen);
    /*
     * We choose to partition the work to be done according to the X coordinate. Seems legit.
     *
     * However, we suppose that all the remaining numbers (Y, Z positions, X, Y, Z velocities)
     * have been generated in a shuffled order. The shuffling follows the min-to-max order of
     * X positions.
     *
     * This does not alter the random generation, as it could be performed by seeding adequately
     * the drand48() generator after having generated fully `genXpos` and right before generating
     * the remaining vectors. Since this pre-seeding has to be performed by adding or removing
     * known integer constants to/from the previous seed (as shown by George Marsaglia), this
     * produces a distribution indistinguishable from a genuine one output by drand48() via
     * single initial seeding (Cook, 2011).
     *
     * We therefore "just" sort in-place `genXpos` and assign to particles according to the i-th
     * position in the array(s).
     *
     * This "trick" has enormous positive performance impact on the following development...
     * So efficient!
     */


    /* DECOMPOSE THE DOMAIN (again, in a smart way!) */

    // Instead of using logical conditions, we assign to MPI nodes this way...
    auto assignments = new int[particles_to_gen];
    // Split loop for further (eventual) OpenMP parallelization!
    for (ull i{0}; i < particles_to_gen; ++i)
    {
        assignments[i] = static_cast<int>((genXpos[i] / (1.0 / P)));  // `1` is the cube side length!
    }

    // Of all the generated positions, determine which node will handle them. We can use ranges!
    // Must be sequential!
    auto genAssignments = new IntegerTwople<ull>[P];
    int curthread = 0;
    genAssignments[curthread].start = 0;        // Start of first
    for (ull i{0}; i < particles_to_gen; ++i)
    {
        if (assignments[i] > curthread)
        {
            genAssignments[curthread].stop = i;     // End of previous
            curthread++;
            genAssignments[curthread].start = i;    // Start of new
        }
    }
    genAssignments[curthread].stop = particles_to_gen;  // End of last


    // Now we need to know the amount of particles that given node has generated

    ull generated_by_node[P];   // Requires C >= 99
    for (int i{0}; i < P; ++i)
    {
        if (i == procId)
        {
            generated_by_node[i] = particles_to_gen;
        }
        else
        {
            generated_by_node[i] = workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(i)).stop
                                   - workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(i)).start;
        }
    }


    // We now need to receive the positions from all other nodes (will be needed for physical evolution anyway!)
    // We need also to share them in a matrix, to give order to everything!

    /*
     * ~~ Another trick, right-off the 'trick hat'... ~~
     */
    std::vector<std::vector<real_t>> PositionMatrix(P);
    for (int i = 0; i < P; i++)
    {
        PositionMatrix[i].resize(generated_by_node[i]);
    }
    /*
     * In a different context, this should be considered to be idiotic. However, in this case...
     *
     * Here we need to build a "thing" that has:
     * - 2 dimensions;
     * - Both dimensions known only at runtime;
     * - Preferrably one dimension has to dependent on the other (i.e. a list of non-equidimensional arrays)
     * - Preferrably each row has to be castable to and from raw arrays.
     *
     * std::vectors are slower than arrays (not too much, though!) and this seems the most reasonable
     * alternative of boost::multiarrays, std::lists of boost::arrays of other similar results
     * of types mumbo-jumbo.
     *
     * The alternative is to write a custom class in order to that. It is surely the best thing to do.
     * However, it will take time, a lot of testing and a sparkle of genial attitude to build
     * something more battle-tested than C++'s std:: library.
     *
     * Probably fair enough.
     */

    // To be continued...





    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    // Testing area...

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    MPI_Finalize();
    return 0;
}
