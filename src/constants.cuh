/*---------------------------------------------------------------------------*\
|                                                                             |
| MULTIC-TS-LBM: CUDA-based multicomponent Lattice Boltzmann Method           |
| Developed at UDESC - State University of Santa Catarina                     |
| Website: https://www.udesc.br                                               |
| Github: https://github.com/brenogemelgo/MULTIC-TS-LBM                       |
|                                                                             |
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\

Copyright (C) 2023 UDESC Geoenergia Lab
Authors: Breno Gemelgo (Geoenergia Lab, UDESC)

License
    This file is part of MULTIC-TS-LBM.

    MULTIC-TS-LBM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

Description
    A header defining the constants used in the simulation

    Namespace

SourceFiles
    constants.cuh

\*---------------------------------------------------------------------------*/

#ifndef CONSTANTS_CUH
#define CONSTANTS_CUH

#include "cuda/utils.cuh"
#include "structs/LBMFields.cuh"
#include "functions/constexprFor.cuh"
#include "velocitySet/velocitySet.cuh"

namespace LBM
{
#if defined(VS_D3Q19)
    using VelocitySet = D3Q19;
#elif defined(VS_D3Q27)
    using VelocitySet = D3Q27;
#endif
}

namespace Phase
{
    using VelocitySet = LBM::D3Q7;
}

#define RUN_MODE
// #define SAMPLE_MODE
// #define PROFILE_MODE

#if defined(RUN_MODE)

static constexpr int MACRO_SAVE = 1000;
static constexpr int NSTEPS = 100000;

#elif defined(SAMPLE_MODE)

static constexpr int MACRO_SAVE = 100;
static constexpr int NSTEPS = 10000;

#elif defined(PROFILE_MODE)

static constexpr int MACRO_SAVE = 1;
static constexpr int NSTEPS = 0;

#endif

namespace mesh
{
    static constexpr label_t res = 128;
    static constexpr label_t nx = res;
    static constexpr label_t ny = res * 2;
    static constexpr label_t nz = res * 2;

    static constexpr int diam_water = 13;
    static constexpr int diam_oil = 13;

    static constexpr int radius_water = diam_water / 2;
    static constexpr int radius_oil = diam_oil / 2;
}

namespace physics
{
    static constexpr scalar_t u_water = static_cast<scalar_t>(0.06);
    static constexpr scalar_t u_oil = static_cast<scalar_t>(0.05);

    static constexpr int reynolds_water = 1400;
    static constexpr int reynolds_oil = 450;

    static constexpr scalar_t rho_water = static_cast<scalar_t>(1);
    static constexpr scalar_t rho_oil = static_cast<scalar_t>(0.82);

    static constexpr int weber = 500;
    static constexpr scalar_t sigma = rho_oil * (u_oil * u_oil * mesh::diam_oil) / weber;

    static constexpr scalar_t width = static_cast<scalar_t>(1);
    static constexpr scalar_t gamma = static_cast<scalar_t>(static_cast<double>(1) / static_cast<double>(width));
}

#endif