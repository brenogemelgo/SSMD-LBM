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
    A class applying boundary conditions

Namespace
    LBM

SourceFiles
    boundaryConditions.cuh

\*---------------------------------------------------------------------------*/

#ifndef BOUNDARYCONDITIONS_CUH
#define BOUNDARYCONDITIONS_CUH

namespace LBM
{
    class BoundaryConditions
    {
    public:
        __host__ __device__ [[nodiscard]] inline consteval BoundaryConditions(){};

        __device__ static inline void applyWaterInflow(LBMFields d) noexcept
        {
            const label_t x = threadIdx.x + blockIdx.x * blockDim.x;
            const label_t z = threadIdx.y + blockIdx.y * blockDim.y;

            if (x >= mesh::nx || z >= mesh::nz)
            {
                return;
            }

            const scalar_t dx = static_cast<scalar_t>(x) - geometry::center_x();
            const scalar_t dz = static_cast<scalar_t>(z) - geometry::z_pos();
            const scalar_t r2 = dx * dx + dz * dz;

            if (r2 > geometry::R2_water())
            {
                return;
            }

            const label_t idx3_bnd = device::global3(x, 0, z);
            const label_t idx3_yp1 = device::global3(x, 1, z);

            const scalar_t p = static_cast<scalar_t>(0);
            const scalar_t phi = static_cast<scalar_t>(0);
            const scalar_t ux = static_cast<scalar_t>(0);
            const scalar_t uy = physics::u_water;
            const scalar_t uz = static_cast<scalar_t>(0);

            d.p[idx3_bnd] = p;
            d.phi[idx3_bnd] = phi;
            d.ux[idx3_bnd] = ux;
            d.uy[idx3_bnd] = uy;
            d.uz[idx3_bnd] = uz;

            const scalar_t uu = static_cast<scalar_t>(1.5) * (ux * ux + uy * uy + uz * uz);

            device::constexpr_for<0, VelocitySet::Q()>(
                [&](const auto Q)
                {
                    if constexpr (VelocitySet::cy<Q>() == 1)
                    {
                        const label_t xx = x + static_cast<label_t>(VelocitySet::cx<Q>());
                        const label_t zz = z + static_cast<label_t>(VelocitySet::cz<Q>());

                        const label_t fluidNode = device::global3(xx, 1, zz);

                        constexpr scalar_t w = VelocitySet::w<Q>();
                        constexpr scalar_t cx = static_cast<scalar_t>(VelocitySet::cx<Q>());
                        constexpr scalar_t cy = static_cast<scalar_t>(VelocitySet::cy<Q>());
                        constexpr scalar_t cz = static_cast<scalar_t>(VelocitySet::cz<Q>());

                        const scalar_t cu = VelocitySet::as2() * (cx * ux + cy * uy + cz * uz);

                        const scalar_t feq = VelocitySet::f_eq<Q>(p, uu, cu);
                        const scalar_t fneq = VelocitySet::f_neq<Q>(d.pxx[fluidNode], d.pyy[fluidNode], d.pzz[fluidNode],
                                                                    d.pxy[fluidNode], d.pxz[fluidNode], d.pyz[fluidNode],
                                                                    d.ux[fluidNode], d.uy[fluidNode], d.uz[fluidNode]);

                        d.f[Q * size::cells() + fluidNode] = to_pop(feq + relaxation::omco_water() * fneq);
                    }
                });

            d.g[3 * size::cells() + idx3_yp1] = Phase::VelocitySet::w<3>() * phi * (static_cast<scalar_t>(1) + Phase::VelocitySet::as2() * uy);
        }

        __device__ static inline void applyOilInflow(LBMFields d) noexcept
        {
            const label_t x = threadIdx.x + blockIdx.x * blockDim.x;
            const label_t y = threadIdx.y + blockIdx.y * blockDim.y;

            if (x >= mesh::nx || y >= mesh::ny)
            {
                return;
            }

            const scalar_t dx = static_cast<scalar_t>(x) - geometry::center_x();
            const scalar_t dy = static_cast<scalar_t>(y) - geometry::y_pos();
            const scalar_t r2 = dx * dx + dy * dy;

            if (r2 > geometry::R2_oil())
            {
                return;
            }

            const label_t idx3_bnd = device::global3(x, y, 0);
            const label_t idx3_zp1 = device::global3(x, y, 1);

            const scalar_t p = static_cast<scalar_t>(0);
            const scalar_t phi = static_cast<scalar_t>(1);
            const scalar_t ux = static_cast<scalar_t>(0);
            const scalar_t uy = static_cast<scalar_t>(0);
            const scalar_t uz = physics::u_oil;

            d.p[idx3_bnd] = p;
            d.phi[idx3_bnd] = phi;
            d.ux[idx3_bnd] = ux;
            d.uy[idx3_bnd] = uy;
            d.uz[idx3_bnd] = uz;

            const scalar_t uu = static_cast<scalar_t>(1.5) * (ux * ux + uy * uy + uz * uz);

            device::constexpr_for<0, VelocitySet::Q()>(
                [&](const auto Q)
                {
                    if constexpr (VelocitySet::cz<Q>() == 1)
                    {
                        const label_t xx = x + static_cast<label_t>(VelocitySet::cx<Q>());
                        const label_t yy = y + static_cast<label_t>(VelocitySet::cy<Q>());

                        const label_t fluidNode = device::global3(xx, yy, 1);

                        constexpr scalar_t w = VelocitySet::w<Q>();
                        constexpr scalar_t cx = static_cast<scalar_t>(VelocitySet::cx<Q>());
                        constexpr scalar_t cy = static_cast<scalar_t>(VelocitySet::cy<Q>());
                        constexpr scalar_t cz = static_cast<scalar_t>(VelocitySet::cz<Q>());

                        const scalar_t cu = VelocitySet::as2() * (cx * ux + cy * uy + cz * uz);

                        const scalar_t feq = VelocitySet::f_eq<Q>(p, uu, cu);
                        const scalar_t fneq = VelocitySet::f_neq<Q>(d.pxx[fluidNode], d.pyy[fluidNode], d.pzz[fluidNode],
                                                                    d.pxy[fluidNode], d.pxz[fluidNode], d.pyz[fluidNode],
                                                                    d.ux[fluidNode], d.uy[fluidNode], d.uz[fluidNode]);

                        d.f[Q * size::cells() + fluidNode] = to_pop(feq + relaxation::omco_oil() * fneq);
                    }
                });

            d.g[5 * size::cells() + idx3_zp1] = Phase::VelocitySet::w<5>() * phi * (static_cast<scalar_t>(1) + Phase::VelocitySet::as2() * uz);
        }

        __device__ static inline void applyOutflowY(LBMFields d) noexcept
        {
            const label_t x = threadIdx.x + blockIdx.x * blockDim.x;
            const label_t z = threadIdx.y + blockIdx.y * blockDim.y;

            if (x == 0 || x == mesh::nx - 1 || z == 0 || z == mesh::nz - 1)
            {
                return;
            }

            const label_t idx3_bnd = device::global3(x, mesh::ny - 1, z);
            const label_t idx3_ym1 = device::global3(x, mesh::ny - 2, z);

            d.p[idx3_bnd] = d.p[idx3_ym1];
            d.phi[idx3_bnd] = d.phi[idx3_ym1];
            d.ux[idx3_bnd] = d.ux[idx3_ym1];
            d.uy[idx3_bnd] = d.uy[idx3_ym1];
            d.uz[idx3_bnd] = d.uz[idx3_ym1];

            const scalar_t p = d.p[idx3_bnd];
            const scalar_t phi = d.phi[idx3_bnd];
            const scalar_t ux = d.ux[idx3_bnd];
            const scalar_t uy = d.uy[idx3_bnd];
            const scalar_t uz = d.uz[idx3_bnd];

            const scalar_t uu = static_cast<scalar_t>(1.5) * (ux * ux + uy * uy + uz * uz);

            device::constexpr_for<0, VelocitySet::Q()>(
                [&](const auto Q)
                {
                    if constexpr (VelocitySet::cy<Q>() == -1)
                    {
                        const label_t xx = x + static_cast<label_t>(VelocitySet::cx<Q>());
                        const label_t zz = z + static_cast<label_t>(VelocitySet::cz<Q>());

                        const label_t fluidNode = device::global3(xx, mesh::ny - 2, zz);

                        constexpr scalar_t w = VelocitySet::w<Q>();
                        constexpr scalar_t cx = static_cast<scalar_t>(VelocitySet::cx<Q>());
                        constexpr scalar_t cy = static_cast<scalar_t>(VelocitySet::cy<Q>());
                        constexpr scalar_t cz = static_cast<scalar_t>(VelocitySet::cz<Q>());

                        const scalar_t cu = VelocitySet::as2() * (cx * ux + cy * uy + cz * uz);

                        const scalar_t feq = VelocitySet::f_eq<Q>(p, uu, cu);
                        const scalar_t fneq = VelocitySet::f_neq<Q>(d.pxx[fluidNode], d.pyy[fluidNode], d.pzz[fluidNode],
                                                                    d.pxy[fluidNode], d.pxz[fluidNode], d.pyz[fluidNode],
                                                                    d.ux[fluidNode], d.uy[fluidNode], d.uz[fluidNode]);

                        d.f[Q * size::cells() + fluidNode] = to_pop(feq + relaxation::omco_ref() * fneq);
                    }
                });

            d.g[4 * size::cells() + idx3_ym1] = Phase::VelocitySet::w<4>() * phi * (static_cast<scalar_t>(1) - Phase::VelocitySet::as2() * physics::u_oil);
        }

        __device__ static inline void applyOutflowZ(LBMFields d) noexcept
        {
            const label_t x = threadIdx.x + blockIdx.x * blockDim.x;
            const label_t y = threadIdx.y + blockIdx.y * blockDim.y;

            if (x >= mesh::nx || y >= mesh::ny)
            {
                return;
            }

            const label_t idx3_bnd = device::global3(x, y, mesh::nz - 1);
            const label_t idx3_zm1 = device::global3(x, y, mesh::nz - 2);

            d.p[idx3_bnd] = d.p[idx3_zm1];
            d.phi[idx3_bnd] = d.phi[idx3_zm1];
            d.ux[idx3_bnd] = d.ux[idx3_zm1];
            d.uy[idx3_bnd] = d.uy[idx3_zm1];
            d.uz[idx3_bnd] = d.uz[idx3_zm1];

            const scalar_t p = d.p[idx3_bnd];
            const scalar_t phi = d.phi[idx3_bnd];
            const scalar_t ux = d.ux[idx3_bnd];
            const scalar_t uy = d.uy[idx3_bnd];
            const scalar_t uz = d.uz[idx3_bnd];

            const scalar_t uu = static_cast<scalar_t>(1.5) * (ux * ux + uy * uy + uz * uz);

            device::constexpr_for<0, VelocitySet::Q()>(
                [&](const auto Q)
                {
                    if constexpr (VelocitySet::cz<Q>() == -1)
                    {
                        const label_t xx = x + static_cast<label_t>(VelocitySet::cx<Q>());
                        const label_t yy = y + static_cast<label_t>(VelocitySet::cy<Q>());

                        const label_t fluidNode = device::global3(xx, yy, mesh::nz - 2);

                        constexpr scalar_t w = VelocitySet::w<Q>();
                        constexpr scalar_t cx = static_cast<scalar_t>(VelocitySet::cx<Q>());
                        constexpr scalar_t cy = static_cast<scalar_t>(VelocitySet::cy<Q>());
                        constexpr scalar_t cz = static_cast<scalar_t>(VelocitySet::cz<Q>());

                        const scalar_t cu = VelocitySet::as2() * (cx * ux + cy * uy + cz * uz);

                        const scalar_t feq = VelocitySet::f_eq<Q>(p, uu, cu);
                        const scalar_t fneq = VelocitySet::f_neq<Q>(d.pxx[fluidNode], d.pyy[fluidNode], d.pzz[fluidNode],
                                                                    d.pxy[fluidNode], d.pxz[fluidNode], d.pyz[fluidNode],
                                                                    d.ux[fluidNode], d.uy[fluidNode], d.uz[fluidNode]);

                        d.f[Q * size::cells() + fluidNode] = to_pop(feq + relaxation::omco_ref() * fneq);
                    }
                });

            d.g[6 * size::cells() + idx3_zm1] = Phase::VelocitySet::w<6>() * phi * (static_cast<scalar_t>(1) - Phase::VelocitySet::as2() * physics::u_oil);
        }

    private:
        // No private methods
    };
}

#endif
