#pragma once

__global__ 
void applyOilInflow(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = 0;

    if (x >= NX || y >= NY) return;

    const float dx = static_cast<float>(x) - CENTER_X;
    const float dy = static_cast<float>(y) - Y_POS;
    const float radialDist = sqrtf(dx*dx + dy*dy);
    const float radius = 0.5f * static_cast<float>(DIAM_OIL);
    if (radialDist > radius) return;

    const idx_t idx3_in = global3(x,y,z);
    
    emitInflowZ<5,
        5,9,11,16,18
    #if defined(D3Q27)
        ,19,22,23,25
    #endif
    >(d, x,y,z, idx3_in);
}

__global__
void applyWaterInflow(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = 0;
    const idx_t z = threadIdx.y + blockIdx.y * blockDim.y;

    if (x >= NX || z >= NZ) return;

    const float dx = static_cast<float>(x) - CENTER_X;
    const float dz = static_cast<float>(z) - Z_POS;
    const float r2 = dx*dx + dz*dz;
    const float R = 0.5f * static_cast<float>(DIAM_WATER);
    if (r2 > R*R) return;

    const idx_t idx3_in = global3(x,y,z);
    
    emitInflowY<3,
        3,7,11,14,17
    #if defined(D3Q27)
        ,19,21,24,25
    #endif
    >(d, x,y,z, idx3_in);
}

__global__ 
void applyOutflowZ(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = NZ-1;

    if (x >= NX || y >= NY) return;

    const idx_t idx3_zm1 = global3(x,y,z-1);
    d.phi[global3(x,y,z)] = d.phi[idx3_zm1];

    const float uxOut = d.ux[idx3_zm1];
    const float uyOut = d.uy[idx3_zm1];
    const float uzOut = d.uz[idx3_zm1];

    emitOutflowZ<6,
        6,10,12,15,17
    #if defined(D3Q27)
        ,20,21,24,26
    #endif
    >(d, x,y,z, uxOut, uyOut, uzOut);
}

__global__ 
void applyOutflowY(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = NY-1;
    const idx_t z = threadIdx.y + blockIdx.y * blockDim.y;

    if (x == 0 || x == NX-1 || z == 0 || z == NZ-1) return;

    const idx_t idx3_ym1 = global3(x,y-1,z);
    d.phi[global3(x,y,z)] = d.phi[idx3_ym1];

    const float uxOut = d.ux[idx3_ym1];
    const float uyOut = d.uy[idx3_ym1];
    const float uzOut = d.uz[idx3_ym1];

    emitOutflowY<4,
        4,8,12,13,18
    #if defined(D3Q27)
        ,20,22,23,26
    #endif
    >(d, x,y,z, uxOut, uyOut, uzOut);
}

__global__ 
void periodicX(
    LBMFields d
) {
    const idx_t y = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t z = threadIdx.y + blockIdx.y * blockDim.y;

    if (y <= 0 || y >= NY-1 || z <= 0 || z >= NZ-1) return;

    const idx_t bL = global3(1,y,z);
    const idx_t bR = global3(NX-2,y,z);

    // positive x contributions
    d.f[     PLANE + bL] = d.f[     PLANE + bR];
    d.f[7  * PLANE + bL] = d.f[7  * PLANE + bR];
    d.f[9  * PLANE + bL] = d.f[9  * PLANE + bR];
    d.f[13 * PLANE + bL] = d.f[13 * PLANE + bR];
    d.f[15 * PLANE + bL] = d.f[15 * PLANE + bR];
    #if defined(D3Q27)
    d.f[19 * PLANE + bL] = d.f[19 * PLANE + bR];
    d.f[21 * PLANE + bL] = d.f[21 * PLANE + bR];
    d.f[23 * PLANE + bL] = d.f[23 * PLANE + bR];
    d.f[26 * PLANE + bL] = d.f[26 * PLANE + bR];
    #endif
    d.g[     PLANE + bL] = d.g[     PLANE + bR];

    // negative x contributions
    d.f[2  * PLANE + bR] = d.f[2  * PLANE + bL];
    d.f[8  * PLANE + bR] = d.f[8  * PLANE + bL];
    d.f[10 * PLANE + bR] = d.f[10 * PLANE + bL];
    d.f[14 * PLANE + bR] = d.f[14 * PLANE + bL];
    d.f[16 * PLANE + bR] = d.f[16 * PLANE + bL];
    #if defined(D3Q27)
    d.f[20 * PLANE + bR] = d.f[20 * PLANE + bL];
    d.f[22 * PLANE + bR] = d.f[22 * PLANE + bL];
    d.f[24 * PLANE + bR] = d.f[24 * PLANE + bL];
    d.f[25 * PLANE + bR] = d.f[25 * PLANE + bL];
    #endif
    d.g[2  * PLANE + bR] = d.g[2  * PLANE + bL];

    // ghost cells
    d.phi[global3(0,y,z)]    = d.phi[bR];
    d.phi[global3(NX-1,y,z)] = d.phi[bL];
}

__global__ 
void periodicY(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t z = threadIdx.y + blockIdx.y * blockDim.y;

    if (x <= 0 || x >= NX-1 || z <= 0 || z >= NZ-1) return;

    const idx_t bB = global3(x,1,z);
    const idx_t bT = global3(x,NY-2,z);

    // positive y contributions
    d.f[3  * PLANE + bB] = d.f[3  * PLANE + bT];
    d.f[7  * PLANE + bB] = d.f[7  * PLANE + bT];
    d.f[11 * PLANE + bB] = d.f[11 * PLANE + bT];
    d.f[14 * PLANE + bB] = d.f[14 * PLANE + bT];
    d.f[17 * PLANE + bB] = d.f[17 * PLANE + bT];
    #if defined(D3Q27)
    d.f[19 * PLANE + bB] = d.f[19 * PLANE + bT];
    d.f[21 * PLANE + bB] = d.f[21 * PLANE + bT];
    d.f[24 * PLANE + bB] = d.f[24 * PLANE + bT];
    d.f[25 * PLANE + bB] = d.f[25 * PLANE + bT];
    #endif
    d.g[3  * PLANE + bB] = d.g[3  * PLANE + bT];

    // negative y contributions
    d.f[4  * PLANE + bT] = d.f[4  * PLANE + bB];
    d.f[8  * PLANE + bT] = d.f[8  * PLANE + bB];
    d.f[12 * PLANE + bT] = d.f[12 * PLANE + bB];
    d.f[13 * PLANE + bT] = d.f[13 * PLANE + bB];
    d.f[18 * PLANE + bT] = d.f[18 * PLANE + bB];
    #if defined(D3Q27)
    d.f[20 * PLANE + bT] = d.f[20 * PLANE + bB];
    d.f[22 * PLANE + bT] = d.f[22 * PLANE + bB];
    d.f[23 * PLANE + bT] = d.f[23 * PLANE + bB];
    d.f[26 * PLANE + bT] = d.f[26 * PLANE + bB];
    #endif
    d.g[4   * PLANE + bT] = d.g[4   * PLANE + bB];

    // ghost cells
    d.phi[global3(x,0,z)]    = d.phi[bT];
    d.phi[global3(x,NY-1,z)] = d.phi[bB];
}


