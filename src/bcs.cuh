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
    const float dy = static_cast<float>(y) - CENTER_Y;
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
    const float dz = static_cast<float>(z) - CENTER_Z;
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
    copyDirs<pop_t,1,7,9,13,15>(d.f,bL,bR);   
    #if defined(D3Q27)
    copyDirs<pop_t,19,21,23,26>(d.f,bL,bR);
    #endif 

    // negative x contributions
    copyDirs<pop_t,2,8,10,14,16>(d.f,bR,bL); 
    #if defined(D3Q27)
    copyDirs<pop_t,20,22,24,25>(d.f,bR,bL);
    #endif 

    d.g[PLANE+bL] = d.g[PLANE+bR];
    d.g[PLANE2+bR] = d.g[PLANE2+bL];
    d.phi[global3(0,y,z)] = d.phi[bR];
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

    d.g[PLANE3+bB] = d.g[PLANE3+bT];
    d.g[PLANE4+bT] = d.g[PLANE4+bB];
    d.phi[global3(x,0,z)] = d.phi[bT];
    d.phi[global3(x,NY-1,z)] = d.phi[bB];

    // positive y contributions
    copyDirs<pop_t,3,7,11,14,17>(d.f,bB,bT);
    #if defined(D3Q27)
    copyDirs<pop_t,19,21,24,25>(d.f,bB,bT);
    #endif 

    // negative y contributions
    copyDirs<pop_t,4,8,12,13,18>(d.f,bT,bB);
    #if defined(D3Q27)
    copyDirs<pop_t,20,22,23,26>(d.f,bT,bB);
    #endif 

}

