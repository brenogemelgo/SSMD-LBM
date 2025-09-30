#pragma once

__global__ 
void setFields(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = threadIdx.z + blockIdx.z * blockDim.z;

    if (x >= NX || y >= NY || z >= NZ) return;

    const idx_t idx3 = global3(x,y,z);

    d.ux[idx3] = 0.0f;
    d.uy[idx3] = 0.0f;
    d.uz[idx3] = 0.0f;
    d.phi[idx3] = 0.0f;
    d.rho[idx3] = 1.0f;
    d.normx[idx3] = 0.0f;
    d.normy[idx3] = 0.0f;
    d.normz[idx3] = 0.0f;
    d.pxx[idx3] = 0.0f;
    d.pyy[idx3] = 0.0f;
    d.pzz[idx3] = 0.0f;
    d.pxy[idx3] = 0.0f;
    d.pxz[idx3] = 0.0f;
    d.pyz[idx3] = 0.0f;
}

__global__ 
void setDistros(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = threadIdx.z + blockIdx.z * blockDim.z;

    if (x >= NX || y >= NY || z >= NZ) return;

    const idx_t idx3 = global3(x,y,z);

    const float uu = d.ux[idx3]*d.ux[idx3] + d.uy[idx3]*d.uy[idx3] + d.uz[idx3]*d.uz[idx3];

    #pragma unroll FLINKS
    for (idx_t Q = 0; Q < FLINKS; ++Q) {
        d.f[global4(x,y,z,Q)] = computeFeq(d.rho[idx3],d.ux[idx3],d.uy[idx3],d.uz[idx3],uu,Q);
    }
    #pragma unroll GLINKS
    for (idx_t Q = 0; Q < GLINKS; ++Q) {
        d.g[global4(x,y,z,Q)] = computeGeq(d.phi[idx3],d.ux[idx3],d.uy[idx3],d.uz[idx3],Q);
    }
} 

__global__ 
void setOilJet(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = 0;

    if (x >= NX || y >= NY) return;

    const float dx = static_cast<float>(x) - CENTER_X;
    const float dy = static_cast<float>(y) - Y_POS;
    const float r2 = dx*dx + dy*dy;
    const float R = 0.5f * static_cast<float>(DIAM_OIL);
    if (r2 > R*R) return;

    const idx_t idx3_in = global3(x,y,z);
    d.uz[idx3_in] = U_OIL;
    d.phi[idx3_in] = 1.0f;
}

__global__ 
void setWaterJet(
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
    d.uy[idx3_in] = U_WATER;
    d.phi[idx3_in] = 0.0f;
}


