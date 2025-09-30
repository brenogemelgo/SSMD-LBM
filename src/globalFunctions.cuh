#pragma once
#include "constants.cuh"

__device__ __forceinline__ 
idx_t global3(
    const idx_t x, 
    const idx_t y, 
    const idx_t z
) {
    return x + y * NX + z * STRIDE;
}

__device__ __forceinline__ 
idx_t global4(
    const idx_t x, 
    const idx_t y, 
    const idx_t z, 
    const idx_t Q
) {
    return Q * PLANE + global3(x,y,z);
}

#include "../include/discontinued/constexprBased/lbmFunctions.cuh"
#include "../include/discontinued/functionBased/lbmFunctions.cuh"

__device__ __forceinline__ 
float smoothstep(
    float edge0, 
    float edge1, 
    float x
) {
    x = __saturatef((x - edge0) / (edge1 - edge0));
    return x * x * (3.0f - 2.0f * x);
}

template<typename T, idx_t... Qs>
__device__ __forceinline__
void copyDirs(
    T* __restrict__ arr, 
    const idx_t dst, 
    const idx_t src
) {
    ((arr[Qs * PLANE + dst] = arr[Qs * PLANE + src]), ...);
}

struct LBMFields {
    float *rho;
    float *phi;
    float *ux;
    float *uy;
    float *uz;
    float *pxx;
    float *pyy;
    float *pzz; 
    float *pxy; 
    float *pxz; 
    float *pyz;
    float *normx; 
    float *normy; 
    float *normz;
    pop_t *f; 
    float *g; 
};
LBMFields lbm;

#if defined(D_FIELDS)
struct DerivedFields {
    float *vorticity_mag;
    float *velocity_mag;
};
DerivedFields dfields;
#endif 