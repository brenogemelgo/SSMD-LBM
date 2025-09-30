#pragma once
#include "../include/cudaUtils.cuh"
#include "../include/velocitySets.cuh"

//#define RUN_MODE
#define SAMPLE_MODE
//#define DEBUG_MODE

#if defined(RUN_MODE)

    constexpr int MACRO_SAVE = 100;
    constexpr int NSTEPS = 50000;
    
#elif defined(SAMPLE_MODE)

    constexpr int MACRO_SAVE = 100;
    constexpr int NSTEPS = 1000;

#elif defined(DEBUG_MODE)

    constexpr int MACRO_SAVE = 1;
    constexpr int NSTEPS = 0;

#endif

constexpr idx_t MESH = 128;
constexpr idx_t NX   = MESH;
constexpr idx_t NY   = MESH*2;
constexpr idx_t NZ   = MESH*2;

constexpr float CENTER_X = (NX-1) * 0.5f;
constexpr float CENTER_Y = (NY-1) * 0.5f;
constexpr float CENTER_Z = (NZ-1) * 0.5f;

constexpr float U_OIL      = 0.05f;
constexpr int   DIAM_OIL   = 20;

constexpr float U_WATER    = 0.05f; 
constexpr int   DIAM_WATER = 20;

constexpr int REYNOLDS_WATER = 5000;
constexpr int REYNOLDS_OIL   = 5000;     
constexpr int WEBER = 500;  

constexpr float VISC_OIL   = (U_OIL * DIAM_OIL) / REYNOLDS_OIL;     
constexpr float VISC_WATER = (U_WATER * DIAM_WATER) / REYNOLDS_WATER; 

constexpr float OMEGA_OIL   = 1.0f / (0.5f + 3.0f * VISC_OIL);
constexpr float OMEGA_WATER = 1.0f / (0.5f + 3.0f * VISC_WATER);

constexpr float OMCO_YMIN = 1.0f - OMEGA_WATER;
constexpr float OMCO_ZMIN = 1.0f - OMEGA_OIL;

constexpr float SIGMA = (U_OIL * U_OIL * DIAM_OIL) / WEBER; 
constexpr float GAMMA = 0.3f * 3.0f; 
constexpr float CSSQ  = 1.0f / 3.0f;  

constexpr float K          = 50.0f;
constexpr float P          = 3.0f;            
constexpr int SPONGE_CELLS = static_cast<int>(NZ/12);      
static_assert(SPONGE_CELLS > 0, "SPONGE_CELLS must be > 0");

constexpr float SPONGE     = static_cast<float>(SPONGE_CELLS) / static_cast<float>(NZ-1);
constexpr float Z_START    = static_cast<float>(NZ-1-SPONGE_CELLS) / static_cast<float>(NZ-1);
constexpr float INV_NZ_M1  = 1.0f / static_cast<float>(NZ-1);
constexpr float INV_SPONGE = 1.0f / SPONGE;

constexpr float OMEGA_REF = OMEGA_WATER;
constexpr float VISC_REF = VISC_WATER; 

constexpr float OMEGA_MAX  = 1.0f / (0.5f + 3.0f * VISC_REF * (K + 1.0f));
constexpr float OMCO_MAX   = 1.0f - OMEGA_MAX; 
constexpr float OMEGA_DELTA = OMEGA_MAX - OMEGA_REF; 

constexpr idx_t STRIDE = NX * NY;
constexpr idx_t PLANE  = NX * NY * NZ;

constexpr idx_t PLANE2  = 2 * PLANE;
constexpr idx_t PLANE3  = 3 * PLANE;
constexpr idx_t PLANE4  = 4 * PLANE;
constexpr idx_t PLANE5  = 5 * PLANE;
constexpr idx_t PLANE6  = 6 * PLANE;
constexpr idx_t PLANE7  = 7 * PLANE;
constexpr idx_t PLANE8  = 8 * PLANE;
constexpr idx_t PLANE9  = 9 * PLANE;
constexpr idx_t PLANE10 = 10 * PLANE;
constexpr idx_t PLANE11 = 11 * PLANE;
constexpr idx_t PLANE12 = 12 * PLANE;   
constexpr idx_t PLANE13 = 13 * PLANE;
constexpr idx_t PLANE14 = 14 * PLANE;
constexpr idx_t PLANE15 = 15 * PLANE;
constexpr idx_t PLANE16 = 16 * PLANE;
constexpr idx_t PLANE17 = 17 * PLANE;
constexpr idx_t PLANE18 = 18 * PLANE;
#if defined(D3Q27)
constexpr idx_t PLANE19 = 19 * PLANE;
constexpr idx_t PLANE20 = 20 * PLANE;
constexpr idx_t PLANE21 = 21 * PLANE;
constexpr idx_t PLANE22 = 22 * PLANE;
constexpr idx_t PLANE23 = 23 * PLANE;
constexpr idx_t PLANE24 = 24 * PLANE;
constexpr idx_t PLANE25 = 25 * PLANE;
constexpr idx_t PLANE26 = 26 * PLANE;
#endif                     

 