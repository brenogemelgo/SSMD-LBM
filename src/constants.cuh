#pragma once
#include "../helpers/cudaUtils.cuh"
#include "../helpers/velocitySets.cuh"

#define RUN_MODE
//#define SAMPLE_MODE
//#define DEBUG_MODE

#if defined(RUN_MODE)

constexpr int MACRO_SAVE = 1000;
constexpr int NSTEPS = 100000;
    
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

// ================= INFLOW PARAMETERS ================== //

constexpr float U_WATER    = 0.06f; 
constexpr int   DIAM_WATER = 13;

constexpr float U_OIL    = 0.05f;
constexpr int   DIAM_OIL = 13;

constexpr int REYNOLDS_WATER = 1400;
constexpr int REYNOLDS_OIL   = 450;     

constexpr int WEBER = 500;  

constexpr float Y_POS = 0.5f * CENTER_Y; // y position of oil inflow
constexpr float Z_POS = 0.8f * CENTER_Z; // z position of water inflow

// ===================================================== //

constexpr float UU_OIL    = U_OIL * U_OIL;
constexpr float UU_WATER  = U_WATER * U_WATER;

constexpr float VISC_OIL   = (U_OIL * DIAM_OIL) / REYNOLDS_OIL;     
constexpr float VISC_WATER = (U_WATER * DIAM_WATER) / REYNOLDS_WATER; 

constexpr float OMEGA_OIL   = 1.0f / (0.5f + 3.0f * VISC_OIL);
constexpr float OMEGA_WATER = 1.0f / (0.5f + 3.0f * VISC_WATER);

constexpr float OMCO_YMIN = 1.0f - OMEGA_WATER;
constexpr float OMCO_ZMIN = 1.0f - OMEGA_OIL;

constexpr float SIGMA = (U_OIL * U_OIL * DIAM_OIL) / WEBER; 
constexpr float GAMMA = 0.9f; 
constexpr float CSSQ  = 1.0f / 3.0f;  

constexpr float OMEGA_REF = OMEGA_WATER;
constexpr float VISC_REF = VISC_WATER; 

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

 