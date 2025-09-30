#pragma once
#include "cudaUtils.cuh"

#if defined(D3Q19) //             0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 

    __constant__ ci_t CIX[19] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0 };
    __constant__ ci_t CIY[19] = { 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1 };
    __constant__ ci_t CIZ[19] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1 };
    __constant__ float W[19] = { 1.0f / 3.0f, 
                                 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f,
                                 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f };
    static constexpr float W_0 = 1.0f / 3.0f;  // 0
    static constexpr float W0_45 = W_0 * 4.5f;
    static constexpr float W_1 = 1.0f / 18.0f; // 1 to 6
    static constexpr float W1_45 = W_1 * 4.5f;
    static constexpr float W_2 = 1.0f / 36.0f; // 7 to 18
    static constexpr float W2_45 = W_2 * 4.5f;
    static constexpr idx_t FLINKS = 19;

#elif defined(D3Q27) //           0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26

    __constant__ ci_t CIX[27] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1 };
    __constant__ ci_t CIY[27] = { 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1 };
    __constant__ ci_t CIZ[27] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1 };
    __constant__ float W[27] = { 8.0f / 27.0f,
                                 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 
                                 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 
                                 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f };
    static constexpr float W_0 = 8.0f / 27.0f;  // 0
    static constexpr float W0_45 = W_0 * 4.5f;
    static constexpr float W_1 = 2.0f / 27.0f;  // 1 to 6
    static constexpr float W1_45 = W_1 * 4.5f;
    static constexpr float W_2 = 1.0f / 54.0f;  // 7 to 18
    static constexpr float W2_45 = W_2 * 4.5f;
    static constexpr float W_3 = 1.0f / 216.0f; // 19 to 26
    static constexpr float W3_45 = W_3 * 4.5f;
    static constexpr idx_t FLINKS = 27;
    
#endif

__constant__ float W_G[7] = { 1.0f / 4.0f, 
                              1.0f / 8.0f, 1.0f / 8.0f, 
                              1.0f / 8.0f, 1.0f / 8.0f, 
                              1.0f / 8.0f, 1.0f / 8.0f };
constexpr float WG_0 = 1.0f / 4.0f; // 0
constexpr float WG_1 = 1.0f / 8.0f; // 1 to 6
constexpr idx_t GLINKS = 7;   
