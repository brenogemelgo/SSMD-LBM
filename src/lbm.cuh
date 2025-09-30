#pragma once

__global__ 
void computePhase(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = threadIdx.z + blockIdx.z * blockDim.z;

    if (x >= NX || y >= NY || z >= NZ || 
        x == 0 || x == NX-1 || 
        y == 0 || y == NY-1 || 
        z == 0 || z == NZ-1) return;

    const idx_t idx3 = global3(x,y,z);

    const float phi = d.g[idx3] + 
                      d.g[PLANE+idx3] + 
                      d.g[PLANE2+idx3] + 
                      d.g[PLANE3+idx3] + 
                      d.g[PLANE4+idx3] + 
                      d.g[PLANE5+idx3] + 
                      d.g[PLANE6+idx3];
                      
    d.phi[idx3] = phi;
}

__global__ 
void forceStreamCollide(
    LBMFields d
) {
    const idx_t x = threadIdx.x + blockIdx.x * blockDim.x;
    const idx_t y = threadIdx.y + blockIdx.y * blockDim.y;
    const idx_t z = threadIdx.z + blockIdx.z * blockDim.z;

    if (x >= NX || y >= NY || z >= NZ || 
        x == 0 || x == NX-1 || 
        y == 0 || y == NY-1 || 
        z == 0 || z == NZ-1) return;

    const idx_t idx3 = global3(x,y,z);

    const idx_t xp1 = idx3 + 1;
    const idx_t xm1 = idx3 - 1;
    const idx_t yp1 = idx3 + NX;
    const idx_t ym1 = idx3 - NX;
    const idx_t zp1 = idx3 + STRIDE;
    const idx_t zm1 = idx3 - STRIDE;

    const idx_t xp1yp1 = xp1 + NX;
    const idx_t xm1ym1 = xm1 - NX;
    const idx_t xp1ym1 = xp1 - NX;
    const idx_t xm1yp1 = xm1 + NX;

    const idx_t xp1zp1 = xp1 + STRIDE;
    const idx_t xm1zm1 = xm1 - STRIDE;
    const idx_t xp1zm1 = xp1 - STRIDE;
    const idx_t xm1zp1 = xm1 + STRIDE;

    const idx_t yp1zp1 = yp1 + STRIDE;
    const idx_t ym1zm1 = ym1 - STRIDE;
    const idx_t yp1zm1 = yp1 - STRIDE;
    const idx_t ym1zp1 = ym1 + STRIDE;


    float sumGradX = W_1 * (d.phi[xp1] - d.phi[xm1]) +
                     W_2 * ((d.phi[xp1yp1] - d.phi[xm1ym1]) +
                            (d.phi[xp1zp1] - d.phi[xm1zm1]) +
                            (d.phi[xp1ym1] - d.phi[xm1yp1]) +
                            (d.phi[xp1zm1] - d.phi[xm1zp1]));

    float sumGradY = W_1 * (d.phi[yp1] - d.phi[ym1]) +
                     W_2 * ((d.phi[xp1yp1] - d.phi[xm1ym1]) +
                            (d.phi[yp1zp1] - d.phi[ym1zm1]) +
                            (d.phi[xm1yp1] - d.phi[xp1ym1]) +
                            (d.phi[yp1zm1] - d.phi[ym1zp1]));

    float sumGradZ = W_1 * (d.phi[zp1] - d.phi[zm1]) +
                     W_2 * ((d.phi[xp1zp1] - d.phi[xm1zm1]) +
                            (d.phi[yp1zp1] - d.phi[ym1zm1]) +
                            (d.phi[xm1zp1] - d.phi[xp1zm1]) +
                            (d.phi[ym1zp1] - d.phi[yp1zm1]));

    #if defined(D3Q27)
    const idx_t xp1yp1zp1 = xp1yp1 + STRIDE;
    const idx_t xm1ym1zm1 = xm1ym1 - STRIDE;
    const idx_t xp1yp1zm1 = xp1yp1 - STRIDE;
    const idx_t xm1ym1zp1 = xm1ym1 + STRIDE;

    const idx_t xp1ym1zp1 = xp1ym1 + STRIDE;
    const idx_t xm1yp1zm1 = xm1yp1 - STRIDE;
    const idx_t xp1ym1zm1 = xp1ym1 - STRIDE;
    const idx_t xm1yp1zp1 = xm1yp1 + STRIDE;
    
    sumGradX += W_3 * ((d.phi[xp1yp1zp1] - d.phi[xm1ym1zm1]) +
                       (d.phi[xp1yp1zm1] - d.phi[xm1ym1zp1]) +
                       (d.phi[xp1ym1zp1] - d.phi[xm1yp1zm1]) +
                       (d.phi[xp1ym1zm1] - d.phi[xm1yp1zp1]));

    sumGradY += W_3 * ((d.phi[xp1yp1zp1] - d.phi[xm1ym1zm1]) +
                       (d.phi[xp1yp1zm1] - d.phi[xm1ym1zp1]) +
                       (d.phi[xm1yp1zm1] - d.phi[xp1ym1zp1]) +
                       (d.phi[xm1yp1zp1] - d.phi[xp1ym1zm1]));

    sumGradZ += W_3 * ((d.phi[xp1yp1zp1] - d.phi[xm1ym1zm1]) +
                       (d.phi[xm1ym1zp1] - d.phi[xp1yp1zm1]) +
                       (d.phi[xp1ym1zp1] - d.phi[xm1yp1zm1]) +
                       (d.phi[xm1yp1zp1] - d.phi[xp1ym1zm1]));
    #endif 
        
    const float gradX = 3.0f * sumGradX;
    const float gradY = 3.0f * sumGradY;
    const float gradZ = 3.0f * sumGradZ;
    
    const float ind = sqrtf(gradX*gradX + gradY*gradY + gradZ*gradZ);
    const float invInd = 1.0f / (ind + 1e-9f);

    const float normX = gradX * invInd;
    const float normY = gradY * invInd;
    const float normZ = gradZ * invInd;

    d.normx[idx3] = normX;
    d.normy[idx3] = normY;
    d.normz[idx3] = normZ;

    float sumCurvX = W_1 * (d.normx[xp1] - d.normx[xm1]) +
                     W_2 * (d.normx[xp1yp1] - d.normx[xm1ym1] +
                            d.normx[xp1zp1] - d.normx[xm1zm1] +
                            d.normx[xp1ym1] - d.normx[xm1yp1] +
                            d.normx[xp1zm1] - d.normx[xm1zp1]);

    float sumCurvY = W_1 * (d.normy[yp1] - d.normy[ym1]) +
                     W_2 * (d.normy[xp1yp1] - d.normy[xm1ym1] +
                            d.normy[yp1zp1] - d.normy[ym1zm1] +
                            d.normy[xm1yp1] - d.normy[xp1ym1] +
                            d.normy[yp1zm1] - d.normy[ym1zp1]);

    float sumCurvZ = W_1 * (d.normz[zp1] - d.normz[zm1]) +
                     W_2 * (d.normz[xp1zp1] - d.normz[xm1zm1] +
                            d.normz[yp1zp1] - d.normz[ym1zm1] +
                            d.normz[xm1zp1] - d.normz[xp1zm1] +
                            d.normz[ym1zp1] - d.normz[yp1zm1]);
    #if defined(D3Q27)
    sumCurvX += W_3 * (d.normx[xp1yp1zp1] - d.normx[xm1ym1zm1] +
                       d.normx[xp1yp1zm1] - d.normx[xm1ym1zp1] +
                       d.normx[xp1ym1zp1] - d.normx[xm1yp1zm1] +
                       d.normx[xp1ym1zm1] - d.normx[xm1yp1zp1]);

    sumCurvY += W_3 * (d.normy[xp1yp1zp1] - d.normy[xm1ym1zm1] +
                       d.normy[xp1yp1zm1] - d.normy[xm1ym1zp1] +
                       d.normy[xm1yp1zm1] - d.normy[xp1ym1zp1] +
                       d.normy[xm1yp1zp1] - d.normy[xp1ym1zm1]);

    sumCurvZ += W_3 * (d.normz[xp1yp1zp1] - d.normz[xm1ym1zm1] +
                       d.normz[xm1ym1zp1] - d.normz[xp1yp1zm1] +
                       d.normz[xp1ym1zp1] - d.normz[xm1yp1zm1] +
                       d.normz[xm1yp1zp1] - d.normz[xp1ym1zm1]);
    #endif 
    float curvature = -3.0f * (sumCurvX + sumCurvY + sumCurvZ);   

    const float stCurv = SIGMA * curvature * ind;
    const float ffx = stCurv * normX;
    const float ffy = stCurv * normY;
    const float ffz = stCurv * normZ;
        
    const float pop0 = from_pop(d.f[idx3]);         
    const float pop1 = from_pop(d.f[PLANE+idx3]);   
    const float pop2 = from_pop(d.f[PLANE2+idx3]);  
    const float pop3 = from_pop(d.f[PLANE3+idx3]);  
    const float pop4 = from_pop(d.f[PLANE4+idx3]);  
    const float pop5 = from_pop(d.f[PLANE5+idx3]);  
    const float pop6 = from_pop(d.f[PLANE6+idx3]);  
    const float pop7 = from_pop(d.f[PLANE7+idx3]);  
    const float pop8 = from_pop(d.f[PLANE8+idx3]);  
    const float pop9 = from_pop(d.f[PLANE9+idx3]);  
    const float pop10 = from_pop(d.f[PLANE10+idx3]);
    const float pop11 = from_pop(d.f[PLANE11+idx3]);
    const float pop12 = from_pop(d.f[PLANE12+idx3]); 
    const float pop13 = from_pop(d.f[PLANE13+idx3]); 
    const float pop14 = from_pop(d.f[PLANE14+idx3]);
    const float pop15 = from_pop(d.f[PLANE15+idx3]); 
    const float pop16 = from_pop(d.f[PLANE16+idx3]); 
    const float pop17 = from_pop(d.f[PLANE17+idx3]);
    const float pop18 = from_pop(d.f[PLANE18+idx3]); 
    #if defined(D3Q27)
    const float pop19 = from_pop(d.f[PLANE19+idx3]); 
    const float pop20 = from_pop(d.f[PLANE20+idx3]); 
    const float pop21 = from_pop(d.f[PLANE21+idx3]); 
    const float pop22 = from_pop(d.f[PLANE22+idx3]); 
    const float pop23 = from_pop(d.f[PLANE23+idx3]); 
    const float pop24 = from_pop(d.f[PLANE24+idx3]);
    const float pop25 = from_pop(d.f[PLANE25+idx3]);
    const float pop26 = from_pop(d.f[PLANE26+idx3]); 
    #endif 

    float rho = pop0 + pop1 + pop2 + pop3 + pop4 + pop5 + pop6 + pop7 + pop8 + pop9 + pop10 + pop11 + pop12 + pop13 + pop14 + pop15 + pop16 + pop17 + pop18;
    #if defined(D3Q27)
    rho += pop19 + pop20 + pop21 + pop22 + pop23 + pop24 + pop25 + pop26;
    #endif
    rho += 1.0f; 
    d.rho[idx3] = rho;

    const float invRho = 1.0f / rho;
    
    #if defined(D3Q19)
    float ux = invRho * (pop1 - pop2 + pop7 - pop8 + pop9 - pop10 + pop13 - pop14 + pop15 - pop16);
    float uy = invRho * (pop3 - pop4 + pop7 - pop8 + pop11 - pop12 + pop14 - pop13 + pop17 - pop18);
    float uz = invRho * (pop5 - pop6 + pop9 - pop10 + pop11 - pop12 + pop16 - pop15 + pop18 - pop17);
    #elif defined(D3Q27)
    float ux = invRho * (pop1 - pop2 + pop7 - pop8 + pop9 - pop10 + pop13 - pop14 + pop15 - pop16 + pop19 - pop20 + pop21 - pop22 + pop23 - pop24 + pop26 - pop25);
    float uy = invRho * (pop3 - pop4 + pop7 - pop8  + pop11 - pop12 + pop14 - pop13 + pop17 - pop18 + pop19 - pop20 + pop21 - pop22 + pop24 - pop23 + pop25 - pop26);
    float uz = invRho * (pop5 - pop6 + pop9 - pop10 + pop11 - pop12 + pop16 - pop15 + pop18 - pop17 + pop19 - pop20 + pop22 - pop21 + pop23 - pop24 + pop25 - pop26);
    #endif
    
    ux += ffx * 0.5f * invRho;
    uy += ffy * 0.5f * invRho;
    uz += ffz * 0.5f * invRho;

    d.ux[idx3] = ux; 
    d.uy[idx3] = uy; 
    d.uz[idx3] = uz;

    #include "../include/momentumFluxAU.cuh" // Agressive Unrolling
    //#include "../include/momentumFluxNT.cuh" // No Temps

    d.pxx[idx3] = pxx;
    d.pyy[idx3] = pyy;
    d.pzz[idx3] = pzz;
    d.pxy[idx3] = pxy;
    d.pxz[idx3] = pxz;   
    d.pyz[idx3] = pyz;

    float omcoLocal;
    #if defined(VISC_CONTRAST)
    {
        const float phi = d.phi[idx3]; 
        const float nuLocal = fmaf(phi, (VISC_OIL - VISC_WATER), VISC_WATER);
        const float omegaPhys = 1.0f / (0.5f + 3.0f * nuLocal);

        omcoLocal = 1.0f - fminf(omegaPhys, cubicSponge(z));
    }
    #else
    {
        omcoLocal = 1.0f - cubicSponge(z);
    }
    #endif

    #include "../include/streamCollide.cuh" 

    { // ====================================== ADVECTION-DIFFUSION ====================================== //
        #if !defined(VISC_CONTRAST)
        const float phi = d.phi[idx3];
        #endif

        // Q0
        d.g[idx3] = WG_0 * phi;

        const float multPhi = WG_1 * phi;
        const float phiNorm = GAMMA * multPhi * (1.0f - phi);
        const float a4      = 4.0f * multPhi;

        // -------------------------------- X+ (Q1)
        float geq = multPhi + a4 * ux;
        float hi = phiNorm * d.normx[idx3];
        d.g[global4(x+1,y,z,1)] = geq + hi;
        
        // -------------------------------- X- (Q2)
        geq = multPhi - a4 * ux;
        d.g[global4(x-1,y,z,2)] = geq - hi;

        // -------------------------------- Y+ (Q3)
        geq = multPhi + a4 * uy;
        hi = phiNorm * d.normy[idx3];
        d.g[global4(x,y+1,z,3)] = geq + hi;

        // -------------------------------- Y- (Q4)
        geq = multPhi - a4 * uy;
        d.g[global4(x,y-1,z,4)] = geq - hi;

        // -------------------------------- Z+ (Q5)
        geq = multPhi + a4 * uz;
        hi = phiNorm * d.normz[idx3];
        d.g[global4(x,y,z+1,5)] = geq + hi;

        // -------------------------------- Z- (Q6)
        geq = multPhi - a4 * uz;
        d.g[global4(x,y,z-1,6)] = geq - hi;
    } // ============================================= END ============================================= //
}



