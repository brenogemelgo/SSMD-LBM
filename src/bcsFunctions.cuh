#pragma once

template<int DIR> struct dir_of;
template<> struct dir_of<0>  { static constexpr int dx= 0, dy= 0, dz= 0; };
template<> struct dir_of<1>  { static constexpr int dx= 1, dy= 0, dz= 0; };
template<> struct dir_of<2>  { static constexpr int dx=-1, dy= 0, dz= 0; };
template<> struct dir_of<3>  { static constexpr int dx= 0, dy= 1, dz= 0; };
template<> struct dir_of<4>  { static constexpr int dx= 0, dy=-1, dz= 0; };
template<> struct dir_of<5>  { static constexpr int dx= 0, dy= 0, dz= 1; };
template<> struct dir_of<6>  { static constexpr int dx= 0, dy= 0, dz=-1; };
template<> struct dir_of<7>  { static constexpr int dx= 1, dy= 1, dz= 0; };
template<> struct dir_of<8>  { static constexpr int dx=-1, dy=-1, dz= 0; };
template<> struct dir_of<9>  { static constexpr int dx= 1, dy= 0, dz= 1; };  
template<> struct dir_of<10> { static constexpr int dx=-1, dy= 0, dz=-1; };
template<> struct dir_of<11> { static constexpr int dx= 0, dy= 1, dz= 1; }; 
template<> struct dir_of<12> { static constexpr int dx= 0, dy=-1, dz=-1; };
template<> struct dir_of<13> { static constexpr int dx= 1, dy=-1, dz= 0; };
template<> struct dir_of<14> { static constexpr int dx=-1, dy= 1, dz= 0; };
template<> struct dir_of<15> { static constexpr int dx= 1, dy= 0, dz=-1; };
template<> struct dir_of<16> { static constexpr int dx=-1, dy= 0, dz= 1; };  
template<> struct dir_of<17> { static constexpr int dx= 0, dy= 1, dz=-1; };
template<> struct dir_of<18> { static constexpr int dx= 0, dy=-1, dz= 1; };  
template<> struct dir_of<19> { static constexpr int dx= 1, dy= 1, dz= 1; };
template<> struct dir_of<20> { static constexpr int dx=-1, dy=-1, dz=-1; };
template<> struct dir_of<21> { static constexpr int dx= 1, dy= 1, dz=-1; };
template<> struct dir_of<22> { static constexpr int dx=-1, dy=-1, dz= 1; };
template<> struct dir_of<23> { static constexpr int dx= 1, dy=-1, dz= 1; }; 
template<> struct dir_of<24> { static constexpr int dx=-1, dy= 1, dz=-1; };
template<> struct dir_of<25> { static constexpr int dx=-1, dy= 1, dz= 1; }; 
template<> struct dir_of<26> { static constexpr int dx= 1, dy=-1, dz=-1; };

// ========================================================================= //

template<int DIR, int GDIR>
__device__ __forceinline__
void emitDirInflowZ(
    LBMFields d,
    const int x, 
    const int y, 
    const int z, 
    const idx_t idx3_in         
) {
    constexpr int dx = dir_of<DIR>::dx;
    constexpr int dy = dir_of<DIR>::dy;
    constexpr int dz = dir_of<DIR>::dz;

    const idx_t nbrIdx = global3(x+dx, y+dy, z+dz);

    float feq = computeFeq(d.rho[nbrIdx], 0.0f, 0.0f, U_OIL, UU_OIL, DIR);
    float fneqReg = computeNeq(
        d.pxx[nbrIdx], d.pyy[nbrIdx], d.pzz[nbrIdx],
        d.pxy[nbrIdx], d.pxz[nbrIdx], d.pyz[nbrIdx],
        d.ux[nbrIdx],  d.uy[nbrIdx],  d.uz[nbrIdx],
        DIR
    );

    d.f[DIR * PLANE + nbrIdx] = to_pop(feq + OMCO_ZMIN * fneqReg);

    if constexpr (DIR == GDIR) {
        feq = computeGeq(d.phi[idx3_in], 0.0f, 0.0f, U_OIL, DIR);
        d.g[DIR * PLANE + nbrIdx] = feq;
    }
}

template<int GDIR, int... DIRS>
__device__ __forceinline__
void emitInflowZ(
    LBMFields d,
    const int x, 
    const int y, 
    const int z,
    const idx_t idx3_in
) {
    (emitDirInflowZ<DIRS, GDIR>(d, x,y,z, idx3_in), ...);
}

// ========================================================================= //

template<int DIR, int GDIR>
__device__ __forceinline__
void emitDirInflowY(
    LBMFields d,
    const int x, 
    const int y, 
    const int z, 
    const idx_t idx3_in      
) {
    constexpr int dx = dir_of<DIR>::dx;
    constexpr int dy = dir_of<DIR>::dy;
    constexpr int dz = dir_of<DIR>::dz;

    const idx_t nbrIdx = global3(x+dx, y+dy, z+dz);

    float feq = computeFeq(d.rho[nbrIdx], 0.0f, U_WATER, 0.0f, UU_WATER, DIR);
    float fneqReg = computeNeq(
        d.pxx[nbrIdx], d.pyy[nbrIdx], d.pzz[nbrIdx],
        d.pxy[nbrIdx], d.pxz[nbrIdx], d.pyz[nbrIdx],
        d.ux[nbrIdx],  d.uy[nbrIdx],  d.uz[nbrIdx],
        DIR
    );

    d.f[DIR * PLANE + nbrIdx] = to_pop(feq + OMCO_YMIN * fneqReg);

    if constexpr (DIR == GDIR) {
        feq = computeGeq(d.phi[idx3_in], 0.0f, U_WATER, 0.0f, DIR);
        d.g[DIR * PLANE + nbrIdx] = feq;
    }
}

template<int GDIR, int... DIRS>
__device__ __forceinline__
void emitInflowY(
    LBMFields d,
    const int x, 
    const int y, 
    const int z,
    const idx_t idx3_in
) {
    (emitDirInflowY<DIRS, GDIR>(d, x,y,z, idx3_in), ...);
}

// ========================================================================= //

template<int DIR, int GDIR>
__device__ __forceinline__
void emitDirOutflowZ(
    LBMFields d,
    const int x, 
    const int y, 
    const int z, 
    const float ux,
    const float uy,
    const float uz
) {
    constexpr int dx = dir_of<DIR>::dx;
    constexpr int dy = dir_of<DIR>::dy;
    constexpr int dz = dir_of<DIR>::dz;

    const idx_t nbrIdx = global3(x+dx, y+dy, z+dz);

    const float uu = ux*ux + uy*uy + uz*uz;

    float feq = computeFeq(d.rho[nbrIdx], ux, uy, uz, uu, DIR);
    float fneqReg = computeNeq(
        d.pxx[nbrIdx], d.pyy[nbrIdx], d.pzz[nbrIdx],
        d.pxy[nbrIdx], d.pxz[nbrIdx], d.pyz[nbrIdx],
        d.ux[nbrIdx],  d.uy[nbrIdx],  d.uz[nbrIdx],
        DIR
    );

    d.f[DIR * PLANE + nbrIdx] = to_pop(feq + OMCO_MAX * fneqReg);

    if constexpr (DIR == GDIR) {
        feq = computeGeq(d.phi[nbrIdx], ux, uy, uz, DIR);
        d.g[DIR * PLANE + nbrIdx] = feq;
    }
}

template<int GDIR, int... DIRS>
__device__ __forceinline__
void emitOutflowZ(
    LBMFields d,
    const int x, 
    const int y, 
    const int z,
    const float ux,
    const float uy,
    const float uz
) {
    (emitDirOutflowZ<DIRS, GDIR>(d, x,y,z, ux, uy, uz), ...);
}

// ========================================================================= //

template<int DIR, int GDIR>
__device__ __forceinline__
void emitDirOutflowY(
    LBMFields d,
    const int x, 
    const int y, 
    const int z, 
    const float ux,
    const float uy,
    const float uz
) {
    constexpr int dx = dir_of<DIR>::dx;
    constexpr int dy = dir_of<DIR>::dy;
    constexpr int dz = dir_of<DIR>::dz;

    const idx_t nbrIdx = global3(x+dx, y+dy, z+dz);

    const float uu = ux*ux + uy*uy + uz*uz;

    float feq = computeFeq(d.rho[nbrIdx], ux, uy, uz, uu, DIR);
    float fneqReg = computeNeq(
        d.pxx[nbrIdx], d.pyy[nbrIdx], d.pzz[nbrIdx],
        d.pxy[nbrIdx], d.pxz[nbrIdx], d.pyz[nbrIdx],
        d.ux[nbrIdx],  d.uy[nbrIdx],  d.uz[nbrIdx],
        DIR
    );

    d.f[DIR * PLANE + nbrIdx] = to_pop(feq + OMCO_MAX * fneqReg);

    if constexpr (DIR == GDIR) {
        feq = computeGeq(d.phi[nbrIdx], ux, uy, uz, DIR);
        d.g[DIR * PLANE + nbrIdx] = feq;
    }
}

template<int GDIR, int... DIRS>
__device__ __forceinline__
void emitOutflowY(
    LBMFields d,
    const int x, 
    const int y, 
    const int z,
    const float ux,
    const float uy,
    const float uz
) {
    (emitDirOutflowY<DIRS, GDIR>(d, x,y,z, ux, uy, uz), ...);
}

// ========================================================================= //


