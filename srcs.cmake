set(NTL deps/ntl/src)
set(NTL_SRC
    ${NTL}/tools.cpp
    ${NTL}/lip.cpp
    ${NTL}/WordVector.cpp
    ${NTL}/ZZ.cpp
    ${NTL}/GF2.cpp
    ${NTL}/GF2X.cpp
    ${NTL}/GF2X1.cpp
    ${NTL}/GF2XVec.cpp
    ${NTL}/GF2E.cpp
    ${NTL}/GF2EX.cpp
    ${NTL}/vec_GF2E.cpp
    ${NTL}/mat_GF2E.cpp)
set(NTL_SRC_TEST
    ${NTL}/BasicThreadPool.cpp
    ${NTL}/thread.cpp
    ${NTL}/fileio.cpp
    ${NTL}/GetPID1.cpp
    ${NTL}/ctools.cpp
    ${NTL}/FFT.cpp
    ${NTL}/ZZX.cpp
    ${NTL}/lzz_p.cpp
    ${NTL}/lzz_pX.cpp
    ${NTL}/lzz_pX1.cpp
    ${NTL}/lzz_pE.cpp
    ${NTL}/lzz_pEX.cpp
    ${NTL}/ZZ_p.cpp
    ${NTL}/ZZ_pX.cpp
    ${NTL}/ZZ_pX1.cpp
    ${NTL}/ZZ_pE.cpp
    ${NTL}/ZZ_pEX.cpp
    ${NTL}/ZZVec.cpp
    ${NTL}/vec_GF2.cpp
    ${NTL}/vec_ZZ_p.cpp
    ${NTL}/vec_ZZ_pE.cpp
    ${NTL}/vec_lzz_p.cpp
    ${NTL}/vec_lzz_pE.cpp
    ${NTL}/MatPrime.cpp
    ${NTL}/mat_ZZ_p.cpp
    ${NTL}/mat_ZZ_pE.cpp
    ${NTL}/mat_lzz_p.cpp
    ${NTL}/mat_lzz_pE.cpp
)

set(GF2X deps/gf2x)
set(GF2X_SRC ${GF2X}/gf2x.c ${GF2X}/toom.c ${GF2X}/toom-gpl.c
             ${GF2X}/fft/gf2x-ternary-fft.c)
