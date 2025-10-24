/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: ARM compatibility wrapper
*****************************************************************************************/

#ifndef SIMD_COMPAT_H
#define SIMD_COMPAT_H

// Detect ARM architecture
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    #define SSE2NEON_PRECISE_DIV 1
    #define SSE2NEON_PRECISE_SQRT 1
    #include "../ext/sse2neon/sse2neon.h"
    #ifndef __mmask8
        #define __mmask8 uint8_t
    #endif
    #ifndef __mmask16
        #define __mmask16 uint16_t
    #endif
    // Map __rdtsc (x86) to _rdtsc (ARM via sse2neon)
    #define __rdtsc _rdtsc
    // On ARM, override _mm_prefetch to use __builtin_prefetch directly
    // since sse2neon's version might have compatibility issues
    #undef _mm_prefetch
    #define _mm_prefetch(addr, hint) __builtin_prefetch((const char*)(addr), 0, 3)
#else
    // x86/x64 architecture - use native intrinsics
    #if (__AVX512BW__ || __AVX2__)
        #include <immintrin.h>
    #else
        #include <smmintrin.h>  // for SSE4.1
        #define __mmask8 uint8_t
        #define __mmask16 uint16_t
    #endif
#endif

#endif // SIMD_COMPAT_H
