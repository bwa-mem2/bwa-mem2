#ifndef SSE2NEON_H
#define SSE2NEON_H

#include <arm_sve.h>

#if defined(__GNUC__) || defined(__clang__)

#pragma push_macro("FORCE_INLINE")
#pragma push_macro("ALIGN_STRUCT")
#define FORCE_INLINE static inline __attribute__((always_inline))
#define ALIGN_STRUCT(x) __attribute__((aligned(x)))

#else

#error "Macro name collisions may happens with unknown compiler"
#ifdef FORCE_INLINE
#undef FORCE_INLINE
#endif
#define FORCE_INLINE static inline
#ifndef ALIGN_STRUCT
#define ALIGN_STRUCT(x) __declspec(align(x))
#endif

#endif

#include <stdint.h>
#include <stdlib.h>
#include <time.h>

typedef svint64_t __m128i;

#define _MM_HINT_NTA 0
#define _MM_HINT_T0 1

FORCE_INLINE void* _mm_malloc(size_t size, size_t align)
{
    //return aligned_alloc(align, size);
    return aligned_alloc(256, size);
}

FORCE_INLINE void _mm_free(void *ptr)
{
    free(ptr);
}

FORCE_INLINE uint64_t __rdtsc()
{
      uint64_t virtual_timer_value;
      asm volatile ("mrs %0, cntvct_el0":"=r" (virtual_timer_value));
      //virtual_timer_value = (unsigned long) time(0);
      return virtual_timer_value;
}

FORCE_INLINE void _mm_prefetch(const char* rseq, int i) {
    __builtin_prefetch((void*) rseq);
}

FORCE_INLINE void _mm_prefetch(const int8_t* rseq, int i) {
    __builtin_prefetch((void*) rseq);
}

FORCE_INLINE void _mm_prefetch(const uint32_t* rseq, int i) {
    __builtin_prefetch((void*) rseq);
}

FORCE_INLINE svint64_t _mm_setzero_si128() {
	return svdup_s64(0);
}

FORCE_INLINE svint64_t _mm_set1_epi8(int8_t data) {
        return svreinterpret_s64(svdup_s8(data));
}

FORCE_INLINE svint64_t _mm_set1_epi16(int16_t data) {
        return svreinterpret_s64(svdup_s16(data));
}

FORCE_INLINE svint64_t _mm_set1_epi32(int32_t data) {
        return svreinterpret_s64(svdup_s32(data));
}

FORCE_INLINE svint64_t _mm_blend_epi8(svint64_t a, svint64_t b, svbool_t mask) {
        svint8_t a_aux = svreinterpret_s8(a);
        svint8_t b_aux = svreinterpret_s8(b);
        svint8_t r_aux = svsel(mask,b_aux,a_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_blend_epi16(svint64_t a, svint64_t b, svbool_t mask) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svsel(mask,b_aux,a_aux);
        return svreinterpret_s64(r_aux);
}

// ----- ARITHMETIC OPS -----

FORCE_INLINE svint64_t _mm_add_epi8(svint64_t a, svint64_t b) {
        svint8_t a_aux = svreinterpret_s8(a);
        svint8_t b_aux = svreinterpret_s8(b);
        svint8_t r_aux = svadd_x(svptrue_b8(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_add_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svadd_x(svptrue_b16(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_adds_epu8(svint64_t a, svint64_t b) {
        svuint8_t a_aux = svreinterpret_u8(a);
        svuint8_t b_aux = svreinterpret_u8(b);
        svuint8_t r_aux = svqadd(a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_adds_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svqadd(a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_sub_epi8(svint64_t a, svint64_t b) {
        svint8_t a_aux = svreinterpret_s8(a);
        svint8_t b_aux = svreinterpret_s8(b);
        svint8_t r_aux = svsub_x(svptrue_b8(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_sub_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svsub_x(svptrue_b16(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_subs_epu8(svint64_t a, svint64_t b) {
        svuint8_t a_aux = svreinterpret_u8(a);
        svuint8_t b_aux = svreinterpret_u8(b);
        svuint8_t r_aux = svqsub(a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_subs_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svqsub(a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_subs_epu16(svint64_t a, svint64_t b) {
        svuint16_t a_aux = svreinterpret_u16(a);
        svuint16_t b_aux = svreinterpret_u16(b);
        svuint16_t r_aux = svqsub(a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_max_epu8(svint64_t a, svint64_t b) {
        svuint8_t a_aux = svreinterpret_u8(a);
        svuint8_t b_aux = svreinterpret_u8(b);
        svuint8_t r_aux = svmax_x(svptrue_b8(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_max_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svmax_x(svptrue_b16(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_max_epu16(svint64_t a, svint64_t b) {
        svuint16_t a_aux = svreinterpret_u16(a);
        svuint16_t b_aux = svreinterpret_u16(b);
        svuint16_t r_aux = svmax_x(svptrue_b16(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_min_epu8(svint64_t a, svint64_t b) {
        svuint8_t a_aux = svreinterpret_u8(a);
        svuint8_t b_aux = svreinterpret_u8(b);
        svuint8_t r_aux = svmin_x(svptrue_b8(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_min_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        svint16_t r_aux = svmin_x(svptrue_b16(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

FORCE_INLINE svint64_t _mm_min_epu16(svint64_t a, svint64_t b) {
        svuint16_t a_aux = svreinterpret_u16(a);
        svuint16_t b_aux = svreinterpret_u16(b);
        svuint16_t r_aux = svmin_x(svptrue_b16(),a_aux,b_aux);
        return svreinterpret_s64(r_aux);
}

// ----- BIT WISE OPS -----

FORCE_INLINE svbool_t _mm_and_si128(svbool_t a, svbool_t b) {
        return svand_z(svptrue_b8(),a,b);
}

FORCE_INLINE svint64_t _mm_and_si128(svint64_t a, svint64_t b) {
        return svand_x(svptrue_b64(),a,b);
}

FORCE_INLINE svbool_t _mm_or_si128(svbool_t a, svbool_t b) {
        return svorr_z(svptrue_b8(),a,b);
}

FORCE_INLINE svint64_t _mm_or_si128(svint64_t a, svint64_t b) {
        return svorr_x(svptrue_b64(),a,b);
}

FORCE_INLINE svint64_t _mm_xor_si128(svint64_t a, svint64_t b) {
        return sveor_x(svptrue_b64(),a,b);
}

FORCE_INLINE svint64_t _mm_andnot_si128(svint64_t a, svint64_t b) {
        return svbic_x(svptrue_b64(),b,a);
}

FORCE_INLINE svbool_t _mm_andnot_si128(svbool_t a, svbool_t b) {
        return svbic_z(svptrue_b8(),b,a);
}

// ----- CMP OPS -----

FORCE_INLINE svbool_t _mm_cmpeq_epi8(svint64_t a, svint64_t b) {
        svint8_t a_aux = svreinterpret_s8(a);
        svint8_t b_aux = svreinterpret_s8(b);
        return svcmpeq(svptrue_b8(),a_aux,b_aux);
}

FORCE_INLINE svbool_t _mm_cmpeq_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        return svcmpeq(svptrue_b16(),a_aux,b_aux);
}

FORCE_INLINE svbool_t _mm_cmpgt_epi8(svint64_t a, svint64_t b) {
        svint8_t a_aux = svreinterpret_s8(a);
        svint8_t b_aux = svreinterpret_s8(b);
        return svcmpgt(svptrue_b8(),a_aux,b_aux);
}

FORCE_INLINE svbool_t _mm_cmpgt_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        return svcmpgt(svptrue_b16(),a_aux,b_aux);
}

FORCE_INLINE svbool_t _mm_cmpge_epi16(svint64_t a, svint64_t b) {
        svint16_t a_aux = svreinterpret_s16(a);
        svint16_t b_aux = svreinterpret_s16(b);
        return svcmpge(svptrue_b16(),a_aux,b_aux);
}

// ----- MEM OPS -----

FORCE_INLINE svint64_t _mm_load_si128(svint64_t * dir) {
	return svld1(svptrue_b64(),(int64_t*)dir);
}

FORCE_INLINE void _mm_store_si128(svint64_t * dir, svint64_t reg) {
	svst1(svptrue_b64(),(int64_t*)dir,reg);
}

#if defined(__GNUC__) || defined(__clang__)
#pragma pop_macro("ALIGN_STRUCT")
#pragma pop_macro("FORCE_INLINE")
#endif

#endif

