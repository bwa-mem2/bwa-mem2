/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef LH3_UTILS_H
#define LH3_UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#ifdef __GNUC__
// Tell GCC to validate printf format string and args
#define ATTRIBUTE(list) __attribute__ (list)
#else
#define ATTRIBUTE(list)
#endif

#define err_fatal_simple(msg) _err_fatal_simple(__func__, msg)
#define err_fatal_simple_core(msg) _err_fatal_simple_core(__func__, msg)

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)

#define xassert(cond, msg) if ((cond) == 0) _err_fatal_simple_core(__func__, msg)

// clang and gcc 11+ define __rdtsc intrinsics on x86/64
#if defined(__GNUC__) && __GNUC__ < 11 && !defined(__clang__)
#if defined(__i386__)
static inline unsigned long long __rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#elif defined(__x86_64__)
static inline unsigned long long __rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif
#endif
// From https://github.com/google/benchmark/blob/37177a84b7e8d33696ea1e1854513cb0de3b4dc3/src/cycleclock.h
// Apache 2.0 license
#if defined(__aarch64__)
  // System timer of ARMv8 runs at a different frequency than the CPU's.
  // The frequency is fixed, typically in the range 1-50MHz.  It can be
  // read at CNTFRQ special register.  We assume the OS has set up
  // the virtual timer properly.
static inline unsigned long long __rdtsc(void)
{
  int64_t virtual_timer_value;
  asm volatile("mrs %0, cntvct_el0" : "=r"(virtual_timer_value));
  return virtual_timer_value;
}
#elif defined(__ARM_ARCH)
  // V6 is the earliest arch that has a standard cyclecount
  // Native Client validator doesn't allow MRC instructions.
#if (__ARM_ARCH >= 6)
static inline unsigned long long __rdtsc(void)
{
  uint32_t pmccntr;
  uint32_t pmuseren;
  uint32_t pmcntenset;
  // Read the user mode perf monitor counter access permissions.
  asm volatile("mrc p15, 0, %0, c9, c14, 0" : "=r"(pmuseren));
  if (pmuseren & 1) {  // Allows reading perfmon counters for user mode code.
    asm volatile("mrc p15, 0, %0, c9, c12, 1" : "=r"(pmcntenset));
    if (pmcntenset & 0x80000000ul) {  // Is it counting?
      asm volatile("mrc p15, 0, %0, c9, c13, 0" : "=r"(pmccntr));
      // The counter is set up to count every 64th cycle
      return static_cast<int64_t>(pmccntr) * 64;  // Should optimize to << 6
    }
  }
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return static_cast<int64_t>(tv.tv_sec) * 1000000 + tv.tv_usec;
}
#else
#error __ARM_ARCH < 6 does not have a standard cyclecount
#endif
#endif

typedef struct {
	uint64_t x, y;
} pair64_t;

typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; pair64_t *a; } pair64_v;

#ifdef __cplusplus
extern "C" {
#endif

	void err_fatal(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void err_fatal_core(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void _err_fatal_simple(const char *func, const char *msg) ATTRIBUTE((noreturn));
	void _err_fatal_simple_core(const char *func, const char *msg) ATTRIBUTE((noreturn));
	FILE *err_xopen_core(const char *func, const char *fn, const char *mode);
	FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp);
	gzFile err_xzopen_core(const char *func, const char *fn, const char *mode);
	size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
	size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream);

	int err_gzread(gzFile file, void *ptr, unsigned int len);
	int err_fseek(FILE *stream, long offset, int whence);
#define err_rewind(FP) err_fseek((FP), 0, SEEK_SET)
	long err_ftell(FILE *stream);
	int err_fprintf(FILE *stream, const char *format, ...)
        ATTRIBUTE((format(printf, 2, 3)));
	int err_printf(const char *format, ...)
        ATTRIBUTE((format(printf, 1, 2)));
	int err_fputc(int c, FILE *stream);
#define err_putchar(C) err_fputc((C), stdout)
	char* err_fgets(char *str, int size, FILE *stream);
	int err_fputs(const char *s, FILE *stream);
	int err_puts(const char *s);
	int err_fflush(FILE *stream);
	int err_fclose(FILE *stream);
	int err_gzclose(gzFile file);

	double cputime();
	double realtime();

	void ks_introsort_64 (size_t n, uint64_t *a);
	void ks_introsort_128(size_t n, pair64_t *a);

#ifdef __cplusplus
}
#endif

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

#endif
