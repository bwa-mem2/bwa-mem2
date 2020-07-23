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

Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
                                Heng Li <hli@jimmy.harvard.edu> 
*****************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <limits.h>
#include <assert.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

#define SIMD_SSE     0x1
#define SIMD_SSE2    0x2
#define SIMD_SSE3    0x4
#define SIMD_SSSE3   0x8
#define SIMD_SSE4_1  0x10
#define SIMD_SSE4_2  0x20
#define SIMD_AVX     0x40
#define SIMD_AVX2    0x80
#define SIMD_AVX512F 0x100
#define SIMD_AVX512BW 0x200

#ifndef _MSC_VER
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void __cpuidex(int cpuid[4], int func_id, int subfunc_id)
{
#if defined(__x86_64__)
	__asm__ volatile ("cpuid"
			: "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#else // on 32bit, ebx can NOT be used as PIC code
	__asm__ volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
			: "=a" (cpuid[0]), "=r" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#endif
}
#endif

static int x86_simd(void)
{
	int flag = 0, cpuid[4], max_id;
	__cpuidex(cpuid, 0, 0);
	max_id = cpuid[0];
	if (max_id == 0) return 0;
	__cpuidex(cpuid, 1, 0);
	if (cpuid[3]>>25&1) flag |= SIMD_SSE;
	if (cpuid[3]>>26&1) flag |= SIMD_SSE2;
	if (cpuid[2]>>0 &1) flag |= SIMD_SSE3;
	if (cpuid[2]>>9 &1) flag |= SIMD_SSSE3;
	if (cpuid[2]>>19&1) flag |= SIMD_SSE4_1;
	if (cpuid[2]>>20&1) flag |= SIMD_SSE4_2;
	if (cpuid[2]>>28&1) flag |= SIMD_AVX;
	if (max_id >= 7) {
		__cpuidex(cpuid, 7, 0);
		if (cpuid[1]>>5 &1) flag |= SIMD_AVX2;
		if (cpuid[1]>>16&1) flag |= SIMD_AVX512F;
		if (cpuid[1]>>30&1) flag |= SIMD_AVX512BW;
	}
	return flag;
}

static int exe_path(const char *exe, int max, char buf[], int *base_st)
{
	int i, len, last_slash, ret = 0;
	if (exe == 0 || max == 0) return -1;
	buf[0] = 0;
	len = strlen(exe);
	for (i = len - 1; i >= 0; --i)
		if (exe[i] == '/') break;
	last_slash = i;
	if (base_st) *base_st = last_slash + 1;
	if (exe[0] == '/') {
		if (max < last_slash + 2) return -1;
		strncpy_s(buf, max, exe, last_slash + 1);
		buf[last_slash + 1] = 0;
	} else if (last_slash >= 0) { // actually, can't be 0
		char *p;
		//int abs_len;
		p = getcwd(buf, max);
		if (p == 0) return -1;
		//abs_len = strlen(buf);
		//if (max < abs_len + 3 + last_slash) return -1;
		//buf[abs_len] = '/';
        strcat_s(buf, max, "/");
		//strncpy_s(buf + abs_len + 1, max - (abs_len + 1), exe, last_slash + 1);
        strncat_s(buf, max, exe, last_slash + 1);
		//buf[abs_len + last_slash + 2] = 0;
	} else {
		char *env, *p, *q, *tmp;
		int env_len, found = 0;
		struct stat st;
		env = getenv("PATH");
        assert(env != NULL);
		env_len = strlen(env);
		if ((tmp = (char*)malloc(env_len + len + 2)) == NULL) { fprintf( stderr, "ERROR: out of memory %s", __func__); exit(EXIT_FAILURE);}
		for (p = q = env;; ++p) {
			if (*p == ':' || *p == 0) {
				strncpy_s(tmp, env_len + len + 2, q, p - q);
				tmp[p - q] = '/';
				strcpy_s(tmp + (p - q + 1), env_len + len + 2 - (p - q + 1), exe);
				if (stat(tmp, &st) == 0 && (st.st_mode & S_IXUSR)) {
					found = 1;
					break;
				}
				if (*p == 0) break;
				q = p + 1;
			}
		}
		if (!found) {
			free(tmp);
			return -2; // shouldn't happen!
		}
		ret = exe_path(tmp, max, buf, 0);
		free(tmp);
	}
	return ret;
}

static void test_and_launch(char *argv[], char *prefix, const char *simd) // we assume prefix is long enough
{
	struct stat st;
    int prefix_len = strlen(prefix);
	strcat_s(prefix, PATH_MAX, simd);
    fprintf(stderr, "Looking to launch executable \"%s\", simd = %s\n", prefix, simd);
    fprintf(stderr, "stat(prefix, &st) = %d\n", stat(prefix, &st));
    if(stat(prefix, &st) == 0)
    {
        fprintf(stderr, "(st.st_mode & S_IXUSR) = %d\n", (st.st_mode & S_IXUSR));
        if ((st.st_mode & S_IXUSR)) {
            fprintf(stderr, "Launching executable \"%s\"\n", prefix);
            execv(prefix, argv);
        }
        else
        {
            fprintf(stderr, "Can not run executable: %s\n", prefix);
        }
    }
    else
    {
        fprintf(stderr, "Can not run executable: %s\n", prefix);
    }
    prefix[prefix_len] = 0;
}

int main(int argc, char *argv[])
{
	char buf[PATH_MAX], *prefix, *argv0 = argv[0];
	int ret, base_st, simd;
	//int prefix_len;
	ret = exe_path(argv0, PATH_MAX, buf, &base_st);
	if (ret != 0) {
		fprintf(stderr, "ERROR: prefix is too long!\n");
		return 1;
	}
	//printf("%s\n", buf);
	if ((prefix = (char*)malloc(PATH_MAX)) == NULL) {
		fprintf(stderr, "ERROR: out of memory.\n");
		return 1;
        }
	strcpy_s(prefix, PATH_MAX, buf);
	strcat_s(prefix, PATH_MAX, &argv0[base_st]);
	//prefix_len = strlen(prefix);
	simd = x86_simd();
	if (simd & SIMD_AVX512BW) test_and_launch(argv, prefix, ".avx512bw");
	if (simd & SIMD_AVX512F) test_and_launch(argv, prefix, ".avx512f");
	if (simd & SIMD_AVX2) test_and_launch(argv, prefix, ".avx2");
	if (simd & SIMD_AVX) test_and_launch(argv, prefix, ".avx");
	if (simd & SIMD_SSE4_2) test_and_launch(argv, prefix, ".sse42");
	if (simd & SIMD_SSE4_1) test_and_launch(argv, prefix, ".sse41");
	if (simd & SIMD_SSE2) test_and_launch(argv, prefix, ".sse2");
	if (simd & SIMD_SSE) test_and_launch(argv, prefix, ".sse");
	free(prefix);
	fprintf(stderr, "ERROR: fail to find the right executable\n");
	return 2;
}
