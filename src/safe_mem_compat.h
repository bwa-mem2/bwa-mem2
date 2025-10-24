/*************************************************************************************
 ARM/macOS compatibility wrapper for safestringlib's safe_mem_lib.h

 - On macOS, the system provides memset_s with a different signature.
   To avoid prototype conflicts when including safe_mem_lib.h, we remap
   memset_s to a private name for the purpose of including the header.

 This header should be included instead of <safe_mem_lib.h> within this project.
*************************************************************************************/
#ifndef SAFE_MEM_COMPAT_H
#define SAFE_MEM_COMPAT_H

#if defined(__APPLE__) && defined(__MACH__)
  // Remap memset_s symbol in the header to avoid conflicting declaration
  #define memset_s memset_s_safeclib_header_only
#endif

#include "safe_mem_lib.h"

#endif // SAFE_MEM_COMPAT_H
