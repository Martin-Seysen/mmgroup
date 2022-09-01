// %%GEN h

/** @file mm_reduce.h

 Yet to be documented
*/

#include <stdint.h>
#include <stdlib.h>
#include "mat24_functions.h"
#include "mmgroup_generators.h"
#include "clifford12.h"   
#include "mm_basics.h"   
#include "mm_op3.h"   
#include "mm_op15.h"   
#include "mm_reduce.h"   



/// @cond DO_NOT_DOCUMENT 


#define MM_REDUCE_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define MM_REDUCE_HELPER_DLL_IMPORT __declspec(dllimport)
  #define MM_REDUCE_HELPER_DLL_EXPORT __declspec(dllexport)
  #define MM_REDUCE_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define MM_REDUCE_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define MM_REDUCE_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define MM_REDUCE_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define MM_REDUCE_HELPER_DLL_IMPORT
    #define MM_REDUCE_HELPER_DLL_EXPORT
    #define MM_REDUCE_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define MM_REDUCE_API 
// and MM_REDUCE_LOCAL.
// MM_REDUCE_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// MM_REDUCE_LOCAL is used for non-api symbols.

#ifdef MM_REDUCE_DLL // defined if MM_BASICS is compiled as a DLL
  #ifdef MM_REDUCE_DLL_EXPORTS // defined if we are building a DLL 
                           // (instead of using it)
    #define MM_REDUCE_API MM_REDUCE_HELPER_DLL_EXPORT
  #else
    #define MM_REDUCE_API MM_REDUCE_HELPER_DLL_IMPORT
  #endif // MM_REDUCE_DLL_EXPORTS
  #define MM_REDUCE_LOCAL MM_REDUCE_HELPER_DLL_LOCAL
#else // MM_REDUCE_DLL is not defined: this means we build a static lib.
  #define MM_REDUCE_API
  #define MM_REDUCE_LOCAL
#endif // MM_REDUCE_DLL




/// @endcond  


