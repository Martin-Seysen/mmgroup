// %%GEN h
// This header has been created automatically, do not edit!

#ifndef MM_OP_H
#define MM_OP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mm_basics.h"


// %%GEN h


/** @file mm_op_sub.h

 The header file ``mm_op_sub.h`` contains basic definitions for the
 C files dealing with  vectors of the 198884-dimensional representation
 of the monster group modulo ``p``,
 as described in  *The C interface of the mmgroup project*,
 section *Description of the mmgroup.mm extension*.
*/

/// @cond DO_NOT_DOCUMENT 


#define MM_OP_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define MM_OP_HELPER_DLL_IMPORT __declspec(dllimport)
  #define MM_OP_HELPER_DLL_EXPORT __declspec(dllexport)
  #define MM_OP_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define MM_OP_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define MM_OP_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define MM_OP_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define MM_OP_HELPER_DLL_IMPORT
    #define MM_OP_HELPER_DLL_EXPORT
    #define MM_OP_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define MM_OP_API 
// and MM_OP_LOCAL.
// MM_OP_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// MM_OP_LOCAL is used for non-api symbols.

#ifdef MM_OP_DLL // defined if MM_BASICS is compiled as a DLL
  #ifdef MM_OP_DLL_EXPORTS // defined if we are building the MM_BASICS DLL 
                           // (instead of using it)
    #define MM_OP_API MM_OP_HELPER_DLL_EXPORT
  #else
    #define MM_OP_API MM_OP_HELPER_DLL_IMPORT
  #endif // MM_OP_DLL_EXPORTS
  #define MM_OP_LOCAL MM_OP_HELPER_DLL_LOCAL
#else // MM_OP_DLL is not defined: this means MM_BASICS is a static lib.
  #define MM_OP_API
  #define MM_OP_LOCAL
#endif // MM_OP_DLL


/// @endcond  


// %%INCLUDE_HEADERS


// %%GEN h
#ifdef __cplusplus
}
#endif
#endif  // #ifndef MM_OP_H


                                  
