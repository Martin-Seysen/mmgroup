// %%GEN h
#ifndef MAT24_FUNCTIONS_H
#define MAT24_FUNCTIONS_H

#include <stdint.h>

/** @file mat24_functions.h
  File ``mat24_functions.h`` is the header file for ``mat24_functions.c``. 
*/

/// @cond DO_NOT_DOCUMENT 

/**************************************************************
#define MAT24_DLL  // We want a DLL!!

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define MAT24_HELPER_DLL_IMPORT __declspec(dllimport)
  #define MAT24_HELPER_DLL_EXPORT __declspec(dllexport)
  #define MAT24_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define MAT24_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define MAT24_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define MAT24_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define MAT24_HELPER_DLL_IMPORT
    #define MAT24_HELPER_DLL_EXPORT
    #define MAT24_HELPER_DLL_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define MAT24_API 
// and MAT24_LOCAL.
// MAT24_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// MAT24_LOCAL is used for non-api symbols.

#ifdef MAT24_DLL // defined if MAT24 is compiled as a DLL
  #ifdef MAT24_DLL_EXPORTS // defined if we are building the MAT24 DLL 
                           // (instead of using it)
    #define MAT24_API MAT24_HELPER_DLL_EXPORT
  #else
    #define MAT24_API MAT24_HELPER_DLL_IMPORT
  #endif // MAT24_DLL_EXPORTS
  #define MAT24_LOCAL MAT24_HELPER_DLL_LOCAL
#else // MAT24_DLL is not defined: this means MAT24 is a static lib.
  #define MAT24_API
  #define MAT24_LOCAL
#endif // MAT24_DLL
****************************************************************/

/// @endcond


// %%INCLUDE_HEADERS


// %%GEN h
#endif  //  ifndef MAT24_FUNCTIONS_H
// %%GEN c




