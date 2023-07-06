// %%EXPORT_KWD CLIFFORD12_API

// %%GEN h

/** @file clifford12.h
 File ``clifford.h`` is the header file for shared library 
 ``mmgroup_clifford12``. This comprises the following C modules:

 {0}.
*/

#ifndef CLIFFORD12_H
#define CLIFFORD12_H

/// @cond DO_NOT_DOCUMENT 

#include <stdint.h>

#define CLIFFORD12_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define CLIFFORD12_HELPER_DLL_IMPORT __declspec(dllimport)
  #define CLIFFORD12_HELPER_DLL_EXPORT __declspec(dllexport)
  #define CLIFFORD12_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define CLIFFORD12_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define CLIFFORD12_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define CLIFFORD12_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define CLIFFORD12_HELPER_DLL_IMPORT
    #define CLIFFORD12_HELPER_DLL_EXPORT
    #define CLIFFORD12_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define CLIFFORD12_API 
// and CLIFFORD12_LOCAL.
// CLIFFORD12_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// CLIFFORD12_LOCAL is used for non-api symbols.

#ifdef CLIFFORD12_DLL // defined if CLIFFORD12 is compiled as a DLL
  #ifdef CLIFFORD12_DLL_EXPORTS // defined if we are building the CLIFFORD12 DLL 
                           // (instead of using it)
    #define CLIFFORD12_API CLIFFORD12_HELPER_DLL_EXPORT
  #else
    #define CLIFFORD12_API CLIFFORD12_HELPER_DLL_IMPORT
  #endif // CLIFFORD12_DLL_EXPORTS
  #define CLIFFORD12_LOCAL CLIFFORD12_HELPER_DLL_LOCAL
#else // CLIFFORD12_DLL is not defined: this means CLIFFORD12 is a static lib.
  #define CLIFFORD12_API
  #define CLIFFORD12_LOCAL
#endif // CLIFFORD12_DLL


/// @endcond 


// %%INCLUDE_HEADERS


// %%GEN h

#ifdef __cplusplus
}
#endif
#endif  // #ifndef CLIFFORD12_H
