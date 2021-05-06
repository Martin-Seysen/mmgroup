// %%GEN h

/** @file mm_basics.h

 This is file mm_basics.h. It is yet to be documented.

*/

#include <stdint.h>
#include "mat24_functions.h"
#include "mmgroup_generators.h"

/// @cond DO_NOT_DOCUMENT 


#define MM_BASICS_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define MM_BASICS_HELPER_DLL_IMPORT __declspec(dllimport)
  #define MM_BASICS_HELPER_DLL_EXPORT __declspec(dllexport)
  #define MM_BASICS_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define MM_BASICS_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define MM_BASICS_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define MM_BASICS_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define MM_BASICS_HELPER_DLL_IMPORT
    #define MM_BASICS_HELPER_DLL_EXPORT
    #define MM_BASICS_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define MM_BASICS_API 
// and MM_BASICS_LOCAL.
// MM_BASICS_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// MM_BASICS_LOCAL is used for non-api symbols.

#ifdef MM_BASICS_DLL // defined if MM_BASICS is compiled as a DLL
  #ifdef MM_BASICS_DLL_EXPORTS // defined if we are building the MM_BASICS DLL 
                           // (instead of using it)
    #define MM_BASICS_API MM_BASICS_HELPER_DLL_EXPORT
  #else
    #define MM_BASICS_API MM_BASICS_HELPER_DLL_IMPORT
  #endif // MM_BASICS_DLL_EXPORTS
  #define MM_BASICS_LOCAL MM_BASICS_HELPER_DLL_LOCAL
#else // MM_BASICS_DLL is not defined: this means MM_BASICS is a static lib.
  #define MM_BASICS_API
  #define MM_BASICS_LOCAL
#endif // MM_BASICS_DLL



typedef uint%{INT_BITS}_t uint_mmv_t;

/// @endcond  


/**********************************************************************
*** Some macros
**********************************************************************/

// %%GEN h
// Return nonzero value if p is a bad modulus,  
// i.e. not p = 2**k - 1 for some 2 <= k <= 8
#define mm_aux_bad_p(p) (((p) & ((p)+1)) | (((p)-3) & ((0UL-256UL))))

// Offsets for tags A,B,C,T,X,Z,Y in the internal representation
#define MM_AUX_OFS_A       0UL
#define MM_AUX_OFS_B     768UL    //    24*32
#define MM_AUX_OFS_C    1536UL    //  2*24*32
#define MM_AUX_OFS_T    2304UL    //  3*24*32
#define MM_AUX_OFS_X   50880UL    //  MM_AUX_OFS_T +    759*64
#define MM_AUX_OFS_Z  116416UL    //  MM_AUX_OFS_X +   2048*32
#define MM_AUX_OFS_Y  181952UL    //  MM_AUX_OFS_X + 2*2048*32
#define MM_AUX_OFS_E  247488UL    //  MM_AUX_OFS_X + 3*2048*32. i.e
                                  //  total length of internal rep

// Offsets for tags A,B,C,T,X,Z,Y in the external representation
#define MM_AUX_XOFS_A      24UL
#define MM_AUX_XOFS_B     300UL    //  24 + 1*276
#define MM_AUX_XOFS_C     576UL    //  24 + 2*276
#define MM_AUX_XOFS_T     852UL    //  24 + 3*276
#define MM_AUX_XOFS_X   49428UL    //  MM_AUX_XOFS_T +    759*64
#define MM_AUX_XOFS_Z   98580UL    //  MM_AUX_XOFS_X +   2048*24
#define MM_AUX_XOFS_Y  147732UL    //  MM_AUX_XOFS_X + 2*2048*24
#define MM_AUX_XOFS_E  196884UL    //  MM_AUX_XOFS_X + 3*2048*24. i.e
                                   //  total length of external rep


// Tags for labels and values of vectors in the representation space
// A multiple of a unit vector with coordinate 'coord' is encoded
// in the bit fields of a 32-bit integers in the form. 
//   coord (tag, par1, par2) 
#define MM_SPACE_TAG_A      0x2000000
#define MM_SPACE_TAG_B      0x4000000
#define MM_SPACE_TAG_C      0x6000000
#define MM_SPACE_TAG_T      0x8000000
#define MM_SPACE_TAG_X      0xA000000
#define MM_SPACE_TAG_Z      0xC000000
#define MM_SPACE_TAG_Y      0xE000000 
// Mask for all tags:
// Use y = (x & MM_SPACE_MASK_PAR1) << MM_SPACE_SHIFT_PAR1
// to set parameter par1 in y to the value x.
#define MM_SPACE_MASK_TAG     0xE000000 
// Mask and shift factor for parameter par1  
// Use y = (x << MM_SPACE_SHIFT_PAR1) & MM_SPACE_MASK_PAR1
// to set parameter par1 in y to the value x.
#define MM_SPACE_MASK_PAR1    0x1FFC000   
#define MM_SPACE_SHIFT_PAR1          14   
// Mask and shift factor for parameter par12  
// Use y = (x << MM_SPACE_SHIFT_PAR2) & MM_SPACE_MASK_PAR2
// to set parameter par2 in y to the value x.
#define MM_SPACE_MASK_PAR2       0x3F00   
#define MM_SPACE_SHIFT_PAR2           8 
// Mask for coordinate:  
// Use y = x  & MM_SPACE_MASK_COORD
// to set the coordiante in y to the value x.
// Caution: some special routines for modulus p = 2**k - 1
// use only the lowest k bits of the coordinate.
#define MM_SPACE_COORD_PAR1    0x1FFC000   



