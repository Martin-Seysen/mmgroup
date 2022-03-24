// %%GEN h

/** @file mm_basics.h

 The header file ``mm_basics.h`` contains basic definitions for dealing with 
 vectors of the 198884-dimensional representation of the monster group,
 as described in  *The C interface of the mmgroup project*, 
 section *Description of the mmgroup.mm extension*.


 It also contains prototypes for the C files in the ``mm`` extension. This 
 extension comprises the  files ``mm_aux.c``,  ``mm_crt.c``,  
``mm_group_word.c``, ``mm_random.c``, ``mm_tables.c``, 
``mm_tables_xi.c``. 
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




/// @endcond  

/** @var typedef uint%{INT_BITS}_t uint_mmv_t
    @brief Used for the representation of the monster group

    Internally, a vector in the 196884-dimensional representation of 
    the monster is stored as an array of integers of 
    type ``uint_mmv_t``. Here several entries are stored in such
    an integer. See ``enum MM_AUX_OFS_type`` for more details.
*/
typedef uint%{INT_BITS}_t uint_mmv_t;


/**********************************************************************
*** Some macros
**********************************************************************/

// %%GEN h

/** 
  This enumeration contains the  offsets for the tags ``A,B,C,T,X,Z,Y``
  in a vector in the 196884-dimensional representation of the monster,
  stored in the internal representation.

  Such an offset counts the number of entries starting at the 
  beginning of th vector. Note that several entries of a vector are 
  stored in a %{INT_BITS}-bit integer. Also there may be duplicate or
  unused entries in a vector, in order to speed up the operation
  of the monster group on a vector.
*/
enum MM_AUX_OFS {
  MM_AUX_OFS_A =      0UL, /**< Offset for tag A */ 
  MM_AUX_OFS_B =    768UL, /**< Offset for tag B */
  MM_AUX_OFS_C =   1536UL, /**< Offset for tag C */
  MM_AUX_OFS_T =   2304UL, /**< Offset for tag T */
  MM_AUX_OFS_X =  50880UL, /**< Offset for tag X */
  MM_AUX_OFS_Z = 116416UL, /**< Offset for tag Z */
  MM_AUX_OFS_Y = 181952UL, /**< Offset for tag Y */
  MM_AUX_LEN_V = 247488UL  /**< Total length of the internal representation */
};


/** 
  This enumeration contains the offsets for the tags ``A,B,C,T,X,Z,Y``
  in a vector in the 196884-dimensional representation of the monster,
  stored in the external representation.

  In external representation, a vector is stored as a contiguous 
  array of bytes.
*/
enum MM_AUX_XOFS {
  MM_AUX_XOFS_D =      0UL, /**< Offset for diagonal entries of tag A */ 
  MM_AUX_XOFS_A =     24UL, /**< Offset for tag A */ 
  MM_AUX_XOFS_B =    300UL, /**< Offset for tag B */ 
  MM_AUX_XOFS_C =    576UL, /**< Offset for tag C */ 
  MM_AUX_XOFS_T =    852UL, /**< Offset for tag T */ 
  MM_AUX_XOFS_X =  49428UL, /**< Offset for tag X */ 
  MM_AUX_XOFS_Z =  98580UL, /**< Offset for tag Z */ 
  MM_AUX_XOFS_Y = 147732UL, /**< Offset for tag Y */ 
  MM_AUX_XLEN_V = 196884UL  /**< Total length of the external representation */
};


/**
  This enumeration defines the values of the tags ``A,B,C,T,X,Z,Y``
  in a vector in the 196884-dimensional representation of the monster,
  stored in the sparse representation.

  In the sparse representation an entry of a vector is stored as a tuple
  of bit  fields ``(tag, par1, par2, value)`` inside an integer of
  type ``uint32_t`` as follows:

      Bits 27,...,25:  tag (as indicated below)
 
      Bits 24,...,14:  par1 (an integer of up to 11 bits)

      Bits 13,..., 8:  par2 (an integer of up to 6 bits)
  
      Bits  7,..., 0:  value (Reserved for the value of an entry)
*/
enum MM_SPACE_TAG {
  MM_SPACE_TAG_A =    0x2000000UL, /**< Encodes tag A */ 
  MM_SPACE_TAG_B =    0x4000000UL, /**< Encodes tag B */
  MM_SPACE_TAG_C =    0x6000000UL, /**< Encodes tag C */
  MM_SPACE_TAG_T =    0x8000000UL, /**< Encodes tag T */
  MM_SPACE_TAG_X =    0xA000000UL, /**< Encodes tag X */
  MM_SPACE_TAG_Z =    0xC000000UL, /**< Encodes tag Z */
  MM_SPACE_TAG_Y =    0xE000000UL  /**< Encodes tag Y */
};



/** @def mm_aux_bad_p(p)
    @brief Return 0 if ``p`` is a good modulus and a nonzero value otherwise 
*/
#define mm_aux_bad_p(p) (((p) & ((p)+1)) | (((p)-3) & ((0UL-256UL))))


/// @cond DO_NOT_DOCUMENT 

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

/// @endcond

