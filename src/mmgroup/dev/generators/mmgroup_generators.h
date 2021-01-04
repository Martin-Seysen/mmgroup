// %%GEN h
/** @file mmgroup_generators.h
*/





/// @cond DO_NOT_DOCUMENT 

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

/// @endcond

/**
 * \defgroup  mmgroup_generators_h Definitions in file mmgroup_generators.h
 * @{
 *
 * The header file ``mmgroup_generators.h`` contains definitions for 
 * the C files in the ``generator`` extension. This extension comprises 
 * files ``mm_group_n.c``, ``gen_xi_functions.c``, and ``gen_leech.c``.
 * 
 * In this header file we also define the tags for the atoms generating 
 * the Monster group. An element of the monster group is represented
 * as an array of integers of type uint32_t, where each integer 
 * represents an atom, i.e. an atomic element of the monster.
 * Bits 31,...,28 of a atom are the tag of that atom, as defined below;
 * bits 27,...,0 make up the data part of that atom.
 * 
 * Each tag is represented as a single lower-case letter as described
 * in section *The Monster group* in the *API reference*. Valid letters
 * for tags used in this header are ``'pdxytl'``. If bit 31 of an atom
 * is set, this means that the atom has to be inversed.
 */


/** Tag field of a monster group atom  */
#define MMGROUP_ATOM_TAG_ALL 0xF0000000UL
/** Tag indicating the neutral element of the group */
#define MMGROUP_ATOM_TAG_1 0x00000000UL
/** Tag indicating the neutral element of the group */
#define MMGROUP_ATOM_TAG_I1 0x80000000UL
/** Tag corresponding to 'd' */
#define MMGROUP_ATOM_TAG_D  0x10000000UL
/** Tag corresponding to inverse of tag 'd' */
#define MMGROUP_ATOM_TAG_ID  0x90000000UL
/** Tag corresponding to 'p' */
#define MMGROUP_ATOM_TAG_P  0x20000000UL
/** Tag corresponding to inverse of tag 'p' */
#define MMGROUP_ATOM_TAG_IP  0xA0000000UL
/** Tag corresponding to 'x' */
#define MMGROUP_ATOM_TAG_X  0x30000000UL
/** Tag corresponding to inverse of tag 'x' */
#define MMGROUP_ATOM_TAG_IX  0xB0000000UL
/** Tag corresponding to 'y' */
#define MMGROUP_ATOM_TAG_Y  0x40000000UL
/** Tag corresponding to inverse of tag 'y' */
#define MMGROUP_ATOM_TAG_IY  0xC0000000UL
/** Tag corresponding to 't' */
#define MMGROUP_ATOM_TAG_T  0x50000000UL
/** Tag corresponding to inverse of tag 't' */
#define MMGROUP_ATOM_TAG_IT  0xD0000000UL
/** Tag corresponding to 'l' */
#define MMGROUP_ATOM_TAG_L  0x60000000UL
/** Tag corresponding to inverse of tag 'l' */
#define MMGROUP_ATOM_TAG_IL  0xE0000000UL

/** Data field of a monster group atom  */
#define MMGROUP_ATOM_DATA 0xFFFFFFFUL


/**@}*/

