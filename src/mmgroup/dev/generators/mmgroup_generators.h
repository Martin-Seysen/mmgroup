// %%GEN h
#ifndef MMGROUP_GENERATORS_H
#define MMGROUP_GENERATORS_H
// %%GEN c


// %%GEN h
/** @file mmgroup_generators.h

 The header file ``mmgroup_generators.h`` contains definitions for 
 the C files in the ``generator`` extension. This extension comprises 
 files ``mm_group_n.c``, ``gen_xi_functions.c``, and ``gen_leech.c``.

 In this header we also define an ``enum MMGROUP_ATOM_TAG_`` that
 specifies the format of an atom that acts as a generator of the
 monster group. 
*/










/** 
 In this header file we also define the tags for the atoms generating 
 the Monster group. An element of the monster group is represented
 as an array of integers of type ``uint32_t``, where each integer 
 represents an atom, i.e. an atomic element of the monster. An atom
 represents a triple ``(sign, tag, value)``, and is encoded in the 
 following bit fields of an unsigned 32-bit integer:

       Bit 31  | Bit 30,...,28  | Bit  27,...,0
       --------|----------------|----------------
        Sign   | Tag            | Value

 Standard tags and values are defined as in the constructor of the 
 Python class ``mmgroup.mm``, see section **The monster group**
 in the **API reference**. If the ``sign`` bit is set, this means
 that the atom bit given by the pair ``(tag, value)`` has to be 
 inverted. In ibid., a tag is given by a small letter. These small
 letters are converted to 3-bit numbers as follows:

       Tag  | Tag number | Range of possible values i
       -----|------------|----------------------------
        'd' |     1      |  0 <= i < 0x1000
        'p' |     2      |  0 <= i < 244823040
        'x' |     3      |  0 <= i < 0x2000
        'y' |     4      |  0 <= i < 0x2000
        't' |     5      |  0 <= i < 3
        'l' |     6      |  0 <= i < 3

 A tag with tag number 0 is interpreted as the neutral element.
 A tag with tag number 7 is illegal (and reserved for future use).
 
 Tags with other letters occuring in the constructor of class ``MM`` 
 are converted to a word of atoms with tags taken from the table 
 above.
*/

enum MMGROUP_ATOM_TAG_ {
  /** Tag indicating the neutral element of the group */
  MMGROUP_ATOM_TAG_1  = 0x00000000UL,
  /** Tag indicating the neutral element of the group */
  MMGROUP_ATOM_TAG_I1  = 0x80000000UL,
  /** Tag corresponding to 'd' */
  MMGROUP_ATOM_TAG_D   = 0x10000000UL,
  /** Tag corresponding to inverse of tag 'd' */
  MMGROUP_ATOM_TAG_ID   = 0x90000000UL,
  /** Tag corresponding to 'p' */
  MMGROUP_ATOM_TAG_P   = 0x20000000UL,
  /** Tag corresponding to inverse of tag 'p' */
  MMGROUP_ATOM_TAG_IP   = 0xA0000000UL,
  /** Tag corresponding to 'x' */
  MMGROUP_ATOM_TAG_X   = 0x30000000UL,
  /** Tag corresponding to inverse of tag 'x' */
  MMGROUP_ATOM_TAG_IX   = 0xB0000000UL,
  /** Tag corresponding to 'y' */
  MMGROUP_ATOM_TAG_Y   = 0x40000000UL,
  /** Tag corresponding to inverse of tag 'y' */
  MMGROUP_ATOM_TAG_IY   = 0xC0000000UL,
  /** Tag corresponding to 't' */
  MMGROUP_ATOM_TAG_T   = 0x50000000UL,
  /** Tag corresponding to inverse of tag 't' */
  MMGROUP_ATOM_TAG_IT   = 0xD0000000UL,
  /** Tag corresponding to 'l' */
  MMGROUP_ATOM_TAG_L   = 0x60000000UL,
  /** Tag corresponding to inverse of tag 'l' */
  MMGROUP_ATOM_TAG_IL   = 0xE0000000UL,
};


/** Tag field of a monster group atom  */
#define MMGROUP_ATOM_TAG_ALL 0xF0000000UL


/** Data field of a monster group atom  */
#define MMGROUP_ATOM_DATA 0xFFFFFFFUL



// %%INCLUDE_HEADERS


// %%GEN h
#endif // ifndef MMGROUP_GENERATORS_H
// %%GEN c

