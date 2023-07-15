// %%GEN h



/** 
  This enumeration contains the  offsets for the tags ``A,B,C,T,X,Z,Y``
  in a vector in the 196884-dimensional representation of the monster,
  stored in the internal representation.

  This is similar to enum MM_AUX_OFS in file ``mm_basics.h``. But 
  here the offsets are given in units of %{INT_BITS}-bit integers
  for a vector of the  representation \f$\rho_{%{P}}\f$ of the
  monster group in characteristic  %{P}.
*/
enum MM_OP%{P}_OFS  {
 MM_OP%{P}_OFS_A = (MM_AUX_OFS_A >> %{LOG_INT_FIELDS}), /**< Offset for tag A */
 MM_OP%{P}_OFS_B = (MM_AUX_OFS_B >> %{LOG_INT_FIELDS}), /**< Offset for tag B */   
 MM_OP%{P}_OFS_C = (MM_AUX_OFS_C >> %{LOG_INT_FIELDS}), /**< Offset for tag C */    
 MM_OP%{P}_OFS_T = (MM_AUX_OFS_T >> %{LOG_INT_FIELDS}), /**< Offset for tag T */  
 MM_OP%{P}_OFS_X = (MM_AUX_OFS_X >> %{LOG_INT_FIELDS}), /**< Offset for tag X */   
 MM_OP%{P}_OFS_Z = (MM_AUX_OFS_Z >> %{LOG_INT_FIELDS}), /**< Offset for tag Z */   
 MM_OP%{P}_OFS_Y = (MM_AUX_OFS_Y >> %{LOG_INT_FIELDS}), /**< Offset for tag Y */    
 MM_OP%{P}_LEN_V = (MM_AUX_LEN_V >> %{LOG_INT_FIELDS}), /**< Total length of the internal representation */    
};
                                  
