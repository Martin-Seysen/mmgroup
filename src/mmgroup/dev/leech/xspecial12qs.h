#ifndef XSPECIAL12QS_H
#define XSPECIAL12QS_H


#define QBSTATE12_MAXCOLS     (64)
#define QBSTATE12_MAXROWS     (QBSTATE12_MAXCOLS+1)


typedef struct {
    uint32_t maxrows;
    uint32_t nrows;
    uint32_t ncols;
    int32_t  factor;
    uint64_t *data;
} qbstate12_type;

int32_t qbstate12_set_mem(qbstate12_type *pqs, uint64_t *data, uint32_t size);
int32_t qbstate12_zero(qbstate12_type *pqs, uint32_t nqb);
int32_t qbstate12_standard_state(qbstate12_type *pqs, uint32_t nqb);
int32_t qbstate12_vector_state(qbstate12_type *pqs, uint32_t nqb, uint64_t *pv);
int32_t qbstate12_copy(qbstate12_type *pqs1, qbstate12_type *pqs2); 
int32_t qbstate12_check(qbstate12_type *pqs);
int32_t qbstate12_get_cols(qbstate12_type *pqs, uint32_t c0, uint32_t len, uint32_t *pa);
int32_t qbstate12_set_cols(qbstate12_type *pqs, uint32_t c0, uint32_t len, uint32_t *pa);
int32_t qbstate12_xch_rows(qbstate12_type *pqs, uint32_t i1, uint32_t i2);
int32_t qbstate12_del_rows(qbstate12_type *pqs, uint64_t *pv);
int32_t qbstate12_insert_rows(qbstate12_type *pqs, uint32_t i, uint32_t nrows);
int32_t qbstate12_row_op(qbstate12_type *pqs, uint32_t i1, uint32_t i2); 
int32_t qbstate12_find_pivot(qbstate12_type *pqs, uint32_t j, uint64_t *pv);
int32_t qbstate12_pivot(qbstate12_type *pqs, uint32_t i, uint64_t *pv);
int32_t qbstate12_reduce(qbstate12_type *pqs);
int32_t qbstate12_set_ref_point(qbstate12_type *pqs, uint32_t v, uint32_t ncols);
int32_t qbstate12_gate_x(qbstate12_type *pqs, uint32_t j);
int32_t qbstate12_gate_cx(qbstate12_type *pqs, uint32_t j, uint32_t j1);
int32_t qbstate12_gate_z(qbstate12_type *pqs, uint32_t j);
int32_t qbstate12_gate_cz(qbstate12_type *pqs, uint32_t j1, uint32_t j2);
int32_t qbstate12_gate_h(qbstate12_type *pqs, uint32_t j);
int32_t qbstate12_gate_h_preserve(qbstate12_type *pqs, uint32_t j);
int32_t qbstate12_gate_xch(qbstate12_type *pqs, uint32_t j1, uint32_t j2);
int32_t qbstate12_prep_mul(qbstate12_type *pqs1, qbstate12_type *pqs2, uint32_t nqb);
int32_t qbstate12_reindex(qbstate12_type *pqs, uint32_t n_del, uint32_t j, uint32_t n_ins);
int32_t qbstate12_rot_index(qbstate12_type *pqs, int32_t rot, uint32_t nrot, uint32_t n0);
int32_t qbstate12_mul_elements(qbstate12_type *pqs1, qbstate12_type *pqs2);
int32_t qbstate12_contract(qbstate12_type *pqs1, qbstate12_type *pqs2, uint32_t nqb);
int32_t qbstate12_factor_modp(int32_t factor, uint32_t p);
extern const uint8_t qbstate12_lsbtab[64];
int32_t qbstate12_expand_b(qbstate12_type *pqs, uint8_t *a, uint32_t p, int32_t f);


#endif // #ifndef XSPECIAL12QS_H
