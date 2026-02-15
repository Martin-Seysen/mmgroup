
from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t
from libc.stdint cimport int64_t, int32_t, int8_t


cdef extern from "clifford12.h":
    enum: QSTATE12_MAXCOLS
    enum: QSTATE12_MAXROWS
    enum: QSTATE12_UNDEF_ROW
    ctypedef struct qstate12_type:
        uint32_t maxrows
        uint32_t nrows
        uint32_t ncols
        int32_t  factor
        uint32_t shape1
        uint32_t reduced
        uint64_t *data 

    ctypedef struct qstate12_support_type:
        uint32_t size         # Total size of state; a power of 2
        uint32_t weight       # No of nonzero entries of of matrix; 0 or a power of 2
        int32_t factor        # Scalar factor, encoded as in the structure ``qstate12``
        uint32_t factor_new   # 1 if factor has changed, 0 otherwise
        uint32_t batchlength  # Size of current batch; 0 or a power of 2
        uint32_t n_batches    # Total number of batches
        uint32_t indices[64]  # Indices contained in the current batch
        uint8_t signs[64]     # Signs of entries of the current batch
        uint32_t c_weight


