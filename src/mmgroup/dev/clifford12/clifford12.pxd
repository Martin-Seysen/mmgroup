
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



