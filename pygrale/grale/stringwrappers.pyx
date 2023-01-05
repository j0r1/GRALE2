from cpython.array cimport array

cdef B(s):
    return bytes(s, 'UTF-8')

cdef S(b):
    return b.decode(encoding='UTF-8')

cdef chararrayfrombytes(b):
    tmp = array('b')
    tmp.frombytes(b)
    return tmp
