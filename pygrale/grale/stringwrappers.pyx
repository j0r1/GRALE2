from cpython.version cimport PY_MAJOR_VERSION
from cpython.array cimport array

cdef B(s):
    if PY_MAJOR_VERSION < 3:
        return s
    return bytes(s, 'UTF-8')

cdef S(b):
    if PY_MAJOR_VERSION < 3:
        return b
    return b.decode(encoding='UTF-8')

cdef fromstring(tmp, b):
    if PY_MAJOR_VERSION < 3:
        tmp.fromstring(b)
    else:
        tmp.frombytes(b)

cdef chararrayfrombytes(b):
    tmp = array('b')
    fromstring(tmp, b)
    return tmp
