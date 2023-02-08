from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from libcpp.memory cimport unique_ptr, make_unique
from cython.operator cimport dereference as deref
from libcpp.pair cimport pair
import cython
import numpy as np
from scipy import sparse
cimport numpy as np

cimport grale.qpmatrix as qpmatrix

include "stringwrappers.pyx"

ctypedef np.ndarray ndarray

class MaskedPotentialValuesException(Exception):
    pass

cdef class MaskedPotentialValues:
    cdef unique_ptr[qpmatrix.MaskedPotentialValues] m_maskedValues

    def __init__(self, np.ndarray[double, ndim=2] potentialValues, np.ndarray[cbool, ndim=2] mask):
        cdef vector[double] pVal
        cdef vector[cbool] cMask
        cdef int i, j, rows, cols
        cdef double unit

        assert potentialValues.shape[0] == mask.shape[0] and potentialValues.shape[1] == mask.shape[1], "Shapes of potentialValues and mask must match"
        rows = potentialValues.shape[0]
        cols = potentialValues.shape[1]
        assert rows > 0 and cols > 0, "Some values must be present"

        pVal.reserve(rows*cols)
        cMask.reserve(rows*cols)

        for i in range(rows):
            for j in range(cols):
                pVal.push_back(potentialValues[i,j])
                cMask.push_back(mask[i,j])

        from grale.constants import ANGLE_ARCSEC
        unit = ANGLE_ARCSEC**2
        self.m_maskedValues = make_unique[qpmatrix.MaskedPotentialValues](pVal, cMask, cols, unit)

    def getNumberOfVariables(self):
        return deref(self.m_maskedValues).getNumberOfVariables()

    @staticmethod
    cdef _returnResults(qpmatrix.MatrixResults *results, int returnType, int N, int M):
        cdef np.ndarray[double, ndim=1] values, b

        b = np.empty([ results.second.size() ], dtype=np.double)
        for i in range(results.second.size()):
            b[i] = results.second[i]
        
        assert results.first.m_values.size() == results.first.m_rows.size() and results.first.m_rows.size() == results.first.m_cols.size(), "Internal error: returned arrays don't have same length"
        values = np.empty([ results.first.m_values.size() ], dtype=np.double)
        rows = np.empty([ results.first.m_rows.size() ], dtype=int)
        cols = np.empty([ results.first.m_cols.size() ], dtype=int)
        for i in range(results.first.m_values.size()):
            values[i] = results.first.m_values[i]
            rows[i] = results.first.m_rows[i]
            cols[i] = results.first.m_cols[i]

        if returnType == 0:
            return (sparse.csc_matrix((values, (rows,cols)), shape=(N,M), dtype=np.double), b)
        if returnType == 1:
            return (sparse.csr_matrix((values, (rows,cols)), shape=(N,M), dtype=np.double), b)
        return ( { "values": values, "rows": rows, "columns": cols }, b )

    @staticmethod
    def _getReturnType(returnType):
        if returnType == "csc":
            return 0
        if returnType == "csr":
            return 1
        if returnType == "raw":
            return 2
        raise MaskedPotentialValuesException("Return type must be 'csc', 'csr' or 'raw'")

    def getLinearConstraintMatrices(self, kernel, returnType = "csc"):
        cdef vector[pair[double, pair[int, int]]] cKernel
        cdef pair[int,int] diff
        cdef qpmatrix.MatrixResults results
        cdef np.ndarray[double, ndim=1] values, b
        cdef np.ndarray[long, ndim=1] rows, cols
        cdef int i, rt, N, M
        cdef double factor

        rt = MaskedPotentialValues._getReturnType(returnType)

        for part in kernel:
            factor = part["factor"]
            diff = pair[int, int](part["di"], part["dj"])
            cKernel.push_back(pair[double, pair[int,int]](factor, diff))

        results = qpmatrix.calculateLinearConstraintMatrices(deref(self.m_maskedValues), cKernel)
        N = results.second.size()
        M = self.getNumberOfVariables()
        return MaskedPotentialValues._returnResults(&results, rt, N, M)

    def getQuadraticMinimizationMatrices(self, kernelList, returnType = "csc"):
        cdef vector[pair[double, vector[pair[double,pair[int, int]]]]] cKernelList
        cdef vector[pair[double,pair[int, int]]] cKernel
        cdef double weight, factor
        cdef pair[int,int] diff
        cdef qpmatrix.MatrixResults results
        cdef int rt, N
        
        rt = MaskedPotentialValues._getReturnType(returnType)

        for d in kernelList:
            weight = d["weight"]
            kernel = d["kernel"]
            cKernel.clear()
            
            for part in kernel:
                factor = part["factor"]
                diff = pair[int, int](part["di"], part["dj"])
                cKernel.push_back(pair[double, pair[int,int]](factor, diff))
            
            cKernelList.push_back(pair[double, vector[pair[double,pair[int, int]]]](weight, cKernel))

        results = qpmatrix.calculateQuadraticMimimizationMatrices(deref(self.m_maskedValues), cKernelList)
        N = self.getNumberOfVariables()
        return MaskedPotentialValues._returnResults(&results, rt, N, N)

