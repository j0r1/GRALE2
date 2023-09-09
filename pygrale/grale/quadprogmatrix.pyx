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

def _getReturnType(returnType):
    if returnType == "csc":
        return 0
    if returnType == "csr":
        return 1
    if returnType == "raw":
        return 2
    raise MaskedPotentialValuesException("Return type must be 'csc', 'csr' or 'raw'")

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

cdef class MaskedPotentialValues:
    cdef unique_ptr[qpmatrix.MaskedPotentialValues] m_maskedValues

    def __init__(self, np.ndarray[double, ndim=2] potentialValues, np.ndarray[cbool, ndim=2] mask, double phiScale):
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

        self.m_maskedValues = make_unique[qpmatrix.MaskedPotentialValues](pVal, cMask, cols, phiScale)

    def getNumberOfVariables(self):
        return deref(self.m_maskedValues).getNumberOfVariables()

    def getLinearConstraintMatrices(self, kernel, returnType = "csc", np.ndarray[cbool, ndim=2] relevantGridPositions=None, double limitingValue=0):
        cdef vector[pair[double, pair[int, int]]] cKernel
        cdef vector[cbool] cRelGridPos
        cdef pair[int,int] diff
        cdef qpmatrix.MatrixResults results
        cdef np.ndarray[double, ndim=1] values, b
        cdef np.ndarray[long, ndim=1] rows, cols
        cdef int i, j, rt, N, M
        cdef double factor

        rt = _getReturnType(returnType)

        for part in kernel:
            factor = part["factor"]
            diff = pair[int, int](part["di"], part["dj"])
            cKernel.push_back(pair[double, pair[int,int]](factor, diff))

        if relevantGridPositions is None:
            if limitingValue != 0:
                raise MaskedPotentialValuesException("'limitingValue' should be zero in this case")

            results = qpmatrix.calculateLinearConstraintMatrices(deref(self.m_maskedValues), cKernel)
        else:
            if deref(self.m_maskedValues).getNX() != relevantGridPositions.shape[1] or deref(self.m_maskedValues).getNY() != relevantGridPositions.shape[0]:
                raise MaskedPotentialValuesException("Shape of 'relevantGridPositions' must match that of the underlying grid")

            for i in range(relevantGridPositions.shape[0]):
                for j in range(relevantGridPositions.shape[1]):
                    cRelGridPos.push_back(relevantGridPositions[i,j])

            # TODO: for now, the smaller than case is forced
            # Note that the limitingValue is passed as is, the relevant scale factor is
            # calculated in the createEquivalentPotentialGridLens function
            results = qpmatrix.calculateLinearConstraintMatrices2(deref(self.m_maskedValues), cKernel, cRelGridPos, limitingValue, False)

        N = results.second.size()
        M = self.getNumberOfVariables()
        return _returnResults(&results, rt, N, M)

    def getQuadraticMinimizationMatrices(self, kernel, returnType = "csc"):
        cdef vector[pair[double,pair[int, int]]] cKernel
        cdef double weight, factor
        cdef pair[int,int] diff
        cdef qpmatrix.MatrixResults results
        cdef int rt, N
        
        rt = _getReturnType(returnType)

        for part in kernel:
            factor = part["factor"]
            diff = pair[int, int](part["di"], part["dj"])
            cKernel.push_back(pair[double, pair[int,int]](factor, diff))
            
        results = qpmatrix.calculateQuadraticMimimizationMatrices(deref(self.m_maskedValues), cKernel)
        N = self.getNumberOfVariables()
        return _returnResults(&results, rt, N, N)

    def getInitialValues(self):
        cdef int N,i
        cdef np.ndarray[double, ndim=1] values

        N = self.getNumberOfVariables()
        values = np.empty([ N ], dtype=np.double)
        for i in range(N):
            values[i] = deref(self.m_maskedValues).getInitialValue(i)

        return values

    def getFullSolution(self, np.ndarray[double, ndim=1] solValues):
        cdef np.ndarray[double, ndim=2] newPhi
        cdef int N, NX, NY, i, j, idx
        cdef const vector[double] *oldPhi = &(deref(self.m_maskedValues).getPotentialValues())
        cdef pair[int,int] rowCol
        cdef double val

        N = self.getNumberOfVariables()
        if N != solValues.shape[0]:
            raise MaskedPotentialValuesException("Expecting an input array of length {}".format(N))

        NX = deref(self.m_maskedValues).getNX()
        NY = deref(self.m_maskedValues).getNY()

        # First create a copy of the old values
        newPhi = np.empty([NY,NX], dtype=np.double)
        idx = 0
        for i in range(NY):
            for j in range(NX):
                newPhi[i,j] = (deref(oldPhi))[idx]
                idx += 1

        # Fill in solution
        for idx in range(N):
            rowCol = deref(self.m_maskedValues).getRowColumn(idx)

            val = solValues[idx]
            val = deref(self.m_maskedValues).unadjustForUnit(val)
            newPhi[rowCol.first, rowCol.second] = val

        return newPhi

cdef class MaskedPotentialValuesOffsetGradient:
    cdef unique_ptr[qpmatrix.MaskedPotentialValuesOffsetGradient] m_maskedValues

    def __init__(self, np.ndarray[double, ndim=2] potentialValues, np.ndarray[int, ndim=2] mask, double phiScale):
        cdef vector[double] pVal
        cdef vector[int] cMask
        cdef int i, j, rows, cols, tmp
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
                tmp = mask[i,j]
                if tmp < 0 or tmp > 2:
                    raise MaskedPotentialValuesException("Mask value must be either 0 (optimize potential value), 1 (fixed potential value) or 2 (fixed potential value, but with offset and gradient)")

                cMask.push_back(tmp)

        self.m_maskedValues = make_unique[qpmatrix.MaskedPotentialValuesOffsetGradient](pVal, cMask, cols, phiScale)

    def getNumberOfVariables(self):
        return deref(self.m_maskedValues).getNumberOfVariables()

    def getLinearConstraintMatrices(self, kernel, returnType = "csc"):
        cdef vector[pair[double, pair[int, int]]] cKernel
        cdef pair[int,int] diff
        cdef qpmatrix.MatrixResults results
        cdef np.ndarray[double, ndim=1] values, b
        cdef np.ndarray[long, ndim=1] rows, cols
        cdef int i, rt, N, M
        cdef double factor

        rt = _getReturnType(returnType)

        for part in kernel:
            factor = part["factor"]
            diff = pair[int, int](part["di"], part["dj"])
            cKernel.push_back(pair[double, pair[int,int]](factor, diff))

        results = qpmatrix.calculateLinearConstraintMatrices(deref(self.m_maskedValues), cKernel)
        N = results.second.size()
        M = self.getNumberOfVariables()
        return _returnResults(&results, rt, N, M)

    def getQuadraticMinimizationMatrices(self, kernelList, returnType = "csc"):
        raise Exception("TODO")
        # qpmatrix.calculateQuadraticMimimizationMatrices now only handles a single
        # kernel, this function should either call that multiple times, or rewrite
        # the rest of the code slightly - since all this is experimental, just
        # postponing that for now

#        cdef vector[pair[double, vector[pair[double,pair[int, int]]]]] cKernelList
#        cdef vector[pair[double,pair[int, int]]] cKernel
#        cdef double weight, factor
#        cdef pair[int,int] diff
#        cdef qpmatrix.MatrixResults results
#        cdef int rt, N
#        
#        rt = _getReturnType(returnType)
#
#        for d in kernelList:
#            weight = d["weight"]
#            kernel = d["kernel"]
#            cKernel.clear()
#            
#            for part in kernel:
#                factor = part["factor"]
#                diff = pair[int, int](part["di"], part["dj"])
#                cKernel.push_back(pair[double, pair[int,int]](factor, diff))
#            
#            cKernelList.push_back(pair[double, vector[pair[double,pair[int, int]]]](weight, cKernel))
#
#        results = qpmatrix.calculateQuadraticMimimizationMatrices(deref(self.m_maskedValues), cKernelList)
#        N = self.getNumberOfVariables()
#        return _returnResults(&results, rt, N, N)

    # NOTE: for now the underlying code just returns zero, should we just return a zero ndarray?
    def getInitialValues(self):
        cdef int N,i
        cdef np.ndarray[double, ndim=1] values

        N = self.getNumberOfVariables()
        values = np.empty([ N ], dtype=np.double)
        for i in range(N):
            values[i] = deref(self.m_maskedValues).getInitialValue(i)

        return values

    def getFullSolution(self, np.ndarray[double, ndim=1] solValues):
        cdef np.ndarray[double, ndim=2] newPhi
        cdef int N, NX, NY, i, j, idx, m
        cdef pair[int,int] rowCol
        cdef double val
        cdef qpmatrix.MaskedPotentialValuesOffsetGradient *pMV = self.m_maskedValues.get()
        cdef const vector[double] *oldPhi = &(pMV.getPotentialValues())
        cdef const vector[int] *mask = &(pMV.getMask())

        N = self.getNumberOfVariables()
        if N != solValues.shape[0]:
            raise MaskedPotentialValuesException("Expecting an input array of length {}".format(N))

        NX = pMV.getNX()
        NY = pMV.getNY()

        # First create a copy of the old values, or with offset and gradient
        newPhi = np.empty([NY,NX], dtype=np.double)
        newPhi[:,:] = float("NaN") # Just to check that we're willing in everything

        solOffset = solValues[0]
        solGradX = solValues[1]
        solGradY = solValues[2]

        idx = 0
        for i in range(NY):
            for j in range(NX):
                m = deref(mask)[idx]
                if m == 0:
                    pass # will fill in later
                elif m == 1: # Just the initial value
                    newPhi[i,j] = (deref(oldPhi))[idx]
                elif m == 2: # The initial value, plus offset and gradient
                    newPhi[i,j] = (deref(oldPhi))[idx] + pMV.unadjustForUnit(solOffset + j*solGradX + i*solGradY)
                else:
                    raise MaskedPotentialValuesException("Internal error: unknown mask value {} detected, should have been detected already".format(m))
                idx += 1

        # Fill in solution
        for idx in range(3, N): # the first three are for offset and gradient
            rowCol = pMV.getRowColumn(idx)

            val = solValues[idx]
            val = pMV.unadjustForUnit(val)
            newPhi[rowCol.first, rowCol.second] = val

        assert np.sum(np.isnan(newPhi)) == 0, "NaN detected in final grid, something is not filled in?"

        return newPhi



