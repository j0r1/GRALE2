from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool as cbool
from libcpp.string cimport string

cdef extern from "qpmatrix.h":
    cdef struct SparseMatrixInfo:
        vector[double] m_values
        vector[int] m_rows
        vector[int] m_cols

    ctypedef pair[SparseMatrixInfo, vector[double]] MatrixResults

    cdef cppclass MaskedPotentialValuesBase:
        int getNX() const
        int getNY() const
        int getNumberOfVariables() const
        cbool isValid() const
        string getInvalidReason() const

        double getInitialValue(int varIdx) const
        void getFullSolution(const vector[double] &sol, vector[double] &newPphiGrid) const

    cdef cppclass MaskedPotentialValues(MaskedPotentialValuesBase):
        MaskedPotentialValues(vector[double] &potentialValues, vector[cbool] &mask, int NX, double scaleUnit)

    cdef cppclass MaskedPotentialValuesOffsetGradient(MaskedPotentialValuesBase):
        MaskedPotentialValuesOffsetGradient(vector[double] &potentialValues, vector[int] &mask, int NX, double scaleUnit)

    MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValuesBase &mpv,
        const vector[pair[double, pair[int, int]]] &kernel)
    MatrixResults calculateLinearConstraintMatrices2(const MaskedPotentialValuesBase &mpv,
        const vector[pair[double, pair[int, int]]] &kernel, const vector[cbool] &relevantGridPositions,
        double limitingValue, cbool greaterThanLimitingValue)
    MatrixResults calculateLinearConstraintMatrices3(const MaskedPotentialValuesBase &mpv,
        const vector[pair[double, pair[int, int]]] &kernel,
        const vector[cbool] &relevantGridPositions,
        const vector[double] &limitingValues,
        cbool greaterThanLimitingValue)
    MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValuesBase &mpv,
		const vector[pair[double, pair[int, int]]] &kernel)


    
