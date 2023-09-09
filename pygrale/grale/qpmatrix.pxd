from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool as cbool

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

    cdef cppclass MaskedPotentialValues(MaskedPotentialValuesBase):
        MaskedPotentialValues(vector[double] &potentialValues, vector[cbool] &mask, int NX, double scaleUnit)

        const vector[double] &getPotentialValues() const
        const vector[cbool] &getMask() const

        double getInitialValue(int varIdx) const
        pair[int,int] getRowColumn(int varIdx) const

        double unadjustForUnit(double x) const

    cdef cppclass MaskedPotentialValuesOffsetGradient(MaskedPotentialValuesBase):
        MaskedPotentialValues(vector[double] &potentialValues, vector[int] &mask, int NX, double scaleUnit)

        const vector[double] &getPotentialValues() const
        const vector[int] &getMask() const

        double getInitialValue(int varIdx) const
        pair[int,int] getRowColumn(int varIdx) const

        double unadjustForUnit(double x) const

    MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValuesBase &mpv,
        const vector[pair[double, pair[int, int]]] &kernel)
    MatrixResults calculateLinearConstraintMatrices2(const MaskedPotentialValuesBase &mpv,
        const vector[pair[double, pair[int, int]]] &kernel, const vector[cbool] &relevantGridPositions,
        double limitingValue, cbool greaterThanLimitingValue)
    MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValuesBase &mpv,
		const vector[pair[double, pair[int, int]]] &kernel)

    
