from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool as cbool

cdef extern from "qpmatrix.h":
    cdef struct SparseMatrixInfo:
        vector[double] m_values
        vector[int] m_rows
        vector[int] m_cols

    ctypedef pair[SparseMatrixInfo, vector[double]] MatrixResults

    cdef cppclass MaskedPotentialValues:
        MaskedPotentialValues(vector[double] &&potentialValues, vector[cbool] &&mask, int NX, double scaleUnit)
        int getNX() const
        int getNY() const
        const vector[double] &getPotentialValues() const
        const vector[cbool] &getMask() const

        pair[int,double] getVariableIndexOrValue(int i, int j) const
        int getNumberOfVariables() const
        double getInitialValue(int varIdx) const
        pair[int,int] getRowColumn(int varIdx) const

        double unadjustForUnit(double x) const

    MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValues &mpv,
        const vector[pair[double, pair[int, int]]] &kernel)
    MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValues &mpv,
        const vector[pair[double,vector[pair[double, pair[int, int]]]]] &kernelList)

    
