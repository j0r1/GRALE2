#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <unordered_map>
#include <limits>
#include <cassert>

struct SparseMatrixInfo
{
    std::vector<double> m_values;
    std::vector<int> m_rows;
    std::vector<int> m_cols;
};

typedef std::pair<SparseMatrixInfo, std::vector<double>> MatrixResults; // Sparse matrix and column matrix

class MaskedPotentialValues
{
public:
	MaskedPotentialValues(std::vector<double> &potentialValues, std::vector<bool> &mask, int NX, double scaleUnit);
	int getNX() const { return m_NX; }
	int getNY() const { return m_NY; }
	const std::vector<double> &getPotentialValues() const { return m_potentialValues; }
	const std::vector<bool> &getMask() const { return m_mask; }

	std::pair<int, double> getVariableIndexOrValue(int i, int j) const;
	int getNumberOfVariables() const { return (int)m_idxMapFwd.size(); }
	double getInitialValue(int varIdx) const;
	std::pair<int,int> getRowColumn(int varIdx) const;

	double unadjustForUnit(double x) const { return x*m_scaleUnit; }
private:
	int m_NX, m_NY;
	std::vector<double> m_potentialValues;
	std::vector<bool> m_mask;
	std::vector<int> m_idxMapInv;
	std::vector<int> m_idxMapFwd;
	double m_scaleUnit;
};

inline std::pair<int, double> MaskedPotentialValues::getVariableIndexOrValue(int i, int j) const
{
	assert(i >= 0 && i < m_NY && j >= 0 && j < m_NX);

	int gridIdx = i*m_NX + j;
	if (m_mask[gridIdx]) // need to keep this value
		return { -1, m_potentialValues[gridIdx]/m_scaleUnit };

	assert(m_idxMapInv[gridIdx] == -1);
	return { m_idxMapInv[gridIdx], std::numeric_limits<double>::quiet_NaN() };
}

inline double MaskedPotentialValues::getInitialValue(int varIdx) const
{
	assert(varIdx >= 0 && varIdx < getNumberOfVariables());
	int gridIdx = m_idxMapFwd[varIdx];
	assert(gridIdx >= 0 && gridIdx < m_potentialValues.size());
	return m_potentialValues[gridIdx]/m_scaleUnit;
}

inline std::pair<int,int> MaskedPotentialValues::getRowColumn(int varIdx) const
{
	assert(varIdx >= 0 && varIdx < getNumberOfVariables());
	int gridIdx = m_idxMapFwd[varIdx];
	assert(gridIdx >= 0 && gridIdx < m_potentialValues.size());
	return { gridIdx/m_NX, gridIdx%m_NX };
}

MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValues &mpv,
		const std::vector<std::pair<double, std::pair<int, int>>> &kernel
		);

MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValues &mpv,
		const std::vector<std::pair<double,std::vector<std::pair<double, std::pair<int, int>>>>> &kernelList
		);

