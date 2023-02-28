#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <unordered_map>
#include <limits>
#include <cassert>
#include <cstdlib>

class MaskedPotentialValuesBase
{
public:
	MaskedPotentialValuesBase(int NX, int NY) : m_NX(NX), m_NY(NY), m_numVariables(-1)
	{
		assert(m_NX > 0 && m_NY > 0);
	}

	virtual ~MaskedPotentialValuesBase() { }

	int getNX() const { return m_NX; }
	int getNY() const { return m_NY; }
	int getNumberOfVariables() const { return m_numVariables; }

	// TODO: rewrite this so that a buffer for the return value is passed, for efficiency
	//       Can't just return a reference as it may be used in a double loop
	virtual std::vector<std::tuple<int, double, double>> getVariableIndexOrValue(int i, int j) const = 0;
protected:
	void setNumberOfVariables(int n) { m_numVariables = n; }
	const int m_NX, m_NY;
private:
	int m_numVariables;
};

struct SparseMatrixInfo
{
    std::vector<double> m_values;
    std::vector<int> m_rows;
    std::vector<int> m_cols;
};

typedef std::pair<SparseMatrixInfo, std::vector<double>> MatrixResults; // Sparse matrix and column matrix

class MaskedPotentialValues : public MaskedPotentialValuesBase
{
public:
	MaskedPotentialValues(std::vector<double> &potentialValues, std::vector<bool> &mask, int NX, double scaleUnit);
	const std::vector<double> &getPotentialValues() const { return m_potentialValues; }
	const std::vector<bool> &getMask() const { return m_mask; }

	std::vector<std::tuple<int, double, double>> getVariableIndexOrValue(int i, int j) const override;
	double getInitialValue(int varIdx) const;
	std::pair<int,int> getRowColumn(int varIdx) const;

	double unadjustForUnit(double x) const { return x*m_scaleUnit; }
private:
	std::vector<double> m_potentialValues;
	std::vector<bool> m_mask;
	std::vector<int> m_idxMapInv;
	std::vector<int> m_idxMapFwd;
	double m_scaleUnit;
};

inline std::vector<std::tuple<int, double, double>> MaskedPotentialValues::getVariableIndexOrValue(int i, int j) const
{
	assert(i >= 0 && i < m_NY && j >= 0 && j < m_NX);

	int gridIdx = i*m_NX + j;
	if (m_mask[gridIdx]) // need to keep this value
		return { { -1, m_potentialValues[gridIdx]/m_scaleUnit, 1.0 } };

	assert(m_idxMapInv[gridIdx] >= 0);
	return { { m_idxMapInv[gridIdx], std::numeric_limits<double>::quiet_NaN(), 1.0 } };
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

MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValuesBase &mpv,
		const std::vector<std::pair<double, std::pair<int, int>>> &kernel
		);

MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValuesBase &mpv,
		const std::vector<std::pair<double,std::vector<std::pair<double, std::pair<int, int>>>>> &kernelList
		);

class MaskedPotentialValuesOffsetGradient : public MaskedPotentialValuesBase
{
public:
	MaskedPotentialValuesOffsetGradient(std::vector<double> &potentialValues, std::vector<int> &mask, int NX, double scaleUnit);
	const std::vector<double> &getPotentialValues() const { return m_potentialValues; }
	const std::vector<int> &getMask() const { return m_mask; }

	std::vector<std::tuple<int, double, double>> getVariableIndexOrValue(int i, int j) const override;

	int getNumberOfVariables() const { return (int)m_idxMapFwd.size(); }
	double getInitialValue(int varIdx) const;
	std::pair<int,int> getRowColumn(int varIdx) const;

	double unadjustForUnit(double x) const { return x*m_scaleUnit; }
private:
	std::vector<double> m_potentialValues;
	std::vector<int> m_mask;
	std::vector<int> m_idxMapInv;
	std::vector<int> m_idxMapFwd;
	double m_scaleUnit;
};

inline std::vector<std::tuple<int, double, double>> MaskedPotentialValuesOffsetGradient::getVariableIndexOrValue(int i, int j) const
{
	assert(i >= 0 && i < m_NY && j >= 0 && j < m_NX);

	int gridIdx = i*m_NX + j;

	assert(gridIdx >= 0);
	assert(gridIdx < m_potentialValues.size());
	assert(gridIdx < m_mask.size());

	int maskValue = m_mask[gridIdx];
	assert(maskValue >= 0 && maskValue <= 2);

	if (maskValue == 0) // this is a position for which the value should be optimized
	{
		return { { m_idxMapInv[gridIdx], std::numeric_limits<double>::quiet_NaN(), 1.0 } };
	}
	else if (maskValue == 1) // this is a position for which the potential value needs to be kept
	{
		assert(m_idxMapInv[gridIdx] == -1);
		return { { -1, m_potentialValues[gridIdx]/m_scaleUnit, 1.0 } };
	}
	else if (maskValue == 2) // the value is kept, but an offset and gradient are added
	{
		assert(m_idxMapInv[gridIdx] == -1);

		int idxOffset = 0;
		int idxGradX = 1;
		int idxGradY = 2;
		return {
			{ -1, m_potentialValues[gridIdx]/m_scaleUnit, 1.0 },
			{ idxOffset, std::numeric_limits<double>::quiet_NaN(), 1.0 },
			{ idxGradX, std::numeric_limits<double>::quiet_NaN(), j },
			{ idxGradY, std::numeric_limits<double>::quiet_NaN(), i }
		};
	}

	std::cerr << "Internal error! shouldn't reach this" << std::endl;
	exit(-1);
	return {};
}

inline double MaskedPotentialValuesOffsetGradient::getInitialValue(int varIdx) const
{
	// TODO: how much sense does this make? We don't know from where to start, so we can't really find
	//       something that's a feasible solution if constraints are used. Don't think we need constraints,
	//       and then the initial value doesn't really matter
	return 0;
}

inline std::pair<int,int> MaskedPotentialValuesOffsetGradient::getRowColumn(int varIdx) const
{
	assert(varIdx >= 0 && varIdx < getNumberOfVariables());
	int gridIdx = m_idxMapFwd[varIdx];
	assert(gridIdx >= 0 && gridIdx < m_potentialValues.size());
	return { gridIdx/m_NX, gridIdx%m_NX };
}
