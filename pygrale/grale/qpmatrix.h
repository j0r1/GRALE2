#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <unordered_map>
#include <limits>
#include <cassert>
#include <cstdlib>

struct SparseMatrixInfo
{
    std::vector<double> m_values;
    std::vector<int> m_rows;
    std::vector<int> m_cols;
};

typedef std::pair<SparseMatrixInfo, std::vector<double>> MatrixResults; // Sparse matrix and column matrix

class MaskedPotentialValuesBase
{
public:
	MaskedPotentialValuesBase(int NX, int NY) : m_NX(NX), m_NY(NY), m_numVariables(-1)
	{
		assert(m_NX > 0 && m_NY > 0);
	}

	virtual ~MaskedPotentialValuesBase() { }
	bool isValid() const { return m_invalidReason.empty(); }
	const std::string getInvalidReason() const { return m_invalidReason; }

	int getNX() const { return m_NX; }
	int getNY() const { return m_NY; }
	int getNumberOfVariables() const { return m_numVariables; }

	// TODO: rewrite this so that a buffer for the return value is passed, for efficiency
	//       Can't just return a reference as it may be used in a double loop
	virtual std::vector<std::tuple<int, double, double>> getVariableIndexOrValue(int i, int j) const = 0;

	virtual double getInitialValue(int varIdx) const = 0;
	virtual void getFullSolution(const std::vector<double> &sol, std::vector<double> &newPhiGrid) const = 0;
protected:
	void setInvalid(const std::string &reason) { m_invalidReason = reason; }
	void setNumberOfVariables(int n) { m_numVariables = n; }
	const int m_NX, m_NY;
private:
	int m_numVariables;
	std::string m_invalidReason;
};


class MaskedPotentialValues : public MaskedPotentialValuesBase
{
public:
	MaskedPotentialValues(std::vector<double> &potentialValues, std::vector<bool> &mask, int NX, double scaleUnit);

	std::vector<std::tuple<int, double, double>> getVariableIndexOrValue(int i, int j) const override;
	double getInitialValue(int varIdx) const override;

	void getFullSolution(const std::vector<double> &sol, std::vector<double> &newPhiGrid) const override;

	double unadjustForUnit(double x) const { return x*m_scaleUnit; }
	std::pair<int,int> getRowColumn(int varIdx) const;
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

inline void MaskedPotentialValues::getFullSolution(const std::vector<double> &sol, std::vector<double> &newPhi) const
{
	assert(sol.size() == getNumberOfVariables());
	assert(newPhi.size() == m_potentialValues.size());

	// First copy old values
	newPhi = m_potentialValues;

	// Fill in solution
	for (size_t idx = 0 ; idx < sol.size() ; idx++)
	{
		double val = sol[idx];
		val = unadjustForUnit(val);

		auto rowCol = getRowColumn(idx);
		assert(rowCol.first < getNY() && rowCol.second < getNX());

		size_t gridIdx = rowCol.first * getNX() + rowCol.second;
		assert(gridIdx < newPhi.size());
		newPhi[gridIdx] = val;
	}
}

MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValuesBase &mpv,
		const std::vector<std::pair<double, std::pair<int, int>>> &kernel
		);

MatrixResults calculateLinearConstraintMatrices2(const MaskedPotentialValuesBase &mpv,
		const std::vector<std::pair<double, std::pair<int, int>>> &kernel,
		const std::vector<bool> &relevantGridPositions,
		double limitingValue,
		bool greaterThanLimitingValue
		);

MatrixResults calculateLinearConstraintMatrices3(const MaskedPotentialValuesBase &mpv,
		const std::vector<std::pair<double, std::pair<int, int>>> &kernel,
		const std::vector<bool> &relevantGridPositions,
		const std::vector<double> &limitingValues,
		bool greaterThanLimitingValue
		);

MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValuesBase &mpv,
		const std::vector<std::pair<double, std::pair<int, int>>> &kernel
		);

class MaskedPotentialValuesOffsetGradient : public MaskedPotentialValuesBase
{
public:
	MaskedPotentialValuesOffsetGradient(std::vector<double> &potentialValues, std::vector<int> &mask, int NX, double scaleUnit);

	std::vector<std::tuple<int, double, double>> getVariableIndexOrValue(int i, int j) const override;
	double getInitialValue(int varIdx) const override;
	void getFullSolution(const std::vector<double> &sol, std::vector<double> &newPhiGrid) const override;

	double unadjustForUnit(double x) const { return x*m_scaleUnit; }
	std::pair<int,int> getRowColumn(int varIdx) const;
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
			// Below, by just using j instead of theta_x,j or j instead of theta_y,i, we're
			// actually absorbing a scale factor into grad_x and grad_y
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
	assert(varIdx >= 3 && varIdx < getNumberOfVariables()); // First three are offset, grad_x and grad_y, don't have a grid position
	int gridIdx = m_idxMapFwd[varIdx];
	assert(gridIdx >= 0 && gridIdx < m_potentialValues.size());
	return { gridIdx/m_NX, gridIdx%m_NX };
}

inline void MaskedPotentialValuesOffsetGradient::getFullSolution(const std::vector<double> &sol, std::vector<double> &newPhi) const
{
	assert(sol.size() == getNumberOfVariables());
	assert(newPhi.size() == m_potentialValues.size());

	assert(sol.size() > 3);
	double solOffset = sol[0];
	double solGradX = sol[1];
	double solGradY = sol[2];

	for (size_t idx = 0 ; idx < newPhi.size() ; idx++)
	{
		int m = m_mask[idx];
		if (m == 0)
		{
			// Nothing to do here, will fill in later
		}
		else if (m == 1)
		{
			// Keep the initial value
			newPhi[idx] = m_potentialValues[idx];
		}
		else if (m == 2)
		{
			// Initial value, corrected by offset and gradient
			size_t i = idx/getNX();
			size_t j = idx%getNX();

			double correction = solOffset + j*solGradX + i*solGradY;
			newPhi[idx] = m_potentialValues[idx] + unadjustForUnit(correction);
		}
		else
		{
			std::cerr << "Internal error: invalid mask value, but this should have been checked before!" << std::endl;
			exit(-1);
		}
	}

	// Fill in the rest of the solution
	for (size_t idx = 3 ; idx < sol.size() ; idx++) // First three for offset, grad_x and grad_y
	{
		double val = sol[idx];
		val = unadjustForUnit(val);

		auto rowCol = getRowColumn(idx);
		assert(rowCol.first < getNY() && rowCol.second < getNX());

		size_t gridIdx = rowCol.first * getNX() + rowCol.second;
		assert(gridIdx < newPhi.size());
		newPhi[gridIdx] = val;
	}
}

