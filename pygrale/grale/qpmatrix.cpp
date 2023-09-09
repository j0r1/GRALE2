#include "qpmatrix_templates.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

class GridPos
{
public:
	GridPos(int i, int j, int N, int M) : m_i (i), m_j(j), m_N(N), m_M(M) { }
	GridPos(GridPos p, pair<int,int> diff)
	{
		m_i = p.m_i + diff.first;
		m_j = p.m_j + diff.second;
		m_N = p.m_N;
		m_M = p.m_M;
	}

	GridPos &operator++()
	{
		assert(m_i < m_N);
		assert(m_j < m_M);
		m_j++;
		if (m_j == m_M)
		{
			m_j = 0;
			m_i++;
		}
		return *this;
	}
	bool operator==(const GridPos &p) const
	{
		assert(m_N == p.m_N);
		assert(m_M == p.m_M);
		return m_i == p.m_i && m_j == p.m_j;
	}

	bool operator!=(const GridPos &p) const
	{
		assert(m_N == p.m_N);
		assert(m_M == p.m_M);
		return m_i != p.m_i || m_j != p.m_j;
	}

	GridPos operator*() const { return { m_i, m_j, m_N, m_M }; }

	string toString() const { return "(" + to_string(m_i) + "," + to_string(m_j) + ")"; }

	int m_i, m_j, m_N, m_M;
};

class Grid
{
public:
	Grid(int N, int M) : m_N(N), m_M(M) { }

	GridPos begin() const { return { 0, 0, m_N, m_M }; }
	GridPos end() const { return { m_N, 0, m_N, m_M }; }

	bool isPositionValid(GridPos x) const
	{
		assert(x.m_N == m_N && x.m_M == m_M);
		return x.m_i >= 0 && x.m_j >= 0 && x.m_i < m_N && x.m_j < m_M;
	}
private:
	int m_N, m_M;
};

class MaskedGridPos
{
public:
	MaskedGridPos(int i, int j, int N, int M, const shared_ptr<vector<bool>> &relGridPos) : m_i (i), m_j(j), m_N(N), m_M(M), m_mask(relGridPos) { }

	MaskedGridPos(const MaskedGridPos &p, pair<int,int> diff)
	{
		m_i = p.m_i + diff.first;
		m_j = p.m_j + diff.second;
		m_N = p.m_N;
		m_M = p.m_M;
		m_mask = p.m_mask;
	}

	MaskedGridPos &operator++()
	{
		assert(m_i < m_N);
		assert(m_j < m_M);

		const auto &mask = *m_mask;
		int idx = m_j + m_i*m_M;
		do
		{
			m_j++;
			if (m_j == m_M)
			{
				m_j = 0;
				m_i++;
			}
			idx++;

			if (idx >= mask.size()) // We should be done
				break;

			if (mask[idx]) // Ok, it's a valid point
				break;
		} while (true);

		return *this;
	}

	bool operator==(const MaskedGridPos &p) const
	{
		assert(m_N == p.m_N);
		assert(m_M == p.m_M);
		return m_i == p.m_i && m_j == p.m_j;
	}

	bool operator!=(const MaskedGridPos &p) const
	{
		assert(m_N == p.m_N);
		assert(m_M == p.m_M);
		return m_i != p.m_i || m_j != p.m_j;
	}

	MaskedGridPos operator*() const { return { m_i, m_j, m_N, m_M, m_mask }; }

	string toString() const { return "(" + to_string(m_i) + "," + to_string(m_j) + ")"; }

	int m_i, m_j, m_N, m_M;
	shared_ptr<vector<bool>> m_mask;
};

class MaskedGrid
{
public:
	MaskedGrid(int N, int M, const vector<bool> &relGridPos) : m_N(N), m_M(M)
	{
		assert(N*M == relGridPos.size());
		m_relGridPos = make_shared<vector<bool>>();
		
		std::copy(relGridPos.begin(), relGridPos.end(), std::back_inserter(*m_relGridPos));
	}

	MaskedGridPos begin() const
	{ 
		int i = 0, j = 0, idx = 0;
		const auto &mask = *m_relGridPos;
		do
		{
			if (mask[idx])
				break;

			j++;
			if (j == m_M)
			{
				j = 0;
				i++;
			}
			idx++;
			if (idx >= mask.size())
				break;
		} while(true);

		return { i, j, m_N, m_M, m_relGridPos };
	}
	MaskedGridPos end() const { return { m_N, 0, m_N, m_M, m_relGridPos }; }

	shared_ptr<vector<bool>> m_relGridPos;
private:
	int m_N, m_M;
};

MaskedPotentialValues::MaskedPotentialValues(vector<double> &potentialValues, vector<bool> &mask, int NX, double scaleUnit)
	: MaskedPotentialValuesBase(NX, potentialValues.size()/NX)
{
	assert(potentialValues.size() == mask.size());
	assert(potentialValues.size() % NX == 0);
	
	m_scaleUnit = scaleUnit;

	m_potentialValues = move(potentialValues);
	m_mask = move(mask);

	m_idxMapInv.resize(m_potentialValues.size(), -1);
	for (size_t i = 0 ; i < m_potentialValues.size() ; i++)
	{
		if (!m_mask[i])
		{
			m_idxMapInv[i] = m_idxMapFwd.size();
			m_idxMapFwd.push_back(i);
		}
	}

	setNumberOfVariables((int)m_idxMapFwd.size());
}

MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValuesBase &mpv,
		const vector<pair<double, pair<int, int>>> &kernel
		)
{
	const int NX = mpv.getNX();
	const int NY = mpv.getNY();
	Grid gridPositions(NY, NX);

	return calculateLinearMatrix_Functors(gridPositions,
			[NX,NY](auto pos) { return pos.m_i >= 0 && pos.m_i < NY && pos.m_j >= 0 && pos.m_j < NX; },
			kernel,
			[](auto pos, auto diff) { return GridPos(pos, diff); },
			[&mpv](auto pos) { return mpv.getVariableIndexOrValue(pos.m_i, pos.m_j); },
			0.0,
			true
			);
}

MatrixResults calculateLinearConstraintMatrices2(const MaskedPotentialValuesBase &mpv,
		const vector<pair<double, pair<int, int>>> &kernel,
		const std::vector<bool> &relevantGridPositions,
		double limitingValue,
		bool greaterThanLimitingValue
		)
{
	const int NX = mpv.getNX();
	const int NY = mpv.getNY();
	assert(relevantGridPositions.size() == NX*NY);

	/*
	// TODO: make this more efficient !!
	Grid gridPositions(NY, NX);
	size_t idx = 0;
	vector<GridPos> maskedGridPositions;
	for (auto pos : gridPositions)
	{
		if (relevantGridPositions[idx])
			maskedGridPositions.push_back(pos);
		idx++;
	}
	*/

	// TODO: is this really better? Will be somewhat slower
	MaskedGrid maskedGridPositions(NY, NX, relevantGridPositions);

	return calculateLinearMatrix_Functors(maskedGridPositions,
			[NX,NY](auto pos) { return pos.m_i >= 0 && pos.m_i < NY && pos.m_j >= 0 && pos.m_j < NX; },
			kernel,
			//[](auto pos, auto diff) { return GridPos(pos, diff); },
			[](auto pos, auto diff) { return MaskedGridPos(pos, diff); },
			[&mpv](auto pos) { return mpv.getVariableIndexOrValue(pos.m_i, pos.m_j); },
			limitingValue,
			greaterThanLimitingValue
			);
}

MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValuesBase &mpv,
		const vector<pair<double, pair<int, int>>> &kernel)
{
	const int NX = mpv.getNX();
	const int NY = mpv.getNY();
	const int N = mpv.getNumberOfVariables();
	Grid gridPositions(NY, NX);

	return calculateQuadraticMatrix_Functors(N, gridPositions,
			[NX,NY](auto pos) { return pos.m_i >= 0 && pos.m_i < NY && pos.m_j >= 0 && pos.m_j < NX; },
			kernel,
			[](auto pos, auto diff) { return GridPos(pos, diff); },
			[&mpv](auto pos) { return mpv.getVariableIndexOrValue(pos.m_i, pos.m_j); }
			);
}

MaskedPotentialValuesOffsetGradient::MaskedPotentialValuesOffsetGradient(vector<double> &potentialValues, 
		                                                                 vector<int> &mask, int NX, double scaleUnit)
	: MaskedPotentialValuesBase(NX, potentialValues.size()/NX)
{
	assert(potentialValues.size() == mask.size());
	assert(potentialValues.size() % NX == 0);

	int maskCounts[3] = { 0, 0, 0 };
	for (auto m : mask)
	{
		if (m < 0 || m > 2)
		{
			cerr << "Internal error: illegal mask value " << m << endl;
			exit(-1);
		}
		maskCounts[m]++;
	}

	for (int i = 0 ; i < 3 ; i++)
	{
		if (!maskCounts[i])
		{
			cerr << "Internal error: not all mask values are present, can't define problem" << endl;
			exit(-1);
		}
	}

	m_scaleUnit = scaleUnit;

	m_potentialValues = move(potentialValues);
	m_mask = move(mask);

	m_idxMapInv.resize(m_potentialValues.size(), -1);
	m_idxMapFwd.push_back(-123); // room for offset in the mask == 2 section, doesn't correspond to a point in the map
	m_idxMapFwd.push_back(-123); // room for x gradient in mask == 2 section, doesn't correspond to a point in the map
	m_idxMapFwd.push_back(-123); // room for y gradient in mask == 2 section, doesn't correspond to a point in the map

	for (size_t i = 0 ; i < m_potentialValues.size() ; i++)
	{
		if (m_mask[i] == 0) // a grid value that needs to be optimized
		{
			m_idxMapInv[i] = m_idxMapFwd.size();
			m_idxMapFwd.push_back(i);
		}
	}

	setNumberOfVariables((int)m_idxMapFwd.size());
}

