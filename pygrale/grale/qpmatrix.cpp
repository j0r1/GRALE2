#include "qpmatrix_templates.h"
#include <cassert>
#include <vector>
#include <iostream>

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

MaskedPotentialValues::MaskedPotentialValues(vector<double> &potentialValues, vector<bool> &mask, int NX, double scaleUnit)
{
	assert(potentialValues.size() == mask.size());
	assert(potentialValues.size() % NX == 0);
	
	m_scaleUnit = scaleUnit;
	m_NY = potentialValues.size()/NX;
	m_NX = NX;
	assert(m_NX > 0 && m_NY > 0);

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
}

MatrixResults calculateLinearConstraintMatrices(const MaskedPotentialValues &mpv,
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
			[&mpv](auto pos) { return array<tuple<int,double,double>,1> { mpv.getVariableIndexOrValue(pos.m_i, pos.m_j) }; }
			);
}

MatrixResults calculateQuadraticMimimizationMatrices(const MaskedPotentialValues &mpv,
		const vector<pair<double,vector<pair<double, pair<int, int>>>>> &kernelList
		)
{
	const int NX = mpv.getNX();
	const int NY = mpv.getNY();
	const int N = mpv.getNumberOfVariables();
	Grid gridPositions(NY, NX);

	return calculateQuadraticMatrix_Functors(N, gridPositions,
			[NX,NY](auto pos) { return pos.m_i >= 0 && pos.m_i < NY && pos.m_j >= 0 && pos.m_j < NX; },
			kernelList,
			[](auto pos, auto diff) { return GridPos(pos, diff); },
			[&mpv](auto pos) { return array<tuple<int,double,double>,1> { mpv.getVariableIndexOrValue(pos.m_i, pos.m_j) }; }
			);
}

