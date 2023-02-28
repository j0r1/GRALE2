#pragma once

#include "qpmatrix.h"

template<class Positions, typename PositionValidatorFunctor, class Kernel, typename PositionCombinerFunctor, typename PositionMapperFunctor>
MatrixResults calculateLinearMatrix_Functors(const Positions &positions, PositionValidatorFunctor posVal, 
		                            const Kernel &coefficients, PositionCombinerFunctor combiner,
									PositionMapperFunctor mapping)
{
    MatrixResults results;
    SparseMatrixInfo &sparse = results.first;
    std::vector<double> &colMatrix = results.second;

    auto &b = colMatrix;
	auto &rows = sparse.m_rows;
	auto &cols = sparse.m_cols;
	auto &data = sparse.m_values;
	int rowIdx = 0;

	for (auto pos : positions)
	{
		double bValue = 0;
		std::unordered_map<int, double> Avalues;

		bool valid = true;
		for (auto coeff_diffPos : coefficients)
		{
			auto coeff = coeff_diffPos.first;
			auto diffPos = coeff_diffPos.second;

			auto newPos = combiner(pos, diffPos);
			if (!posVal(newPos))
			{
				valid = false;
				break;
			}

			for (auto idx_value_extra_factor : mapping(newPos))
			{
				auto idx = std::get<0>(idx_value_extra_factor);
				auto value = std::get<1>(idx_value_extra_factor);
				auto extra_factor = std::get<2>(idx_value_extra_factor);

				if (idx < 0) // Don't optimize, move to constraint
					bValue += -coeff * extra_factor * value;
				else
				{
					auto it = Avalues.find(idx);
					if (it == Avalues.end())
						Avalues[idx] = coeff*extra_factor;
					else
						it->second += coeff*extra_factor;
				}
			}
		}

		if (valid && Avalues.size() > 0)
		{
			b.push_back(bValue);

			for (auto k_v : Avalues)
			{
				auto k = k_v.first;
				auto v = k_v.second;

				rows.push_back(rowIdx);
				cols.push_back(k);
				data.push_back(v);
			}

			rowIdx++;
		}
	}

    return results;
}

template<class Positions, class PositionValidator, class Kernel, class PositionCombiner, class PositionMapper>
MatrixResults calculateLinearMatrix(const Positions &positions, const PositionValidator &posVal, 
		                            const Kernel &coefficients, const PositionCombiner &combiner,
									const PositionMapper &mapping)
{
	return calculateLinearMatrix_Functors(positions, 
			[&posVal](auto pos){ return posVal.isPositionValid(pos); },
			coefficients,
			[&combiner](auto pos, auto diff) { return combiner.combinePositionAndDiff(pos, diff); },
			[&mapping](auto pos) { return mapping.getVariableIndexOrValue(pos); });
}

template<class Positions, typename PositionValidatorFunctor, class KernelList, typename PositionCombinerFunctor, typename PositionMapperFunctor>
MatrixResults calculateQuadraticMatrix_Functors(int N, const Positions &positions, PositionValidatorFunctor posVal, 
		                            const KernelList &kernelList, PositionCombinerFunctor combiner,
									PositionMapperFunctor mapping)
{
    MatrixResults results;
	auto &sparse = results.first;
	auto &q = results.second;
	q.assign(N, 0);

	std::vector<std::unordered_map<int, double>> P(N);

	std::vector<std::pair<int, double>> qBuffer;
	std::vector<std::tuple<int,int, double>> PBuffer;

	for (auto pos : positions)
	{
		for (auto &weigthKernel : kernelList)
		{
			auto &weight = weigthKernel.first;
			auto &kernel = weigthKernel.second;

			qBuffer.resize(0);
			PBuffer.resize(0);

			bool valid = true;

			for (auto factor1_diffPos1 : kernel)
			{
				auto factor1 = factor1_diffPos1.first;
				auto diffPos1 = factor1_diffPos1.second;

				//cout << factor1 << " (" << diffPos1.first << "," << diffPos1.second << ")" << endl;
				auto fullPos1 = combiner(pos, diffPos1);
				if (!posVal(fullPos1))
				{
					valid = false;
					break;
				}

				for (auto idx1_value1_extra_factor : mapping(fullPos1))
				{
					auto idx1 = std::get<0>(idx1_value1_extra_factor);
					auto value1 = std::get<1>(idx1_value1_extra_factor);
					auto extra_factor1 = std::get<2>(idx1_value1_extra_factor);

					for (auto factor2_diffPos2 : kernel)
					{
						auto factor2 = factor2_diffPos2.first;
						auto diffPos2 = factor2_diffPos2.second;

						auto fullPos2 = combiner(pos, diffPos2);
						if (!posVal(fullPos2))
						{
							valid = false;
							break;
						}

						for (auto idx2_value2_extra_factor : mapping(fullPos2))
						{
							auto idx2 = std::get<0>(idx2_value2_extra_factor);
							auto value2 = std::get<1>(idx2_value2_extra_factor);
							auto extra_factor2 = std::get<2>(idx2_value2_extra_factor);

							if (idx1 < 0)
							{
								if (idx2 < 0)
								{
									// Nothing to do
								}
								else
									qBuffer.push_back({idx2,weight*factor1*extra_factor1*factor2*extra_factor2*value1});
							}
							else // idx1 >= 0
							{
								if (idx2 < 0)
									qBuffer.push_back({idx1,weight*factor1*extra_factor1*factor2*extra_factor2*value2});
								else
									PBuffer.push_back({idx1, idx2, weight*factor1*extra_factor1*factor2*extra_factor2});
							}
						}
					}

					if (!valid)
						break;
				}
				if (!valid)
					break;
			}
			
			if (valid)
			{
				for (auto k_v : qBuffer)
				{
					auto k = k_v.first;
					auto v = k_v.second;

					assert(k >= 0 && k < q.size());
					q[k] += v;
				}

				for (auto i_j_v : PBuffer)
				{
					auto i = std::get<0>(i_j_v);
					auto j = std::get<1>(i_j_v);
					auto v = std::get<2>(i_j_v);

					assert(i >= 0 && i < P.size());
					assert(j >= 0 && j < P.size());

					auto it = P[i].find(j);
					if (it == P[i].end())
						P[i][j] = v;
					else
						it->second += v;
				}
			}
		}
	}

	auto &data = sparse.m_values;
	auto &cols = sparse.m_cols;
	auto &rows = sparse.m_rows;

	for (int r = 0 ; r < (int)P.size() ; r++)
	{
		for (auto idx_value : P[r])
		{
			auto idx = idx_value.first;
			auto value = idx_value.second;

			data.push_back(value * 2.0); //  Will use 1/2 * x^T * P * x
			rows.push_back(r);
			cols.push_back(idx);
		}
	}

	return results;
}

template<class Positions, class PositionValidator, class KernelList, class PositionCombiner, class PositionMapper>
MatrixResults calculateQuadraticMatrix(int N, const Positions &positions, const PositionValidator &posVal, 
		                            const KernelList &kernelList, const PositionCombiner &combiner,
									const PositionMapper &mapping)
{
	return calculateQuadraticMatrix_Functors(N, positions, 
			[&posVal](auto pos){ return posVal.isPositionValid(pos); },
			kernelList,
			[&combiner](auto pos, auto diff) { return combiner.combinePositionAndDiff(pos, diff); },
			[&mapping](auto pos) { return mapping.getVariableIndexOrValue(pos); });
}


