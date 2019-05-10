#ifndef MULTIPLANECUDA_H

#define MULTIPLANECUDA_H

#include <errut/errorbase.h>
#include <vector>
#include "vector2d.h"

namespace grale
{

class MultiPlaneCUDA : public errut::ErrorBase
{
public:
	struct PlummerInfo
	{
	public:
		PlummerInfo(Vector2Df position = Vector2Df(), float widthInAngularUnit = 0, double initialMass = 0) : 
			m_position(position), m_widthInAngularUnit(widthInAngularUnit), m_initialMass(initialMass) { }

		Vector2Df m_position;
		float m_widthInAngularUnit;
		double m_initialMass;
	};

	MultiPlaneCUDA();
	~MultiPlaneCUDA();

	bool init(const std::string &libraryPath,
		double angularUnit,
		double h, double W_m, double W_r, double W_v, double w,
		const std::vector<float> &lensRedshifts,
		const std::vector<std::vector<PlummerInfo>> &fixedPlummerParameters, 
		const std::vector<float> &sourceRedshifts,
		const std::vector<std::vector<Vector2Df>> &theta);

	bool calculateSourcePositions(const std::vector<std::vector<float>> &massFactors);
	const std::vector<Vector2Df> *getSourcePositions(int srcIdx);
private:
	void zero();
	void cleanup();

	int (*mpcuInitMultiPlaneCalculation)(
		double angularUnit,
		double h, double W_m, double W_r, double W_v, double w,
		const std::vector<float> &lensRedshifts,
		const std::vector<std::vector<PlummerInfo>> &fixedPlummerParameters, 
		const std::vector<float> &sourceRedshifts,
		const std::vector<std::vector<Vector2Df>> &theta, 
		void **pCtx);

	int (*mpcuCalculateSourcePositions)(void *ctx, const std::vector<std::vector<float>> &massFactors);
	const std::vector<Vector2Df> & (*mpcuGetSourcePositions)(void *ctx, int srcIdx);
	void (*mpcuClearContext)(void *ctx);

	void *m_pLibrary;
	void *m_pContext;
	int m_numSourcePlanes;
};

} // end namespace

#endif // MULTIPLANECUDA_H
