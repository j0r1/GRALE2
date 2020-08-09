#ifndef MULTIPLANECUDA_H

#define MULTIPLANECUDA_H

#include "graleconfig.h"
#include "vector2d.h"
#include <errut/errorbase.h>
#include <vector>
#include <memory>

namespace grale
{

class PerNodeCounter;

class GRALE_IMPORTEXPORT MultiPlaneCUDA : public errut::ErrorBase
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

	// Device idx -1 means get next device
	bool init(const std::string &libraryPath, int deviceIdx,
		double angularUnit,
		double h, double W_m, double W_r, double W_v, double w,
		const std::vector<float> &lensRedshifts,
		const std::vector<std::vector<PlummerInfo>> &fixedPlummerParameters, 
		const std::vector<float> &sourceRedshifts,
		const std::vector<std::vector<Vector2Df>> &theta);

	bool calculateSourcePositions(const std::vector<std::vector<float>> &massFactors,
	                              const std::vector<float> &sheetDensities = std::vector<float>());
	const std::vector<Vector2Df> *getSourcePositions(int srcIdx);
	int getDeviceIndex() const { return m_deviceIndex; }
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
		void **pCtx,
		int deviceIdx);

	int (*mpcuCalculateSourcePositions)(void *ctx, const std::vector<std::vector<float>> &massFactors,
	                                    const std::vector<float> &sheetDensities);
	const std::vector<Vector2Df> & (*mpcuGetSourcePositions)(void *ctx, int srcIdx);
	void (*mpcuClearContext)(void *ctx);
	int (*mpcuGetNumberOfDevices)(void);

	void *m_pLibrary;
	void *m_pContext;
	int m_numSourcePlanes;
	std::unique_ptr<PerNodeCounter> m_perNodeCounter;
	int m_deviceIndex;
};

} // end namespace

#endif // MULTIPLANECUDA_H
