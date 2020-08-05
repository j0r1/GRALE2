#pragma once

/**
 * \file precalculatedbackprojector.h
 */

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include "vector2d.h"
#include "imagesdataextended.h"
#include "plummerlensinfo.h"
#include <errut/errorbase.h>
#include <assert.h>
#include <vector>
#include <limits>
#include <string>
#include <memory>

namespace grale
{

class MultiPlaneCUDA;
class Cosmology;

class GRALE_IMPORTEXPORT MPCUDABackProjector : public ProjectedImagesInterface, public errut::ErrorBase
{
public:
	MPCUDABackProjector();
	~MPCUDABackProjector();

    bool init(const std::string &libraryPath, int deviceIndex,
		const Cosmology &cosmology,
		const std::vector<float> &lensRedshifts,
		const std::vector<std::vector<PlummerLensInfo>> &lenses, 
		const std::vector<float> &sourceRedshifts,
		const std::vector<ImagesDataExtended *> &images);

	bool calculateSourcePositions(const std::vector<std::vector<float>> &massFactors,
	                              const std::vector<float> &sheetDensities);

    // No single lens distance is available in the multi-plane case
	double getLensDistance() const														{ return std::numeric_limits<double>::quiet_NaN(); }
    // No single lens redshift is available in the multi-plane case
	double getLensRedshift() const														{ return std::numeric_limits<double>::quiet_NaN(); }

	double getAngularScale() const														{ return m_angularScale; }
	const Vector2D<float> *getBetas(int sourceNumber) const;
	const Vector2D<float> *getBetas(int sourceNumber, int imageNumber) const;
	const Vector2D<float> *getThetas(int sourceNumber) const;
	const Vector2D<float> *getThetas(int sourceNumber, int imageNumber) const;
	const Vector2D<float> *getAlphas(int sourceNumber) const							{ return nullptr; }
	const Vector2D<float> *getAlphas(int sourceNumber, int imageNumber) const			{ return nullptr; }
	const float *getDerivativesXX(int sourceNumber) const								{ return nullptr; }
	const float *getDerivativesXX(int sourceNumber, int imageNumber) const				{ return nullptr; }
	const float *getDerivativesYY(int sourceNumber) const								{ return nullptr; }
	const float *getDerivativesYY(int sourceNumber, int imageNumber) const				{ return nullptr; }
	const float *getDerivativesXY(int sourceNumber) const								{ return nullptr; }
	const float *getDerivativesXY(int sourceNumber, int imageNumber) const				{ return nullptr; }
	const float *getInverseMagnifications(int sourceNumber) const						{ return nullptr; }
	const float *getInverseMagnifications(int sourceNumber, int imageNumber) const		{ return nullptr; }
	const float *getShearComponents1(int sourceNumber) const							{ return nullptr; }
	const float *getShearComponents1(int sourceNumber, int imageNumber) const			{ return nullptr; }
	const float *getShearComponents2(int sourceNumber) const							{ return nullptr; }
	const float *getShearComponents2(int sourceNumber, int imageNumber) const			{ return nullptr; }
	const float *getConvergence(int sourceNumber) const									{ return nullptr; }
	const float *getConvergence(int sourceNumber, int imageNumber) const				{ return nullptr; }

	const float *getLensPotential(int sourceNumber) const								{ return nullptr; }
	const float *getLensPotential(int sourceNumber, int imageNumber) const				{ return nullptr; }
	float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const { return std::numeric_limits<float>::quiet_NaN(); }
private:
    static double getAngularScale(const std::vector<ImagesDataExtended*> &images);

    MultiPlaneCUDA *m_pMPCU;
	std::vector<std::vector<Vector2D<float>>> m_thetas;
	double m_angularScale;
};

inline const Vector2D<float> *MPCUDABackProjector::getThetas(int sourceNumber) const
{
	assert(sourceNumber < m_thetas.size());
	return &(m_thetas[sourceNumber][0]);
}

inline const Vector2D<float> *MPCUDABackProjector::getThetas(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber < m_thetas.size() && sourceNumber < m_offsets.size());
	assert(imageNumber < m_offsets[sourceNumber].size());

	int offset = m_offsets[sourceNumber][imageNumber];
	assert(offset < m_thetas[sourceNumber].size());
	return &(m_thetas[sourceNumber][offset]);
}

} // end namespace

