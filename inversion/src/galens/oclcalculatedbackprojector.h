#pragma once

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include "vector2d.h"
#include <errut/booltype.h>
#include <assert.h>
#include <vector>
#include <limits>
#include <string>

namespace grale
{

class ImagesData;

class GRALE_IMPORTEXPORT OclCalculatedBackProjector : public ProjectedImagesInterface
{
public:
	OclCalculatedBackProjector();
	~OclCalculatedBackProjector();

	errut::bool_t init(const std::vector<ImagesDataExtended *> &images, double angularScale); 
	void setBetaBuffer(const float *pBetas, size_t s)									{ m_pBetas = pBetas; m_betasSize = s; }
	
	double getLensDistance() const														{ return std::numeric_limits<double>::quiet_NaN(); }
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
	const float *getSecondDerivativesXXX(int sourceNumber) const						{ return nullptr; }
	const float *getSecondDerivativesXXX(int sourceNumber, int imageNumber) const		{ return nullptr; }
	const float *getSecondDerivativesYYY(int sourceNumber) const						{ return nullptr; }
	const float *getSecondDerivativesYYY(int sourceNumber, int imageNumber) const		{ return nullptr; }
	const float *getSecondDerivativesXXY(int sourceNumber) const						{ return nullptr; }
	const float *getSecondDerivativesXXY(int sourceNumber, int imageNumber) const		{ return nullptr; }
	const float *getSecondDerivativesYYX(int sourceNumber) const						{ return nullptr; }
	const float *getSecondDerivativesYYX(int sourceNumber, int imageNumber) const		{ return nullptr; }
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
	std::vector<std::vector<Vector2D<float> > > m_thetas;
	std::vector<int> m_sourceOffsets;
	const float *m_pBetas;
	size_t m_betasSize;
	double m_angularScale;
};

inline const Vector2D<float> *OclCalculatedBackProjector::getBetas(int sourceNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_sourceOffsets.size());
	assert(m_pBetas);
	int offset = m_sourceOffsets[sourceNumber];
	assert(offset < m_betasSize);
	return reinterpret_cast<const Vector2Df*>(m_pBetas + offset);
}

inline const Vector2D<float> *OclCalculatedBackProjector::getBetas(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_sourceOffsets.size() && sourceNumber < m_offsets.size());
	assert(imageNumber < m_offsets[sourceNumber].size());
	assert(m_pBetas);
	int offset = m_sourceOffsets[sourceNumber] + m_offsets[sourceNumber][imageNumber] * 2;
	assert(offset < m_betasSize);
	return reinterpret_cast<const Vector2Df*>(m_pBetas + offset);
}

inline const Vector2D<float> *OclCalculatedBackProjector::getThetas(int sourceNumber) const
{
	assert(sourceNumber < m_thetas.size());
	return &(m_thetas[sourceNumber][0]);
}

inline const Vector2D<float> *OclCalculatedBackProjector::getThetas(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber < m_thetas.size() && sourceNumber < m_offsets.size());
	assert(imageNumber < m_offsets[sourceNumber].size());

	int offset = m_offsets[sourceNumber][imageNumber];
	assert(offset < m_thetas[sourceNumber].size());
	return &(m_thetas[sourceNumber][offset]);
}

} // end namespace

