#pragma once

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include "vector2d.h"
#include <errut/booltype.h>
#include <assert.h>
#include <vector>
#include <limits>
#include <string>
#include <memory>

namespace grale
{

class ImagesData;

class GRALE_IMPORTEXPORT OclCalculatedBackProjector : public ProjectedImagesInterface
{
public:
	OclCalculatedBackProjector();
	~OclCalculatedBackProjector();

	errut::bool_t init(const std::vector<ImagesDataExtended *> &images, double angularScale); 

	// TODO: use this if image position randomization is enabled; randomized positions need
	//       to be obtained from GPU - is this overkill? Is currently only used for theta
	//       differences with the retraced thetas, perhaps just calculating the retraced
	//       theta diff is better?
	//       Still, what if at some point the actual changed thetas are needed, and they
	//       will not be updated
	void setAdjustedThetas(const std::vector<Vector2Df> &adjustedThetas);
	void setRetraceInfo(const std::shared_ptr<std::vector<bool>> &tracedSourcesFlags,
			            const std::shared_ptr<std::vector<std::vector<Vector2Df>>> &tracedSourcesPoints,
						const std::shared_ptr<std::vector<std::vector<int>>> &tracedSourcesConvergedFlags)
	{
		m_tracedSourcesFlags = tracedSourcesFlags;
		m_tracedSourcesPoints = tracedSourcesPoints;
		m_tracedSourcesConvergedFlags = tracedSourcesConvergedFlags;
	}

	void setBetaBuffer(const float *pBetas, size_t s)									{ m_pBetas = pBetas; m_betasSize = s; }
	void setAlphaBuffer(const float *pAlphas, size_t s)									{ m_pAlphas = pAlphas; m_alphasSize = s; }
	void setDerivBuffers(const float *pAxx, const float *pAyy, const float *pAxy, size_t s) { m_pAxx = pAxx; m_pAyy = pAyy; m_pAxy = pAxy;  m_derivSize = s; }
	void setPotentialBuffers(const float *pPot, size_t s) { m_pPotentials = pPot; m_potSize = s; }
	
	double getLensDistance() const														{ return std::numeric_limits<double>::quiet_NaN(); }
	double getLensRedshift() const														{ return std::numeric_limits<double>::quiet_NaN(); }
	double getAngularScale() const														{ return m_angularScale; }
	const Vector2D<float> *getBetas(int sourceNumber) const;
	const Vector2D<float> *getBetas(int sourceNumber, int imageNumber) const;
	const Vector2D<float> *getThetas(int sourceNumber) const;
	const Vector2D<float> *getThetas(int sourceNumber, int imageNumber) const;
	bool hasRetracedThetas(int sourceNum) const override;
	const Vector2D<float> *getRetracedThetas(int sourceNum) const override;
	const Vector2D<float> *getRetracedThetas(int sourceNum, int imageNum) const override;
	const int *getRetracingConvergedFlags(int sourceNum) const override;
	const int *getRetracingConvergedFlags(int sourceNum, int imageNum) const override;
	const Vector2D<float> *getAlphas(int sourceNumber) const;
	const Vector2D<float> *getAlphas(int sourceNumber, int imageNumber) const;
	const float *getDerivativesXX(int sourceNumber) const;
	const float *getDerivativesXX(int sourceNumber, int imageNumber) const;
	const float *getDerivativesYY(int sourceNumber) const;
	const float *getDerivativesYY(int sourceNumber, int imageNumber) const;
	const float *getDerivativesXY(int sourceNumber) const;
	const float *getDerivativesXY(int sourceNumber, int imageNumber) const;
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

	const float *getLensPotential(int sourceNumber) const;
	const float *getLensPotential(int sourceNumber, int imageNumber) const;
	float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const { return std::numeric_limits<float>::quiet_NaN(); }
private:
	std::vector<std::vector<Vector2D<float> > > m_thetas;
	std::vector<int> m_sourceOffsets;
	const float *m_pBetas = nullptr, *m_pAlphas = nullptr;
	const float *m_pAxx = nullptr, *m_pAyy = nullptr, *m_pAxy = nullptr;
	const float *m_pPotentials = nullptr;
	size_t m_betasSize = 0, m_alphasSize = 0, m_derivSize = 0, m_potSize = 0;
	double m_angularScale = 0;

	std::shared_ptr<std::vector<bool>> m_tracedSourcesFlags;
	std::shared_ptr<std::vector<std::vector<Vector2Df>>> m_tracedSourcesPoints;
	std::shared_ptr<std::vector<std::vector<int>>> m_tracedSourcesConvergedFlags;

	template <int ptMultiplier>
	const float *getOffsetInArray(const float *pBasePtr, size_t bufSize, int sourceNumber) const
	{
		assert(sourceNumber >= 0 && sourceNumber < m_sourceOffsets.size());
		assert(pBasePtr);
		int offset = m_sourceOffsets[sourceNumber]*ptMultiplier;
		assert(offset < bufSize);
		return pBasePtr + offset;
	}

	template <int ptMultiplier>
	const float *getOffsetInArray(const float *pBasePtr, size_t bufSize, int sourceNumber, int imageNumber) const
	{
		assert(sourceNumber >= 0 && sourceNumber < m_sourceOffsets.size() && sourceNumber < m_offsets.size());
		assert(imageNumber < m_offsets[sourceNumber].size());
		assert(pBasePtr);
		int offset = (m_sourceOffsets[sourceNumber] + m_offsets[sourceNumber][imageNumber]) * ptMultiplier;
		assert(offset < bufSize);
		return pBasePtr + offset;
	}
};

inline const Vector2D<float> *OclCalculatedBackProjector::getBetas(int sourceNumber) const
{
	return reinterpret_cast<const Vector2Df*>(getOffsetInArray<2>(m_pBetas, m_betasSize, sourceNumber));
}

inline const Vector2D<float> *OclCalculatedBackProjector::getBetas(int sourceNumber, int imageNumber) const
{
	return reinterpret_cast<const Vector2Df*>(getOffsetInArray<2>(m_pBetas, m_betasSize, sourceNumber, imageNumber));
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

inline const Vector2D<float> *OclCalculatedBackProjector::getAlphas(int sourceNumber) const
{
	return reinterpret_cast<const Vector2Df*>(getOffsetInArray<2>(m_pAlphas, m_alphasSize, sourceNumber));
}

inline const Vector2D<float> *OclCalculatedBackProjector::getAlphas(int sourceNumber, int imageNumber) const
{
	return reinterpret_cast<const Vector2Df*>(getOffsetInArray<2>(m_pAlphas, m_alphasSize, sourceNumber, imageNumber));
}

inline const float *OclCalculatedBackProjector::getDerivativesXX(int sourceNumber) const
{
	return getOffsetInArray<1>(m_pAxx, m_derivSize, sourceNumber);
}

inline const float *OclCalculatedBackProjector::getDerivativesXX(int sourceNumber, int imageNumber) const
{
	return getOffsetInArray<1>(m_pAxx, m_derivSize, sourceNumber, imageNumber);
}

inline const float *OclCalculatedBackProjector::getDerivativesYY(int sourceNumber) const
{
	return getOffsetInArray<1>(m_pAyy, m_derivSize, sourceNumber);
}

inline const float *OclCalculatedBackProjector::getDerivativesYY(int sourceNumber, int imageNumber) const
{
	return getOffsetInArray<1>(m_pAyy, m_derivSize, sourceNumber, imageNumber);
}

inline const float *OclCalculatedBackProjector::getDerivativesXY(int sourceNumber) const
{
	return getOffsetInArray<1>(m_pAxy, m_derivSize, sourceNumber);
}

inline const float *OclCalculatedBackProjector::getDerivativesXY(int sourceNumber, int imageNumber) const
{
	return getOffsetInArray<1>(m_pAxy, m_derivSize, sourceNumber, imageNumber);
}

inline const float *OclCalculatedBackProjector::getLensPotential(int sourceNumber) const
{
	return getOffsetInArray<1>(m_pPotentials, m_potSize, sourceNumber);
}

inline const float *OclCalculatedBackProjector::getLensPotential(int sourceNumber, int imageNumber) const
{
	return getOffsetInArray<1>(m_pPotentials, m_potSize, sourceNumber, imageNumber);
}

inline bool OclCalculatedBackProjector::hasRetracedThetas(int sourceNum) const
{
	if (!m_tracedSourcesPoints)
		return false;

	assert(sourceNum < m_tracedSourcesFlags->size());
	return (*m_tracedSourcesFlags)[sourceNum];
}

inline const Vector2D<float> *OclCalculatedBackProjector::getRetracedThetas(int sourceNum) const
{
	assert(m_tracedSourcesPoints.get());
	assert(sourceNum < m_tracedSourcesPoints->size());
	const std::vector<Vector2Df> &allSrcPoints = (*m_tracedSourcesPoints)[sourceNum];
	assert(allSrcPoints.size() == getNumberOfImagePoints(sourceNum));
	return allSrcPoints.data();
}

inline const Vector2D<float> *OclCalculatedBackProjector::getRetracedThetas(int sourceNum, int imageNum) const
{
	assert(m_tracedSourcesPoints.get());
	assert(sourceNum < m_tracedSourcesPoints->size());
	const std::vector<Vector2Df> &allSrcPoints = (*m_tracedSourcesPoints)[sourceNum];
	assert(allSrcPoints.size() == getNumberOfImagePoints(sourceNum));

	size_t imgOff = m_offsets[sourceNum][imageNum];
	assert(imgOff + getNumberOfImagePoints(sourceNum, imageNum) <= allSrcPoints.size());

	return allSrcPoints.data() + imgOff;
}

inline const int *OclCalculatedBackProjector::getRetracingConvergedFlags(int sourceNum) const
{
	assert(m_tracedSourcesConvergedFlags.get());
	assert(sourceNum < m_tracedSourcesConvergedFlags->size());
	const std::vector<int> &allSrcConvFlags= (*m_tracedSourcesConvergedFlags)[sourceNum];
	assert(allSrcConvFlags.size() == getNumberOfImagePoints(sourceNum));
	return allSrcConvFlags.data();
}

inline const int *OclCalculatedBackProjector::getRetracingConvergedFlags(int sourceNum, int imageNum) const
{
	assert(m_tracedSourcesConvergedFlags.get());
	assert(sourceNum < m_tracedSourcesConvergedFlags->size());
	const std::vector<int> &allSrcConvFlags = (*m_tracedSourcesConvergedFlags)[sourceNum];
	assert(allSrcConvFlags.size() == getNumberOfImagePoints(sourceNum));

	size_t imgOff = m_offsets[sourceNum][imageNum];
	assert(imgOff + getNumberOfImagePoints(sourceNum, imageNum) <= allSrcConvFlags.size());

	return allSrcConvFlags.data() + imgOff;
}

inline void OclCalculatedBackProjector::setAdjustedThetas(const std::vector<Vector2Df> &adjustedThetas)
{
	size_t idx = 0;
	for (auto &src : m_thetas)
	{
		for (auto &t : src)
		{
			assert(idx < adjustedThetas.size());
			t = adjustedThetas[idx];
			idx++;
		}
	}

	assert(idx == adjustedThetas.size());
}

} // end namespace

