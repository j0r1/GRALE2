#pragma once

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include <errut/booltype.h>
#include <memory>
#include <stdexcept>

namespace grale
{

class PositionRandomizationBackprojectWrapper : public ProjectedImagesInterface
{
public:
	PositionRandomizationBackprojectWrapper(const std::shared_ptr<ProjectedImagesInterface> &baseBackProjector) : m_baseBp(baseBackProjector) { }
	~PositionRandomizationBackprojectWrapper() { }

	errut::bool_t initRandomization(const std::vector<bool> srcHasRandomization, size_t &numRandomizablePoints);
	errut::bool_t setRandomOffsets(const std::vector<Vector2Df> &offsets);
	void clearCachedValues();

	int getNumberOfSources() const override { return m_baseBp->getNumberOfSources(); }
	int getNumberOfImages(int sourceNumber) const override { return m_baseBp->getNumberOfImages(sourceNumber); }
	int getNumberOfImagePoints(int sourceNumber) const override { return m_baseBp->getNumberOfImagePoints(sourceNumber); }
	int getNumberOfImagePoints(int sourceNumber, int imageNumber) const override { return m_baseBp->getNumberOfImagePoints(sourceNumber, imageNumber); }
	bool hasOriginalProperty(ImagesData::PropertyName n, int sourceNumber) const override { return m_baseBp->hasOriginalProperty(n, sourceNumber); }
	const float *getOriginalProperties(ImagesData::PropertyName n, int sourceNumber) const override { return m_baseBp->getOriginalProperties(n, sourceNumber); }
	const float *getOriginalProperties(ImagesData::PropertyName n, int sourceNumber, int imageNumber) const override { return m_baseBp->getOriginalProperties(n, sourceNumber, imageNumber); }
	int getOriginalNumberOfTimeDelays(int sourceNumber) const override { return m_baseBp->getOriginalNumberOfTimeDelays(sourceNumber); }
	void getOriginalTimeDelay(int sourceNumber, int index, int *pImg, int *pPoint, float *pDelay) const override { return m_baseBp->getOriginalTimeDelay(sourceNumber, index, pImg, pPoint, pDelay); }
	float getDistanceFraction(int sourcenum) const override { return m_baseBp->getDistanceFraction(sourcenum); }

	double getLensDistance() const override { return m_baseBp->getLensDistance(); }
	double getLensRedshift() const override { return m_baseBp->getLensRedshift(); }
	double getAngularScale() const override { return m_baseBp->getAngularScale(); };

	bool hasRandomization(int sourceNum) const
	{
		assert(sourceNum >= 0 && sourceNum < m_srcHaveUncerts.size());
		return m_srcHaveUncerts[sourceNum];
	}
	
	void checkHaveNoRandomization(int sourceNum) const
	{
		if (hasRandomization(sourceNum))
			throw std::runtime_error("Cannot use call with image positions that have random offsets");
	}

	const Vector2D<float> *getBetas(int sourcenum) const override;
	const Vector2D<float> *getBetas(int sourcenum, int imagenum) const override;
	const Vector2D<float> *getAlphas(int sourcenum) const override;
	const Vector2D<float> *getAlphas(int sourcenum, int imagenum) const override;
	const Vector2D<float> *getThetas(int sourcenum) const override;
	const Vector2D<float> *getThetas(int sourcenum, int imagenum) const override;
	const float *getLensPotential(int sourceNumber) const override;
	const float *getLensPotential(int sourceNumber, int imageNumber) const override;

	// These should not be used with randomized positions, just forward them
	const float *getDerivativesXX(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getDerivativesXX(sourceNumber); }
	const float *getDerivativesXX(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getDerivativesXX(sourceNumber, imageNumber); }
	const float *getDerivativesYY(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getDerivativesYY(sourceNumber); }
	const float *getDerivativesYY(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getDerivativesYY(sourceNumber, imageNumber); }
	const float *getDerivativesXY(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getDerivativesXY(sourceNumber); }
	const float *getDerivativesXY(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getDerivativesXY(sourceNumber, imageNumber); }
	const float *getSecondDerivativesXXX(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesXXX(sourceNumber); }
	const float *getSecondDerivativesXXX(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesXXX(sourceNumber, imageNumber); }
	const float *getSecondDerivativesYYY(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesYYY(sourceNumber); }
	const float *getSecondDerivativesYYY(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesYYY(sourceNumber, imageNumber); }
	const float *getSecondDerivativesXXY(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesXXY(sourceNumber); }
	const float *getSecondDerivativesXXY(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesXXY(sourceNumber, imageNumber); }
	const float *getSecondDerivativesYYX(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesYYX(sourceNumber); }
	const float *getSecondDerivativesYYX(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getSecondDerivativesYYX(sourceNumber, imageNumber); }
	const float *getInverseMagnifications(int sourcenum) const override { checkHaveNoRandomization(sourcenum); return m_baseBp->getInverseMagnifications(sourcenum); }
	const float *getInverseMagnifications(int sourcenum, int imagenum) const override { checkHaveNoRandomization(sourcenum); return m_baseBp->getInverseMagnifications(sourcenum, imagenum); }
	const float *getShearComponents1(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getShearComponents1(sourceNumber); }
	const float *getShearComponents1(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getShearComponents1(sourceNumber, imageNumber); }
	const float *getShearComponents2(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getShearComponents2(sourceNumber); }
	const float *getShearComponents2(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getShearComponents2(sourceNumber, imageNumber); }
	const float *getConvergence(int sourceNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getConvergence(sourceNumber); }
	const float *getConvergence(int sourceNumber, int imageNumber) const override { checkHaveNoRandomization(sourceNumber); return m_baseBp->getConvergence(sourceNumber, imageNumber); }

	// TODO: remove this, or at least as virtual function, can be calculated from rest?
	float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const override
	{
		throw std::runtime_error("getTimedelay is currently not supported for case with random image offsets");
	}
private:
	const Vector2D<float> *getOriginalAlphas(int sourcenum) const { return m_baseBp->getAlphas(sourcenum); }
	const Vector2D<float> *getOriginalAlphas(int sourcenum, int imagenum) const { return m_baseBp->getAlphas(sourcenum, imagenum); }

	void checkAdjustedThetas(int sourceNum) const
	{
		assert(sourceNum >= 0 && sourceNum < m_srcAdjustedThetas.size());
		assert(m_srcHaveUncerts[sourceNum]);
		if (m_srcAdjustedThetas[sourceNum].size() != 0) // already calculated, reuse
			return;
		
		assert(sourceNum < m_srcImagePositionDifferences.size());
		int numPts = (int)m_srcImagePositionDifferences[sourceNum].size();
		assert(numPts == m_baseBp->getNumberOfImagePoints(sourceNum));

		const Vector2Df *pOrigThetas = m_baseBp->getThetas(sourceNum);
		const Vector2Df *pDiffThetas = m_srcImagePositionDifferences[sourceNum].data();
		for (size_t i = 0 ; i < numPts ; i++)
			m_srcAdjustedThetas[sourceNum].push_back(pDiffThetas[i] + pOrigThetas[i]);
	}

	void checkAdjustedAlphas(int sourceNum) const
	{
		assert(sourceNum >= 0 && sourceNum < m_srcAdjustedAlphas.size());
		assert(m_srcHaveUncerts[sourceNum]);
		if (m_srcAdjustedAlphas[sourceNum].size() != 0) // already calculated, reuse
			return;
		
		assert(sourceNum < m_srcImagePositionDifferences.size());
		int numPts = (int)m_srcImagePositionDifferences[sourceNum].size();
		assert(numPts == m_baseBp->getNumberOfImagePoints(sourceNum));

		const Vector2Df *pOrigAlphas = m_baseBp->getAlphas(sourceNum);
		const Vector2Df *pDiffThetas = m_srcImagePositionDifferences[sourceNum].data();
		const float *pAxx = m_baseBp->getDerivativesXX(sourceNum);
		const float *pAyy = m_baseBp->getDerivativesYY(sourceNum);
		const float *pAxy = m_baseBp->getDerivativesXY(sourceNum);

		for (size_t i = 0 ; i < numPts ; i++)
		{
			Vector2Df alpha = pOrigAlphas[i];
			float dx = pDiffThetas[i].getX();
			float dy = pDiffThetas[i].getY();
			float axx = pAxx[i];
			float ayy = pAyy[i];
			float axy = pAxy[i];

			float newAx = alpha.getX() + axx*dx + axy*dy;
			float newAy = alpha.getY() + axy*dx + ayy*dy;

			m_srcAdjustedAlphas[sourceNum].push_back({newAx, newAy});
		}
	}

	void checkAdjustedBetas(int sourceNum) const
	{
		assert(sourceNum >= 0 && sourceNum < m_srcAdjustedBetas.size());
		assert(m_srcHaveUncerts[sourceNum]);
		if (m_srcAdjustedBetas[sourceNum].size() != 0) // already calculated, reuse
			return;
		
		assert(sourceNum < m_srcImagePositionDifferences.size());
		int numPts = (int)m_srcImagePositionDifferences[sourceNum].size();
		assert(numPts == m_baseBp->getNumberOfImagePoints(sourceNum));

		const Vector2Df *pOrigBetas = m_baseBp->getBetas(sourceNum);
		const Vector2Df *pDiffThetas = m_srcImagePositionDifferences[sourceNum].data();
		const float *pAxx = m_baseBp->getDerivativesXX(sourceNum);
		const float *pAyy = m_baseBp->getDerivativesYY(sourceNum);
		const float *pAxy = m_baseBp->getDerivativesXY(sourceNum);

		for (size_t i = 0 ; i < numPts ; i++)
		{
			Vector2Df beta = pOrigBetas[i];
			float dx = pDiffThetas[i].getX();
			float dy = pDiffThetas[i].getY();
			float axx = pAxx[i];
			float ayy = pAyy[i];
			float axy = pAxy[i];

			float newBx = beta.getX() + (1.0f-axx)*dx - axy*dy;
			float newBy = beta.getY() - axy*dx + (1.0f-ayy)*dy;

			m_srcAdjustedBetas[sourceNum].push_back({newBx, newBy});
		}
	}

	void checkAdjustedPotentials(int sourceNum) const
	{
		assert(sourceNum >= 0 && sourceNum < m_srcAdjustedPotentials.size());
		assert(m_srcHaveUncerts[sourceNum]);
		if (m_srcAdjustedPotentials[sourceNum].size() != 0) // already calculated, reuse
			return;
		
		assert(sourceNum < m_srcImagePositionDifferences.size());
		int numPts = (int)m_srcImagePositionDifferences[sourceNum].size();
		assert(numPts == m_baseBp->getNumberOfImagePoints(sourceNum));

		const float *pOrigPotentials = m_baseBp->getLensPotential(sourceNum);
		const Vector2Df *pOrigAlphas = m_baseBp->getAlphas(sourceNum);
		const Vector2Df *pDiffThetas = m_srcImagePositionDifferences[sourceNum].data();
		const float *pAxx = m_baseBp->getDerivativesXX(sourceNum);
		const float *pAyy = m_baseBp->getDerivativesYY(sourceNum);
		const float *pAxy = m_baseBp->getDerivativesXY(sourceNum);

		for (size_t i = 0 ; i < numPts ; i++)
		{
			float phi = pOrigPotentials[i];
			Vector2Df alpha = pOrigAlphas[i];
			float dx = pDiffThetas[i].getX();
			float dy = pDiffThetas[i].getY();
			float axx = pAxx[i];
			float ayy = pAyy[i];
			float axy = pAxy[i];

			float newPhi = phi + alpha.getX()*dx + alpha.getY()*dy
			               + 0.5f*(axx*dx*dx + ayy*dy*dy + 2.0f*axy*dx*dy);
			m_srcAdjustedPotentials[sourceNum].push_back(newPhi);
		}
	}

	std::shared_ptr<ProjectedImagesInterface> m_baseBp;
	std::vector<bool> m_srcHaveUncerts;
	
	std::vector<std::vector<Vector2Df>> m_srcImagePositionDifferences;
	std::vector<std::vector<int>> m_srcImageOffsets;

	mutable std::vector<std::vector<Vector2Df>> m_srcAdjustedThetas;
	mutable std::vector<std::vector<Vector2Df>> m_srcAdjustedAlphas;
	mutable std::vector<std::vector<Vector2Df>> m_srcAdjustedBetas;
	mutable std::vector<std::vector<float>> m_srcAdjustedPotentials;
};

inline errut::bool_t PositionRandomizationBackprojectWrapper::initRandomization(const std::vector<bool> srcHasRandomization, size_t &numRandomizablePoints)
{
	if (srcHasRandomization.size() != m_baseBp->getNumberOfSources())
		return "Vector length (" + std::to_string(srcHasRandomization.size()) + ") is incompatible with number of sources (" + std::to_string(m_baseBp->getNumberOfSources()) + ")";

	m_srcHaveUncerts = srcHasRandomization;
	size_t numSources = srcHasRandomization.size();

	auto clearEntries = [numSources](auto &vec)
	{
		vec.resize(numSources);
		for (auto &x: vec)
			x.clear();
	};

	clearEntries(m_srcImagePositionDifferences);
	clearEntries(m_srcImageOffsets);
	clearEntries(m_srcAdjustedThetas);
	clearEntries(m_srcAdjustedAlphas);
	clearEntries(m_srcAdjustedBetas);
	clearEntries(m_srcAdjustedPotentials);

	numRandomizablePoints = 0;
	for (size_t s = 0 ; s < numSources ; s++)
	{
		if (!m_srcHaveUncerts[s])
			continue;

		size_t numPoints = m_baseBp->getNumberOfImagePoints(s);
		m_srcImagePositionDifferences[s].resize(numPoints, {0.0f, 0.0f});
		
		float nan = std::numeric_limits<float>::quiet_NaN();
		m_srcAdjustedThetas[s].resize(numPoints, {nan, nan});
		m_srcAdjustedAlphas[s].resize(numPoints, {nan, nan});
		m_srcAdjustedBetas[s].resize(numPoints, {nan, nan});
		m_srcAdjustedPotentials[s].resize(numPoints, nan);

		size_t numImgs = m_baseBp->getNumberOfImages(s);
		assert(numImgs > 0);
		const Vector2Df *pStart = m_baseBp->getThetas(s, 0);
		for (size_t i = 0 ; i < numImgs ; i++)
		{
			const Vector2Df *pVec = m_baseBp->getThetas(s, i);
			size_t offset = pVec - pStart;

			m_srcImageOffsets[s].push_back(offset);

			if (i == numImgs-1)
				assert(offset + m_baseBp->getNumberOfImagePoints(s, i) == numPoints);
		}

		assert(m_srcImageOffsets[s].size() == numImgs);

		numRandomizablePoints += numPoints;
	}

	return true;
}

inline errut::bool_t PositionRandomizationBackprojectWrapper::setRandomOffsets(const std::vector<Vector2Df> &offsets)
{
	// Fill in the data
	size_t off = 0;
	for (size_t s = 0 ; s < m_srcImagePositionDifferences.size() ; s++)
	{
		if (!m_srcHaveUncerts[s])
		{
			assert(m_srcImagePositionDifferences[s].size() == 0);
			continue;
		}

		size_t numEntries = m_srcImagePositionDifferences[s].size();
		if (off + numEntries > offsets.size())
			return "Not enough random positions to fill in required data";

		for (size_t i = 0 ; i < numEntries ; i++)
			m_srcImagePositionDifferences[s][i] = offsets[off++];
	}
	if (off != offsets.size())
		return "Not all provided offsets were used";

	clearCachedValues();
	return true;
}

inline const Vector2D<float> *PositionRandomizationBackprojectWrapper::getThetas(int sourcenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedThetas(sourcenum);
		return m_srcAdjustedThetas[sourcenum].data();
	}
	else
		return m_baseBp->getThetas(sourcenum);
}

inline const Vector2D<float> *PositionRandomizationBackprojectWrapper::getThetas(int sourcenum, int imagenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedThetas(sourcenum);
		assert(sourcenum < m_srcImageOffsets.size() && imagenum >= 0 && imagenum < m_srcImageOffsets[sourcenum].size());
		int offset = m_srcImageOffsets[sourcenum][imagenum];
		assert(offset < m_srcAdjustedThetas[sourcenum].size());
		assert(m_srcAdjustedThetas[sourcenum].size() == m_baseBp->getNumberOfImagePoints(sourcenum));
		return &m_srcAdjustedThetas[sourcenum][offset];
	}
	else
		return m_baseBp->getThetas(sourcenum, imagenum);
}

inline const Vector2D<float> *PositionRandomizationBackprojectWrapper::getAlphas(int sourcenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedAlphas(sourcenum);
		return m_srcAdjustedAlphas[sourcenum].data();
	}
	else
		return m_baseBp->getAlphas(sourcenum);
}

inline const Vector2D<float> *PositionRandomizationBackprojectWrapper::getAlphas(int sourcenum, int imagenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedAlphas(sourcenum);
		assert(sourcenum < m_srcImageOffsets.size() && imagenum >= 0 && imagenum < m_srcImageOffsets[sourcenum].size());
		int offset = m_srcImageOffsets[sourcenum][imagenum];
		assert(offset < m_srcAdjustedAlphas[sourcenum].size());
		assert(m_srcAdjustedAlphas[sourcenum].size() == m_baseBp->getNumberOfImagePoints(sourcenum));
		return &m_srcAdjustedAlphas[sourcenum][offset];
	}
	else
		return m_baseBp->getAlphas(sourcenum, imagenum);
}

inline const Vector2D<float> *PositionRandomizationBackprojectWrapper::getBetas(int sourcenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedBetas(sourcenum);
		return m_srcAdjustedBetas[sourcenum].data();
	}
	else
		return m_baseBp->getBetas(sourcenum);
}

inline const Vector2D<float> *PositionRandomizationBackprojectWrapper::getBetas(int sourcenum, int imagenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedBetas(sourcenum);
		assert(sourcenum < m_srcImageOffsets.size() && imagenum >= 0 && imagenum < m_srcImageOffsets[sourcenum].size());
		int offset = m_srcImageOffsets[sourcenum][imagenum];
		assert(offset < m_srcAdjustedBetas[sourcenum].size());
		assert(m_srcAdjustedBetas[sourcenum].size() == m_baseBp->getNumberOfImagePoints(sourcenum));
		return &m_srcAdjustedBetas[sourcenum][offset];
	}
	else
		return m_baseBp->getBetas(sourcenum, imagenum);
}

inline const float *PositionRandomizationBackprojectWrapper::getLensPotential(int sourcenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedPotentials(sourcenum);
		return m_srcAdjustedPotentials[sourcenum].data();
	}
	else
		return m_baseBp->getLensPotential(sourcenum);
}

inline const float *PositionRandomizationBackprojectWrapper::getLensPotential(int sourcenum, int imagenum) const
{
	if (hasRandomization(sourcenum))
	{
		checkAdjustedPotentials(sourcenum);
		assert(sourcenum < m_srcImageOffsets.size() && imagenum >= 0 && imagenum < m_srcImageOffsets[sourcenum].size());
		int offset = m_srcImageOffsets[sourcenum][imagenum];
		assert(offset < m_srcAdjustedPotentials[sourcenum].size());
		assert(m_srcAdjustedPotentials[sourcenum].size() == m_baseBp->getNumberOfImagePoints(sourcenum));
		return &m_srcAdjustedPotentials[sourcenum][offset];
	}
	else
		return m_baseBp->getLensPotential(sourcenum, imagenum);
}

inline void PositionRandomizationBackprojectWrapper::clearCachedValues()
{
	auto clearEntries = [this](auto &vec)
	{
		assert(vec.size() == m_srcHaveUncerts.size());
		for (auto &x: vec)
			x.clear();
	};

	clearEntries(m_srcAdjustedThetas);
	clearEntries(m_srcAdjustedAlphas);
	clearEntries(m_srcAdjustedBetas);
	clearEntries(m_srcAdjustedPotentials);
}

} // end namespace
