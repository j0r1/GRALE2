#pragma once

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include <memory>

namespace grale
{

class PositionRandomizationBackprojectWrapper : public ProjectedImagesInterface
{
public:
	PositionRandomizationBackprojectWrapper(const std::shared_ptr<ProjectedImagesInterface> &baseBackProjector) : m_baseBp(baseBackProjector) { }
	~PositionRandomizationBackprojectWrapper() { }

	double getLensDistance() const override { return m_baseBp->getLensDistance(); }
	double getLensRedshift() const override { return m_baseBp->getLensRedshift(); }
	double getAngularScale() const override { return m_baseBp->getAngularScale(); };

	// TODO: These should be modified
	void checkHaveRandomization(int sourceNum) const
	{
		// TODO!
	}

	void checkHaveNoRandomization(int sourceNum) const
	{
		// TODO!
	}

	const Vector2D<float> *getBetas(int sourcenum) const override
	{
		checkHaveRandomization(sourcenum);
		return m_baseBp->getBetas(sourcenum);
	}
	const Vector2D<float> *getBetas(int sourcenum, int imagenum) const override
	{
		checkHaveRandomization(sourcenum);
		return m_baseBp->getBetas(sourcenum, imagenum);
	}
	const Vector2D<float> *getAlphas(int sourcenum) const override
	{
		checkHaveRandomization(sourcenum);
		return m_baseBp->getAlphas(sourcenum);
	}
	const Vector2D<float> *getAlphas(int sourcenum, int imagenum) const override
	{
		checkHaveRandomization(sourcenum);
		return m_baseBp->getAlphas(sourcenum, imagenum);
	}
	const Vector2D<float> *getThetas(int sourcenum) const override
	{
		checkHaveRandomization(sourcenum);
		return m_baseBp->getThetas(sourcenum);
	}
	const Vector2D<float> *getThetas(int sourcenum, int imagenum) const override
	{
		checkHaveRandomization(sourcenum);
		return m_baseBp->getThetas(sourcenum, imagenum);
	}
	const float *getLensPotential(int sourceNumber) const override
	{
		checkHaveRandomization(sourceNumber);
		return m_baseBp->getLensPotential(sourceNumber);
	}
	const float *getLensPotential(int sourceNumber, int imageNumber) const override
	{
		checkHaveRandomization(sourceNumber);
		return m_baseBp->getLensPotential(sourceNumber, imageNumber);
	}

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
		assert(0);
	}
private:
	std::shared_ptr<ProjectedImagesInterface> m_baseBp;
};

} // end namespace
