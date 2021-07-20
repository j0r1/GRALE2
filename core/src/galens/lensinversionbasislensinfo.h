#pragma once

#include "graleconfig.h"
#include "gravitationallens.h"
#include "vector2d.h"
#include <memory>

namespace grale
{

class GRALE_IMPORTEXPORT LensInversionBasisLensInfo
{
public:
	LensInversionBasisLensInfo(const std::shared_ptr<GravitationalLens> &lens, Vector2Dd center, double relevantLensingMass);
	LensInversionBasisLensInfo(const LensInversionBasisLensInfo &src);
	LensInversionBasisLensInfo &operator=(const LensInversionBasisLensInfo &src);

	std::shared_ptr<GravitationalLens> m_pLens;
	Vector2Dd m_center;
	double m_relevantLensingMass;
private:
	void copyFrom(const LensInversionBasisLensInfo &src);
};

}
