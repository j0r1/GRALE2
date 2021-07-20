#include "lensinversionbasislensinfo.h"

namespace grale
{

LensInversionBasisLensInfo::LensInversionBasisLensInfo(const std::shared_ptr<GravitationalLens> &lens, Vector2Dd center, double relevantLensingMass)
	: m_pLens(lens), m_center(center), m_relevantLensingMass(relevantLensingMass)
{
}

LensInversionBasisLensInfo::LensInversionBasisLensInfo(const LensInversionBasisLensInfo &src)
{
	copyFrom(src);
}

LensInversionBasisLensInfo &LensInversionBasisLensInfo::operator=(const LensInversionBasisLensInfo &src)
{
	copyFrom(src);
	return *this;
}

void LensInversionBasisLensInfo::copyFrom(const LensInversionBasisLensInfo &src)
{
	m_pLens = src.m_pLens;
	m_center = src.m_center;
	m_relevantLensingMass = src.m_relevantLensingMass;
}

} // end namespace
