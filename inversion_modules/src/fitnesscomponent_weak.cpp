#include "fitnesscomponent_weak.h"
#include "fitnesscomponent_overlap.h"
#include "fitnessutil.h"
#include <grale/imagesdataextended.h>
#include <grale/configurationparameters.h>
#include <grale/cosmology.h>
#include <limits>

#include <iostream>

using namespace std;
using namespace errut;

namespace grale
{

// FitnessComponent_WeakLensing

FitnessComponent_WeakLensing::FitnessComponent_WeakLensing(FitnessComponentCache *pCache) 
	: FitnessComponent("weaklensing", pCache)
{
	addRecognizedTypeName("sheardata");
	m_wlType = AveragedEllipticities;
}

FitnessComponent_WeakLensing::~FitnessComponent_WeakLensing()
{
}

bool FitnessComponent_WeakLensing::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "sheardata")
		return true; // ignore

	if (imgDat.getNumberOfImages() != 1)
	{
		setErrorString("Each images data instance can only contain one 'image' for shear info");
		return false;
	}

	if (!imgDat.hasProperty(ImagesData::ShearComponent1) ||
	    !imgDat.hasProperty(ImagesData::ShearComponent2))
	{
		setErrorString("Not all components are present");
		return false;
	}

	if (!imgDat.hasProperty(ImagesData::ShearWeight))
	{
		setErrorString("Shear weights are not present");
		return false;
	}

	needCalcDeflDeriv = true;
	needCalcShear = true;
	needCalcConvergence = true;

	double threshold = 0;
	if (!imgDat.getExtraParameter("threshold", threshold))
	{
		setErrorString("Shear data instance does not contain a (real valued) 'threshold' parameter for |1-kappa|");
		return false;
	}
	if (threshold < 0)
	{
		setErrorString("The threshold value for |1-kappa| must be positive or zero");
		return false;
	}

	m_thresholds.push_back((float)threshold);

	addImagesDataIndex(idx);

	return true;
}

bool FitnessComponent_WeakLensing::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "type")
	{
		string valueStr = value.getStringValue();

		if (valueStr == "AveragedEllipticities")
		{
			cerr << "Setting WL fitness type to AveragedEllipticities" << endl;
			setFitnessType(AveragedEllipticities);
			return true;
		}
		
		if (valueStr == "RealShear")
		{
			cerr << "Setting WL fitness type to RealShear" << endl;
			setFitnessType(RealShear);
			return true;
		}

		if (valueStr == "RealReducedShear")
		{
			cerr << "Setting WL fitness type to RealReducedShear" << endl;
			setFitnessType(RealReducedShear);
			return true;
		}

		setErrorString("Invalid value for '" + optionName + "': '" + valueStr + "'");
		return false;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_WeakLensing::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateWeakLensingFitness(iface, getUsedImagesDataIndices(), m_wlType, m_thresholds);
	return true;
}

// FitnessComponent_WeakLensing_Bayes

FitnessComponent_WeakLensing_Bayes::FitnessComponent_WeakLensing_Bayes(FitnessComponentCache *pCache) 
	: FitnessComponent("bayesweaklensing", pCache)
{
	addRecognizedTypeName("bayesellipticities");
	addRecognizedTypeName("bayesaveragedensityprior");
	addRecognizedTypeName("bayesstronglensing");
	m_redshiftDistributionNeeded = false;
	m_howManySigmaFactor = 3.0f;
	m_numSigmaSamplePoints = 7;
	m_zLens = 0;
	m_maxZ = 0;
	m_zDistSampleCount = 0;
	m_slSigmaArcsec = 2; // TODO: make this configurable!
}

FitnessComponent_WeakLensing_Bayes::~FitnessComponent_WeakLensing_Bayes()
{
}

bool FitnessComponent_WeakLensing_Bayes::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "bayesellipticities" && 
	    typeName != "bayesaveragedensityprior" &&
		typeName != "bayesstronglensing")
		return true; // ignore

	if (typeName == "bayesstronglensing")
	{
		// We're going to reuse some code of the point group rms fitness
		bool_t r = FitnessComponent_PointGroupOverlap::extendedOrPointImageDataToPointGroups(imgDat, m_pointGroups);
		if (!r)
		{
			setErrorString(r.getErrorString());
			return false;
		}

		m_slImages.push_back(idx);
		needCalcDeflections = true;
	}
	else // ellipticities or avg dens prior
	{
		if (imgDat.getNumberOfImages() != 1)
		{
			setErrorString("Each images data instance can only contain one 'image' for ellipticity info or density points");
			return false;
		}

		const int numPoints = imgDat.getNumberOfImagePoints(0);
		if (numPoints == 0)
		{
			setErrorString("No points present");
			return false;
		}

		if (typeName == "bayesellipticities")
		{
			if (!imgDat.hasProperty(ImagesData::ShearComponent1) || !imgDat.hasProperty(ImagesData::ShearComponent2))
			{
				setErrorString("No ellipticity info is present");
				return false;
			}

			if (!imgDat.hasProperty(ImagesData::Redshift) || !imgDat.hasProperty(ImagesData::RedshiftUncertainty))
			{
				setErrorString("The points must have redshift and uncertainty settings (z=0 z_sigma=0 is completely unknown, z=X z_sigma=0 is accurate redshift, z=X z_sigma=Y is redshift with uncertainty)");
				return false;
			}

			// Check Dds ==1 and Ds == 1, so that Dds/Ds == 1
			if (imgDat.getDs() != 1 || imgDat.getDds() != 1)
			{
				setErrorString("For this type of image, both Dds and Ds need to be set to exactly 1");
				return false;
			}

			// Check distance fraction
			m_distanceFractionsForZ.push_back(vector<float>());
			vector<float> &distFracForZ = m_distanceFractionsForZ[m_distanceFractionsForZ.size()-1];

			for (int i = 0 ; i < numPoints ; i++)
			{
				// We don't have the lens redshift here, not the cosmology
				// For now just save the redshift if it needs to be converted
				float z = (float)imgDat.getImagePointProperty(ImagesData::Redshift, 0, i);
				float dz = (float)imgDat.getImagePointProperty(ImagesData::RedshiftUncertainty, 0, i);
				if (z < 0 || dz < 0)
				{
					setErrorString("All redshifts and uncertainties must be non-negative (z=" + to_string(z) + " z_sigma=" + to_string(dz) + ")");
					return false;
				}
				
				float dfToCalculate = -1; // negative signals that it does not need to be calculated
				if (z == 0) // indicates unknown distance fraction, will need probability info
				{
					if (dz != 0)
					{
						setErrorString("For unknown redshifts, the uncertainty must be set to zero (z_sigma=" + to_string(dz) + ")");
						return false;
					}
					m_redshiftDistributionNeeded = true;
				}
				else
				{
					// Keep track of maximum redshift
					// TODO: is 5 sigma a good upper limit to add?
					m_maxZ = std::max(m_maxZ, z + 5.0f*dz);
					if (dz == 0)
						dfToCalculate = z;
				}

				distFracForZ.push_back(dfToCalculate);
			}
			m_elliptImgs.push_back(idx);
		}
		else // avgdensity
		{
			// Nothing else needs to be checked here
		
			m_priorDensImages.push_back(idx);
			needCalcConvergence = true;
		}
	}
	
	needCalcDeflDeriv = true;

	addImagesDataIndex(idx);
	return true;
}

// This is called before inspect
bool FitnessComponent_WeakLensing_Bayes::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "sigmafactor")
	{
		if (value.isArray() || !value.isReal())
		{
			setErrorString("This should be a real number");
			return false;
		}
		m_howManySigmaFactor = value.getRealValue();
		return true;
	}

	if (optionName == "sigmasteps")
	{
		if (value.isArray() || !value.isInteger())
		{
			setErrorString("This should be an integer");
			return false;
		}

		m_numSigmaSamplePoints = value.getIntegerValue();
		if (m_numSigmaSamplePoints < 2)
		{
			setErrorString("At least two points are needed");
			return false;
		}
		return true;
	}

	// TODO: add b/a distribution, check between 0 and 1, calculate area and normalize
	if (optionName == "b_over_a_distribution")
	{
		if (value.isEmpty()) // Allow empty default, check later if it's been set if needed
			return true;

		if (!value.isArray() || !value.isReal() || value.getNumberOfEntries() < 2)
		{
			setErrorString("This should either be an array with at least two real numbers");
			return false;
		}

		// Normalize it ; b/a goes from 0 to 1
		const vector<double> &probdist = value.getRealValues();
		double area = 0;
		for (int i = 0 ; i < probdist.size()-1 ; i++)
		{
			double x0 = (double)i/(double)(probdist.size()-1);
			double x1 = (double)(i+1)/(double)(probdist.size()-1);
			double y0 = probdist[i];
			double y1 = probdist[i+1];
			if (y0 < 0 || y1 < 0)
			{
				setErrorString("All values in b/a distribution must be non-negative");
				return false;
			}

			area += 0.5*(y0+y1)*(x1-x0);
		}
		vector<float> floatProbdist;
		for (auto p : probdist)
			floatProbdist.push_back((float)(p / area));

		m_baDistFunction = make_unique<DiscreteFunction<float>>();
		auto r = m_baDistFunction->init(0.0f, 1.0f, floatProbdist);
		if (!r)
		{
			setErrorString("Unexpected: unable to initialize b/a prob dist function: " + r.getErrorString());
			return false;
		}

		return true;
	}

	if (optionName == "zdist_values")
	{
		if (value.isEmpty()) // Ok, nothing set, but check later if we need it
			return true;

		if (!value.isArray() || !value.isReal() || value.getNumberOfEntries() < 2)
		{
			setErrorString("This should either be empty or an array of at least two real numbers");
			return false;
		}

		// For now we'll just save this, we need range info from another parameter
		for (auto v : value.getRealValues())
			m_zDistValues.push_back((float)v);

		return true;
	}

	if (optionName == "zdist_range")
	{
		if (value.isEmpty()) // Ok, nothing set, but check later if we need it
			return true;

		if (!value.isArray() || !value.isReal() || value.getNumberOfEntries() != 2)
		{
			setErrorString("This should either be empty or an array of two real numbers");
			return false;
		}

		// Just store this for now, initialize it later
		m_zDistMinMax.resize(2);
		m_zDistMinMax[0] = (float)value.getRealValue(0);
		m_zDistMinMax[1] = (float)value.getRealValue(1);

		return true;
	}

	if (optionName == "zdist_numsamples")
	{
		if (value.isArray() || !value.isInteger())
		{
			setErrorString("This should be an integer");
			return false;
		}
		m_zDistSampleCount = value.getIntegerValue();
		if (m_zDistSampleCount < 2)
		{
			setErrorString("At least two samples are needed");
			return false;
		}

		return true;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_WeakLensing_Bayes::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateWeakLensingFitness_Bayes(iface,
		m_pointGroups, m_slImages, m_elliptImgs, m_priorDensImages,
		m_distanceFractionsForZ, m_distFracFunction,
		m_zDistDistFracAndProb,
		*m_baDistFunction.get(),
		m_howManySigmaFactor, m_numSigmaSamplePoints, m_zLens,
		m_slSigmaArcsec);

	return true;
}

bool FitnessComponent_WeakLensing_Bayes::finalize(double zd, const Cosmology *pCosm)
{
	if (pCosm == nullptr)
	{
		setErrorString("Cosmological model needs to be set for this component");
		return false;
	}

	if (zd <= 0) // Should not happen of course, just for safety an extra check
	{
		setErrorString("Lens redshift should be positive");
		return false;
	}

	m_zLens = zd;

	// If a z dist is needed as well, we need to take that maximum into account:
	// The approx distfrac function is used there, so it must be defined for those
	// redshifts as well!
	if (m_redshiftDistributionNeeded)
		for (auto z: m_zDistMinMax)
			m_maxZ = std::max(m_maxZ, z);

	// Determine a discrete distance fraction function
	if (m_maxZ < zd)
		m_maxZ = zd*1.1f; // TODO: should probably not happen anyway
	int distFracPoints = 8192; // TODO: what's a good value here? get this from config option?
	vector<float> distFracs(distFracPoints);
	distFracs[0] = 0; // Angular diameter distance from zd to zd is zero
	for (int i = 1 ; i < distFracPoints ; i++)
	{
		float frac = (float)i/(float)(distFracPoints-1);
		float zs = (1.0f-frac)*zd + frac*m_maxZ;
		distFracs[i] = (float)(pCosm->getAngularDiameterDistance(zd, zs)/pCosm->getAngularDiameterDistance(zs));
	}

	auto r = m_distFracFunction.init(zd, m_maxZ, distFracs);
	if (!r)
	{
		setErrorString("Unexpected: couldn't init discrete distance fraction function: " + r.getErrorString());
		return false;
	}

	// Pre-calculate Dds/Ds for points where it's useful
	for (auto &sources : m_distanceFractionsForZ)
	{
		for (auto &zsToDf : sources)
		{
			// If the value is positive, it can be pre-calculated to a distance fraction
			if (zsToDf > 0)
				zsToDf = m_distFracFunction(zsToDf); // TODO: use exact calculation here instead?
		}
	}

	//cerr << "DEBUG: m_redshiftDistributionNeeded = " << (int)m_redshiftDistributionNeeded << endl;
	if (m_redshiftDistributionNeeded)
	{
		if (m_zDistValues.size() == 0)
		{
			setErrorString("Not all redshifts are known, a general redshift distribution is needed");
			return false;
		}
		if (m_zDistMinMax.size() != 2)
		{
			setErrorString("Start/end of the redshift distribution range is needed");
			return false;
		}

		m_zDistFunction = make_unique<DiscreteFunction<float>>();
		auto r = m_zDistFunction->init(m_zDistMinMax[0], m_zDistMinMax[1], m_zDistValues);
		if (!r)
		{
			setErrorString("Unable to initialize the redshift distribution function: " + r.getErrorString());
			return false;
		}

		float x0 = zd;
		float x1 = m_zDistMinMax[1];
		// cerr << "DEBUG: m_zDistMinMax = " << m_zDistMinMax[0] << " " << m_zDistMinMax[1] << endl;

		// The distance ratio for zd will be zero, we don't need to consider this
		// If the probability for the highest z is also zero, we don't need to
		// consider that either
		if ((*m_zDistFunction)(x1) == 0)
		{
			float diff = (x1-x0)/(float)(m_zDistSampleCount+1);
			x0 += diff;
			x1 -= diff;
		}
		else
		{
			float diff = (x1-x0)/(float)m_zDistSampleCount;
			x0 += diff;
		}

		//cerr << "DEBUG: x0 = " << x0 << " x1 = " << x1 << endl;
		float l0, l1;
		m_distFracFunction.getLimits(l0, l1);
		//cerr << "DEBUG: l0 = " << l0 << " l1 = " << l1 << endl;

		//cerr << "DISTFRACS" << endl;
		m_zDistDistFracAndProb.clear();
		for (int i = 0 ; i < m_zDistSampleCount ; i++)
		{
			float frac = (float)i/(float)(m_zDistSampleCount-1);
			float z = x0*(1.0f-frac) + x1*frac;
			//cerr << "DEBUG: z = " << z << endl;

			if (z > l1 || z < l0)
			{
				setErrorString("Internal error: z out of bounds of approximate distfrac function");
				return false;
			}

			float zProb = (*m_zDistFunction)(z);
			float distFrac = m_distFracFunction(z);

			m_zDistDistFracAndProb.push_back({ distFrac, zProb });
			//cerr << distFrac << "," << zProb << "," << endl;
		}
	}

	if (m_baDistFunction.get() == nullptr)
	{
		setErrorString("No b/a distribution has been set");
		return false;
	}

	return true;
}

} // end namespace
