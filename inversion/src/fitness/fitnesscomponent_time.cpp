#include "fitnesscomponent_time.h"
#include "fitnessutil.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include <sstream>

#include <iostream>

using namespace std;

namespace grale
{

// FitnessComponent_TimeDelay

FitnessComponent_TimeDelay::FitnessComponent_TimeDelay(const std::shared_ptr<FitnessComponentCache> &pCache) 
	: FitnessComponent("timedelay", pCache),
	  m_fitnessType(NoSrc)
{
	addRecognizedTypeName("pointimages");
	addRecognizedTypeName("extendedimages");
	addRecognizedTypeName("pointgroupimages");
}

FitnessComponent_TimeDelay::~FitnessComponent_TimeDelay()
{
}

bool FitnessComponent_TimeDelay::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "extendedimages" && typeName != "pointimages")
		return true; // not relevant, ignore

	if (!imgDat.hasTimeDelays())
	{
		bool useTd = false;
		if (imgDat.hasExtraParameter("timedelay"))
		{
			if (!imgDat.getExtraParameter("timedelay", useTd))
			{
				setErrorString("Extra parameter 'timedelay' was specified, but is not a boolean");
				return false;
			}
			if (useTd)
			{
				setErrorString("Extra parameter 'timedelay' is set, but no time delays are present");
				return false;
			}
		}
		return true; // no time delays present, and usage was not requested. ignore
	}
	
	// Time delays are present in image, by default they should be used
	if (imgDat.hasExtraParameter("timedelay"))
	{
		bool useTd;

		if (!imgDat.getExtraParameter("timedelay", useTd))
		{
			setErrorString("Extra parameter 'timedelay' was specified, but is not a boolean");
			return false;
		}
		if (!useTd) // time delay info was explicitly specified to be skipped
			return true; // just ignore this
	}

	// Check that there are no negative time delay values: in the past this
	// was used to ignore a specific value, we'll generate an error if such
	// a value is still present
	int numTimeDelays = imgDat.getNumberOfTimeDelays();
	for (int tIdx = 0 ; tIdx < numTimeDelays ; tIdx++)
	{
		int img, point;
		double delay;
		imgDat.getTimeDelay(tIdx, &img, &point, &delay);
		if (delay < 0)
		{
			stringstream ss;
			ss << "A negative time delay (" << delay << ") was present for image " << img << ", point " << point;
			setErrorString(ss.str());
			return false;
		}
	}

	float tdScaleFactor = 0;
	if (imgDat.hasExtraParameter("timedelay_scalefactor"))
	{
		double v;

		if (!imgDat.getExtraParameter("timedelay_scalefactor", v))
		{
			setErrorString("Extra parameter 'timedelay_scalefactor' was specified, but is not a floating point value");
			return false;
		}
		if (v <= 0)
		{
			setErrorString("Extra parameter 'timedelay_scalefactor' should be positive");
			return false;
		}
		tdScaleFactor = (float)v;
	}
	
	needCalcPotential = true;

	//cerr << "Added TD fitness for " << idx << endl;
	addImagesDataIndex(idx);
	m_tdScaleFactors.push_back(tdScaleFactor);

	return true;
}

bool FitnessComponent_TimeDelay::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "type")
	{
		string valueStr = value.getStringValue();
		if (valueStr == "Paper2009")
		{
			cerr << "Regular Time Delay Fitness (from 2009 article)" << endl;
			setFitnessType(Paper2009);
			return true;
		}
		
		if (valueStr == "NoSrc")
		{
			cerr << "EXPERIMENTAL TIME DELAY FITNESS II - NoSrc" << endl;
			setFitnessType(NoSrc);
			return true;
		}

		setErrorString("Invalid value for '" + optionName + "': '" + valueStr + "'");
		return false;
	}

	if (optionName == "nosrc_cutoff")
	{
		if (!value.isReal())
		{
			setErrorString("Value for '" + optionName + "' should be a real number"); 
			return false;
		}

		double v = value.getRealValue();
		if (v < 0)
		{
			setErrorString("Value for '" + optionName + "' should be zero or positive, but is " + to_string(v));
			return false;
		}

		m_nosrcCutoff = v;
		return true;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_TimeDelay::finalize(double zd, const Cosmology *pCosm)
{
	if (zd <= 0)
	{
		setErrorString("Lens redshift is not positive");
		return false;
	}

	if (m_nosrcCutoff > 0 && m_fitnessType != NoSrc)
	{
		setErrorString("Fitness cutoff was set to " + to_string(m_nosrcCutoff) + ", but this is only allowed for the NoSrc calculation type");
		return false;
	}

	return true;
}

bool FitnessComponent_TimeDelay::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	if (m_fitnessType == Paper2009)
		fitness = calculateTimeDelayFitnessPaper2009(iface, getUsedImagesDataIndices(), m_tdScaleFactors);
	else if (m_fitnessType == NoSrc)
	{
		fitness = calculateTimeDelayFitnessNoSrc(iface, getUsedImagesDataIndices(), m_tdScaleFactors);
		if (fitness < m_nosrcCutoff)
			fitness = 0;
	}
	else
	{
		setErrorString("Unknown TD fitness type");
		return false;
	}
	return true;
}

} // end namespace
