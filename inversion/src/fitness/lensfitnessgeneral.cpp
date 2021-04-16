/*

  This file is a part of the GRALE Lens inversion modules, which determine
  the fitness of a lens model for the genetic algorithm used by the
  GRALE library.

  Copyright (C) 2008-2012 Jori Liesenborgs

  Contact: jori.liesenborgs@gmail.com
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  
*/

#include "lensfitnessgeneral.h"
#include "fitnesscomponent_overlap.h"
#include "fitnesscomponent_null.h"
#include "fitnesscomponent_weak.h"
#include "fitnesscomponent_time.h"
#include "fitnesscomponent_caustic.h"
#include "fitnesscomponent_dens.h"
#include "fitnesscomponent_defl.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include <limits>
#include <algorithm>
#include <sstream>
#include <memory>

using namespace std;

namespace grale
{

LensFitnessGeneral::LensFitnessGeneral()
{
	m_initialized = false;
	m_pShortComponent = 0;
	m_deleteShortComponent = true;
	m_numFitnessComponents = 0;
	m_pCache = 0;
}

LensFitnessGeneral::~LensFitnessGeneral()
{
	clear();
}

void LensFitnessGeneral::clear()
{
	if (m_deleteShortComponent)
		delete m_pShortComponent;
	m_pShortComponent = 0;
	m_deleteShortComponent = true;

	for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
		delete m_totalComponents[i];
	m_totalComponents.clear();

	m_totalInverse = false;
	m_totalShear = false; 
	m_totalConvergence = false;
	m_totalDeflectionFlags.clear();
	m_totalDerivativeFlags.clear();
	m_totalPotentialFlags.clear();
	m_totalInverseFlags.clear();
	m_totalShearFlags.clear();
	m_totalConvergenceFlags.clear();

	m_shortInverse = false;
	m_shortShear = false; 
	m_shortConvergence = false;
	m_shortDeflectionFlags.clear();
	m_shortDerivativeFlags.clear();
	m_shortPotentialFlags.clear();
	m_shortInverseFlags.clear();
	m_shortShearFlags.clear();
	m_shortConvergenceFlags.clear();

	delete m_pCache;
	m_pCache = 0;

	m_cosmology = nullptr;
}

inline bool reduceFlags(const vector<bool> &flags)
{
	for (size_t i = 0 ; i < flags.size() ; i++)
	{
		if (flags[i])
			return true;
	}
	return false;
}

void removeEmpty(vector<FitnessComponent *> &comp)
{
	vector<FitnessComponent *> newComp;

	for (size_t i = 0 ; i < comp.size() ; i++)
	{
		FitnessComponent *pComp = comp[i];
		if (pComp)
			newComp.push_back(pComp);
	}

	comp = newComp;
}

#define COMPONENT_POINTOVERLAP_IDX				0
#define COMPONENT_EXTENDEDOVERLAP_IDX			1
#define COMPONENT_POINTGROUPOVERLAP_IDX			2
#define COMPONENT_WEAK_IDX						3
#define COMPONENT_POINTNULL_IDX					4
#define COMPONENT_EXTENDEDNULL_IDX				5
#define COMPONENT_TIMEDELAY_IDX					6
#define COMPONENT_KAPPATHRESHOLD_IDX			7
#define COMPONENT_CAUSTICPENALTY_IDX			8
#define COMPONENT_KAPPAGRADIENT_IDX			    9
#define COMPONENT_BAYESWEAK_IDX					10
#define COMPONENT_DEFLECTIONANGLE_IDX			11
#define COMPONENT_IDX_MAX						12

static const vector<string> componentNames {
	"pointimageoverlap",
	"extendedimageoverlap",
	"pointgroupoverlap",
	"weaklensing",
	"pointimagenull",
	"extendedimagenull",
	"timedelay",
	"kappathreshold",
	"causticpenalty",
	"kappagradient",
	"bayesweaklensing",
	"deflectionangle"
};

vector<FitnessComponent*> getAllComponents(FitnessComponentCache *pCache)
{
	return {
		new FitnessComponent_PointImagesOverlap(pCache),
		new FitnessComponent_ExtendedImagesOverlap(pCache),
		new FitnessComponent_PointGroupOverlap(pCache),
		new FitnessComponent_WeakLensing(pCache),
		new FitnessComponent_NullSpacePointImages(pCache),
		new FitnessComponent_NullSpaceExtendedImages(pCache),
		new FitnessComponent_TimeDelay(pCache),
		new FitnessComponent_KappaThreshold(pCache),
		new FitnessComponent_CausticPenalty(pCache),
		new FitnessComponent_KappaGradient(pCache),
		new FitnessComponent_WeakLensing_Bayes(pCache),
		new FitnessComponent_DeflectionAngle(pCache)
	};
}

FitnessComponent *LensFitnessGeneral::totalToShort(const vector<FitnessComponent *> &total, vector<int> &shortImageIndices,
							   const ConfigurationParameters &params)
{
	assert(componentNames.size() == COMPONENT_IDX_MAX);
	assert(componentNames.size() == total.size());

	map<int, vector<FitnessComponent *>> priorityMap;

	// Build a map of the priorities
	for (int compIdx = 0 ; compIdx < componentNames.size() ; compIdx++)
	{
		string compName = componentNames[compIdx];
		FitnessComponent *pComp = total[compIdx];
		string keyName = "scalepriority_" + compName;

		int priority = 0;
		if (!params.getParameter(keyName, priority))
		{
			setErrorString("Can't find (integer) parameter '" + keyName + "': " + params.getErrorString());
			return nullptr;
		}

		if (!pComp)
			continue;

		if (pComp->getObjectName() != compName)
		{
			setErrorString("Internal error: component name '" + pComp->getObjectName() + "' is not the expected '" + compName + "'");
			return nullptr;
		}
	
		if (priority < 0) // Skip this
			continue;

		priorityMap[priority].push_back(pComp);
	}

	// Sort according to lowest priority
	vector<int> priorityKeys;
	for (auto it = priorityMap.begin() ; it != priorityMap.end() ; ++it)
		priorityKeys.push_back(it->first);
	sort(priorityKeys.begin(), priorityKeys.end());

	if (priorityKeys.size() == 0)
	{
		setErrorString("No suitable component could be found to use mass scale");
		return nullptr;
	}

	auto components = priorityMap[priorityKeys[0]]; // Use the components with the lowest priority
	sort(components.begin(), components.end(), [](FitnessComponent *c1, FitnessComponent *c2)
	{
		if (c1 == nullptr)
			return false;
		if (c2 == nullptr)
			return true;
		return c1->getNumberOfUsedImages() > c2->getNumberOfUsedImages();
	});

	assert(components.size() > 0 && components[0] != nullptr);
	
	FitnessComponent *pShortComp = components[0];

	shortImageIndices = pShortComp->getUsedImagesDataIndices();
	return pShortComp->createShortCopy();
}

bool componentSortFunction(const FitnessComponent *pComp1, const FitnessComponent *pComp2)
{
	assert(pComp1);
	assert(pComp2);

	if (pComp1->getPriority() < pComp2->getPriority())
		return true;
	if (pComp1->getPriority() == pComp2->getPriority())
	{
		// For the same priority, use the number of images involved (that can actually be
		// used as a constraint) as a tie breaker
		if (pComp1->getNumberOfUsedImages() > pComp2->getNumberOfUsedImages())
			return true;
	}
	return false;
}

ConfigurationParameters *LensFitnessGeneral::getDefaultParametersInstance() const
{
	ConfigurationParameters *pParams = new ConfigurationParameters();

	pParams->setParameterEmpty("general_cosmology");

	pParams->setParameter("priority_pointimageoverlap", 300);
	pParams->setParameter("priority_extendedimageoverlap", 300);
	pParams->setParameter("priority_pointgroupoverlap", 250);
	pParams->setParameter("priority_pointimagenull", 200);
	pParams->setParameter("priority_extendedimagenull", 200);
	pParams->setParameter("priority_weaklensing", 500);
	pParams->setParameter("priority_timedelay", 400);
	pParams->setParameter("priority_kappathreshold", 600);
	pParams->setParameter("priority_causticpenalty", 100);
	pParams->setParameter("priority_kappagradient", 700);
	pParams->setParameter("priority_bayesweaklensing", 500);
	pParams->setParameter("priority_deflectionangle", 800);

	pParams->setParameter("scalepriority_pointimageoverlap", 100);
	pParams->setParameter("scalepriority_extendedimageoverlap", 100);
	pParams->setParameter("scalepriority_pointgroupoverlap", 200);
	pParams->setParameter("scalepriority_pointimagenull", -1);
	pParams->setParameter("scalepriority_extendedimagenull", -1);
	pParams->setParameter("scalepriority_weaklensing", 300);
	pParams->setParameter("scalepriority_timedelay", -1);
	pParams->setParameter("scalepriority_kappathreshold", -1);
	pParams->setParameter("scalepriority_causticpenalty", -1);
	pParams->setParameter("scalepriority_kappagradient", 1000);
	pParams->setParameter("scalepriority_bayesweaklensing", 300);
	pParams->setParameter("scalepriority_deflectionangle", 1100);

	pParams->setParameter("fitness_pointimageoverlap_scaletype", string("MinMax"));
	pParams->setParameter("fitness_pointgroupoverlap_rmstype", string("AllBetas"));
	pParams->setParameter("fitness_timedelay_type", string("NoSrc"));
	pParams->setParameter("fitness_weaklensing_type", string("AveragedEllipticities"));

	pParams->setParameterEmpty("fitness_bayesweaklensing_zdist_values");
	pParams->setParameterEmpty("fitness_bayesweaklensing_zdist_range");
	pParams->setParameter("fitness_bayesweaklensing_zdist_numsamples", 16);
	pParams->setParameterEmpty("fitness_bayesweaklensing_b_over_a_distribution");
	pParams->setParameter("fitness_bayesweaklensing_sigmafactor", 3.0);
	pParams->setParameter("fitness_bayesweaklensing_sigmasteps", 7);
	pParams->setParameter("fitness_bayesweaklensing_stronglenssigma", 0.0);
	return pParams;
}

bool LensFitnessGeneral::setFitnessOptions(FitnessComponent *pComp, const ConfigurationParameters *pParams)
{
	vector<string> allKeys;
	pParams->getAllKeys(allKeys);

	string fitnessOptionStart = "fitness_" + pComp->getObjectName() + "_";
	for (auto &k : allKeys)
	{
		if (k.find(fitnessOptionStart) == 0) // key starts with this
		{
			string optionName = k.substr(fitnessOptionStart.length());
			TypedParameter tp;

			pParams->getParameter(k, tp);
			if (!pComp->processFitnessOption(optionName, tp))
			{
				setErrorString("Unable to process fitness option '" + k + "': " + pComp->getErrorString());
				return false;
			}
		}
	}
	return true;
};

bool LensFitnessGeneral::processGeneralParameters(const ConfigurationParameters *pParams)
{
	TypedParameter tp;

	// Check cosmology
	if (!pParams->getParameter("general_cosmology", tp))
	{
		setErrorString("No 'general_cosmology' parameter present");
		return false;
	}
	if (tp.isEmpty())
	{
		m_cosmology = nullptr; // allow this, perhaps it's not needed
		cerr << "DEBUG: cleared cosmology" << endl;
	}
	else
	{
		if (!(tp.isReal() && tp.isArray() && tp.getNumberOfEntries() == 5))
		{
			setErrorString("Cosmology setting in 'general_cosmology' must be an array of five real values (h, Omega_m, Omega_r, Omega_v, w)");
			return false;
		}
		const vector<double> &v = tp.getRealValues();
		m_cosmology = make_shared<Cosmology>(v[0], v[1], v[2], v[3], v[4]);
		cerr << "DEBUG: created cosmology with parameters ("
		     << m_cosmology->getH() << ","
			 << m_cosmology->getOmegaM() << ","
			 << m_cosmology->getOmegaR() << ","
			 << m_cosmology->getOmegaV() << ","
			 << m_cosmology->getW() << ")" << endl;
	}
	
	return true;
}

bool LensFitnessGeneral::processComponentParameters(const ConfigurationParameters *pParams)
{
	// Set priorities and component specific fitness options
	for (FitnessComponent *pComp : m_totalComponents)
	{
		assert(pComp);
		string keyName = "priority_" + pComp->getObjectName();

		int priority = 0;
		if (!pParams->getParameter(keyName, priority))
		{
			setErrorString("Can't find (integer) parameter '" + keyName + "': " + pParams->getErrorString());
			return false;
		}

		pComp->setPriority(priority);
		if (!setFitnessOptions(pComp, pParams))
			return false;
	}

	return true;
}

set<string> LensFitnessGeneral::getSupportedTypeNames()
{
	set<string> supportedNames;
	for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
	{
		FitnessComponent *pComp = m_totalComponents[i];
		assert(pComp);

		vector<string> names = pComp->getRecognizedTypeNames();
		for (size_t n = 0 ; n < names.size() ; n++)
		{
			supportedNames.insert(names[n]);
			//cerr << "SUPPORTED NAME: " << names[n] << endl;
		}
	}
	return supportedNames;
}

bool LensFitnessGeneral::init(double z_d, std::list<ImagesDataExtended *> &images, 
                              std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
	if (m_initialized)
	{
		setErrorString("Already initialized");
		return false;
	}

	if (pParams == 0)
	{
		setErrorString("No parameters with priorities of the different fitness measures were specified");
		return false;
	}

	clear();

	// Using unique_ptr make sure it gets deleted on return (in case of an error)
	unique_ptr<FitnessComponentCache> cacheSmart(new FitnessComponentCache(images.size()));
	FitnessComponentCache *pCache = cacheSmart.get();

	// Clear the config markers (to check if there are unused config options)
	for (auto pImg : images)
		pImg->clearRetrievalMarkers();
	
	// Populate m_totalComponents
	m_totalComponents = getAllComponents(pCache);
	assert(m_totalComponents.size() == COMPONENT_IDX_MAX);

	if (!processGeneralParameters(pParams))
		return false;

	if (!processComponentParameters(pParams))
		return false;

	if (!checkImagesDataParameters(images))
		return false;

	// Let the components process the images data
	if (!inspectImagesByComponents(images, m_totalComponents,
			m_totalDeflectionFlags,
			m_totalDerivativeFlags,
			m_totalPotentialFlags,
			m_totalInverseFlags,
			m_totalShearFlags,
			m_totalConvergenceFlags,
			m_totalInverse,
			m_totalShear,
			m_totalConvergence))
		return false;

	// Clear the components that are not used, finalize others
	if (!finalizeAndCleanUnusedComponents(z_d, m_cosmology.get()))
		return false;

	vector<int> shortIndices;
	// This function creates a copy, but further on we'll check if we can actually
	// reuse an already existing component, and discard the copy
	m_pShortComponent = totalToShort(m_totalComponents, shortIndices, *pParams);
	// Check that we have something to base the mass scale on
	if (m_pShortComponent == nullptr)
		return false; // error string has been set

	// TODO: find something better to report
	cerr << "DBG: shortComponent = " << m_pShortComponent->getObjectName() << endl;

	// Also set the possibly different fitness options
	setFitnessOptions(m_pShortComponent, pParams);

	buildShortImagesList(images, shortIndices, shortImages);

	// Run the shortImages through the newly created m_pShortComponent
	vector<FitnessComponent *> shortComponentVector { m_pShortComponent };
	if (!inspectImagesByComponents(shortImages, shortComponentVector,
			m_shortDeflectionFlags,
			m_shortDerivativeFlags,
			m_shortPotentialFlags,
			m_shortInverseFlags,
			m_shortShearFlags,
			m_shortConvergenceFlags,
			m_shortInverse,
			m_shortShear,
			m_shortConvergence))
		return false;

	//cerr << "Short finalize for " << m_pShortComponent->getObjectName() << endl;
	if (!m_pShortComponent->finalize(z_d, m_cosmology.get()))
	{
		setErrorString("Unable to finalize (short) component '" + m_pShortComponent->getObjectName() + "':" + m_pShortComponent->getErrorString());
		return false;
	}

	// Clear unused components from the list, and determine the number of fitness components
	removeEmpty(m_totalComponents);
	m_numFitnessComponents = m_totalComponents.size();

	if (m_numFitnessComponents == 0)
	{
		setErrorString("No fitness components could be used");
		return false;
	}

	if (m_numFitnessComponents == 1) 
	{
		// We can use the same 'short' version as the total one
		// (we've already checked that a short version can be created, so
		// this can now only be the remaining fitness components)
		if (m_deleteShortComponent)
			delete m_pShortComponent;

		assert(m_totalComponents.size() == 1);
		m_pShortComponent = m_totalComponents[0];
		assert(m_pShortComponent);

		m_deleteShortComponent = false; // don't delete things twice!
		
		shortImages.clear(); // Make sure the caller knows we're not using a separate list
	}

	// Save the order things should be calculated (for possible caching
	// dependencies) and order the components according to priority
	// This may be necessary for null space for example: there the estimated
	// source shape is calculated by the extended overlap component
	m_calculationOrderComponents = m_totalComponents;

	sort(m_totalComponents.begin(), m_totalComponents.end(), componentSortFunction);
	for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
		m_totalComponents[i]->setPriorityOrder(i);

	if (!checkAllImagesUsed(images))
		return false;

	// Check for parameter names that are specified but not used
	if (!checkUnusedImagesDataParameters(images))
		return false;

	// TODO: should find something better than stdout
	cerr << "FITNESS COMPONENT ORDER: ";
	for (const FitnessComponent *pComp : m_totalComponents)
	{
		assert(pComp);
		cerr << pComp->getObjectName() << " ";
	}
	cerr << endl;

	stringstream ss;
	ss << m_totalComponents[0]->getObjectName();
	for (size_t i = 1 ; i < m_totalComponents.size() ; i++)
		ss << " " << m_totalComponents[i]->getObjectName();

	m_fitnessComponentDescription = ss.str();

	// We don't want the cache to be deleted, extract it from the smart pointer
	m_pCache = cacheSmart.release();

	m_initialized = true;

	return true;
}

bool LensFitnessGeneral::checkImagesDataParameters(list<ImagesDataExtended *> &images)
{
	// Get the supported type names
	set<string> supportedNames = getSupportedTypeNames();

	// Let the components process the images data
	int imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		string numStr = to_string(imgDatCount+1);

		ImagesDataExtended *pImgDat = *it;
		assert(pImgDat);
	
		// Note that Dds is set to zero for multi-plane case
		if (pImgDat->getDds() < 0 || pImgDat->getDs() <= 0)
		{
			setErrorString("Source/lens and source/observer distances must be positive for images data set " + numStr);
			return false;
		}
		if (pImgDat->getDds() == 0) // Probably multi-plane case, we need a 'z'
		{
			if (!pImgDat->hasExtraParameter("z"))
			{
				setErrorString("Dds is set to zero for images data " + numStr + ", which indicates a multi-plane scenario, but 'z' is not specified");
				return false;
			}

			double z;
			if (!pImgDat->getExtraParameter("z", z) || z <= 0)
			{
				setErrorString("Redshift parameter 'z' is present for images data set " + numStr + ", but is not positive, or of wrong type");
				return false;
			}
		}

		// Check the type name
		string typeName;
		if (!pImgDat->getExtraParameter("type", typeName))
		{
			setErrorString("Can't find a 'type' string in images data set " + numStr + ": " + pImgDat->getErrorString());
			return false;
		}

		if (supportedNames.find(typeName) == supportedNames.end())
		{
			setErrorString("Images data set " + numStr + " type name '" + typeName + "' is not supported");
			return false;
		}
	}
	return true;
}

bool LensFitnessGeneral::inspectImagesByComponents(list<ImagesDataExtended *> &images,
	vector<FitnessComponent *> &components,
	vector<bool> &deflectionFlags,
	vector<bool> &derivativeFlags,
	vector<bool> &potentialFlags,
	vector<bool> &inverseFlags,
	vector<bool> &shearFlags,
	vector<bool> &convergenceFlags,
	bool &calcInverse, bool &calcShear, bool &calcConvergence)
{
	deflectionFlags.clear();
	derivativeFlags.clear();
	potentialFlags.clear();
	inverseFlags.clear();
	shearFlags.clear();
	convergenceFlags.clear();

	// Let the components process the images data
	int imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		string numStr = to_string(imgDatCount+1);

		ImagesDataExtended *pImgDat = *it;
		assert(pImgDat);
	
		// Process image in this component (component should ignore image if unknown type name)

		bool totNeedDefl = false, totNeedDeflDeriv = false, totNeedPotential = false,
			 totNeedInvMag = false, totNeedShear = false, totNeedConv = false;

		for (auto pComp : components)
		{
			assert(pComp);

			bool needDefl = false, needDeflDeriv = false, needPotential = false,
				 needInvMag = false, needShear = false, needConv = false;

			if (!pComp->inspectImagesData(imgDatCount, *pImgDat, needDefl, needDeflDeriv, needPotential,
						                  needInvMag, needShear, needConv))
			{
				setErrorString("Unable to process images data entry " + numStr + " by component " +
						       pComp->getObjectName() + ": " + pComp->getErrorString());
				return false;
			}
			
			// Sanity check
			if ((needInvMag || needShear || needConv) && !needDeflDeriv)
			{
				setErrorString("Internal error processing imagesdata " + numStr + " by component " +
							   pComp->getObjectName() + ": didn't set need for deflection derivatives, "
							   "but inverse magnifications, shear or convergence was requested");
				return false;
			}

			if (needDefl) totNeedDefl = true;
			if (needDeflDeriv) totNeedDeflDeriv = true;
			if (needPotential) totNeedPotential = true;
			if (needInvMag) totNeedInvMag = true;
			if (needShear) totNeedShear = true;
			if (needConv) totNeedConv = true;
		}

		deflectionFlags.push_back(totNeedDefl);
		derivativeFlags.push_back(totNeedDeflDeriv);
		potentialFlags.push_back(totNeedPotential);
		inverseFlags.push_back(totNeedInvMag);
		shearFlags.push_back(totNeedShear);
		convergenceFlags.push_back(totNeedConv);
	}

	calcInverse = reduceFlags(inverseFlags);
	calcShear = reduceFlags(shearFlags);
	calcConvergence = reduceFlags(convergenceFlags);

	return true;
}

bool LensFitnessGeneral::finalizeAndCleanUnusedComponents(float z_d, const Cosmology *pCosmology)
{
	for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
	{
		FitnessComponent *pComp = m_totalComponents[i];
		assert(pComp);

		if (pComp->getUsedImagesDataIndices().size() == 0)
		{
			delete pComp;
			m_totalComponents[i] = 0;
		}
		else
		{
			//cerr << "Long finalize for " << pComp->getObjectName() << endl;
			if (!pComp->finalize(z_d, pCosmology))
			{
				setErrorString("Unable to finalize component '" + pComp->getObjectName() + "':" + pComp->getErrorString());
				return false;
			}
		}
	}
	return true;
}

bool LensFitnessGeneral::checkUnusedImagesDataParameters(list<ImagesDataExtended *> &images)
{
	int imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		string numStr = to_string(imgDatCount+1);

		ImagesDataExtended *pImgDat = *it;
		assert(pImgDat);

		vector<string> unretrievedKeys;
		pImgDat->getUnretrievedKeys(unretrievedKeys);

		if (unretrievedKeys.size() > 0)
		{
			stringstream ss;

			ss << "Unused keys in images data set " + numStr + ":";
			for (auto &k : unretrievedKeys)
				ss << " " << k;

			setErrorString(ss.str());
			return false;
		}
	}
	return true;
}

bool LensFitnessGeneral::checkAllImagesUsed(list<ImagesDataExtended *> &images)
{
	// Check that each images data set is used
	vector<bool> imageCheckFlags(images.size(), false);
	for (auto pComp : m_totalComponents)
	{
		assert(pComp);
		vector<int> indices = pComp->getUsedImagesDataIndices();
		for (auto idx : indices)
		{
			assert(idx >= 0 && idx < imageCheckFlags.size());
			imageCheckFlags[idx] = true;
		}
	}

	for (size_t i = 0 ; i < imageCheckFlags.size() ; i++)
	{
		if (!imageCheckFlags[i])
		{
			setErrorString("Images data set " + to_string(i+1) + " isn't used by any fitness component");
			return false;
		}
	}
	return true;
}

void LensFitnessGeneral::buildShortImagesList(list<ImagesDataExtended *> &images, 
	const vector<int> &shortIndices,
	list<ImagesDataExtended *> &shortImages)
{
	shortImages.clear();

	// Store the relevant image indices in a set
	set<int> shortIndicesSet(shortIndices.begin(), shortIndices.end());

	// Store the relevant images data instances in shortImages
	int imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		if (shortIndicesSet.find(imgDatCount) != shortIndicesSet.end()) // This is a useful image
			shortImages.push_back(*it);
	}

	assert(shortImages.size() > 0);
}

bool LensFitnessGeneral::calculateMassScaleFitness(const ProjectedImagesInterface &iface, float &fitness) const
{
	assert(m_pShortComponent);

	// it's possible that this cache is used by the short component
	// (it may be a pointer to one of the m_totalComponents), so
	// clear the cache
	m_pCache->clear(); 

	float f = numeric_limits<float>::max();
	if (!m_pShortComponent->calculateFitness(iface, f))
	{
		setErrorString("Unable to calculate mass scale fitness in '" + m_pShortComponent->getObjectName()  + "': " + m_pShortComponent->getErrorString());
		return false;
	}
	fitness = f;
	return true;
}

bool LensFitnessGeneral::calculateOverallFitness(const ProjectedImagesInterface &iface, float *pFitnessValues) const
{
	assert(m_numFitnessComponents == m_totalComponents.size());
	assert(m_numFitnessComponents == m_calculationOrderComponents.size());

	// clear the cache
	m_pCache->clear(); 

	for (int i = 0 ; i < m_numFitnessComponents ; i++)
	{
		FitnessComponent *pComp = m_calculationOrderComponents[i];
		assert(pComp);

		float f = numeric_limits<float>::max();

		if (!pComp->calculateFitness(iface, f))
		{
			setErrorString("Unable to calculate fitness in '" + pComp->getObjectName()  + "': " + pComp->getErrorString());
			return false;
		}

		int idx = pComp->getPriorityOrder();
		assert(idx >= 0 && idx < m_numFitnessComponents);
		pFitnessValues[idx] = f;
	}

	return true;
}

string LensFitnessGeneral::getUsage() const
{
	return "The usage documentation for this module has been migrated to https://research.edm.uhasselt.be/jori/grale2/usage_general.html";
}

} // end namespace

