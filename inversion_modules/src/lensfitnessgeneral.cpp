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
#include <grale/imagesdataextended.h>
#include <grale/configurationparameters.h>
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
#define COMPONENT_IDX_MAX						11

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
	"bayesweaklensing"
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
		new FitnessComponent_WeakLensing_Bayes(pCache)
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

	pParams->setParameter("scalepriority_pointimageoverlap", 100);
	pParams->setParameter("scalepriority_extendedimageoverlap", 100);
	pParams->setParameter("scalepriority_pointgroupoverlap", 200);
	pParams->setParameter("scalepriority_pointimagenull", -1);
	pParams->setParameter("scalepriority_extendedimagenull", -1);
	pParams->setParameter("scalepriority_weaklensing", 300);
	pParams->setParameter("scalepriority_timedelay", -1);
	pParams->setParameter("scalepriority_kappathreshold", -1);
	pParams->setParameter("scalepriority_causticpenalty", -1);
	pParams->setParameter("scalepriority_kappagradient", -1);
	pParams->setParameter("scalepriority_bayesweaklensing", 300);

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
#ifndef WIN32
	return R"TEXT(
Inversion module 'general'
==========================

This lens inversion module can handle a variety of input data: image data,
null space grids, time delay info, etc. Based on the precise input that was
provided, one or more fitness measures will be used and the genetic
algorithm will try to optimize these. The currently implemented fitness 
measures are:

 - `extendedimageoverlap`: for multiply imaged extended sources, this fitness
   measure calculates the fractional overlap of back-projected images as 
   described in [Liesenborgs et al (2006)][1]. Point sources can be included
   as well, in that case the average size of the extended images is used
   as the distance scale for determining how well the back-projected points
   overlap.

 - `pointimageoverlap`: for strongly lensed point sources, this is a measure
   for how well the back-projected point images overlap in the source plane,
   as described in [Zitrin et al (2010)][2]. By default the size of the
   area encompassing all backprojected points is used as a distance scale
   in determining the overlap, but the median of absolute deviations
   [MAD](https://en.wikipedia.org/wiki/Median_absolute_deviation) can be
   used as well (see the section about *'module parameters'*).

 - `pointgroupoverlap`: an experimental fitness measure, which backprojects
   specific points and estimates the RMS in the image plane. Either the
   average source position is used, or each backprojected image is used
   separately as a source position estimate (see the section about *'module parameters'*).
   The magnification matrix is used to map differences in the source plane to 
   differences in the image plane.

 - `extendedimagenull`: this fitness measure takes the null space into 
   account, the region where no images should be found, as described in
   [Liesenborgs et al (2007)][3]. An extension was made so allow it to add
   a penalty when multiple images are generated when only a single image
   should be allowed.

 - `pointimagenull`: similar to the null space fitness measure for extended
   images but now for point images, also described in [Zitrin et al (2010)][2].
   Here too the extension was made to allow singly imaged sources to be used
   as a contraint.

 - `timedelay`: when time delay information is available for one or more
   multiply-imaged sources, a fitness measure will be calculated as described
   in [Liesenborgs et al (2009)][4]. Other, still experimental time delay
   fitness measures can be used as well by changing the parameter
   `fitness_timedelay_type` (see the section about *'module parameters'*).

 - `causticpenalty`: also described in [Liesenborgs et al (2009)][4], it is
   also possible to add a fitness measure to try to prevent a caustic from
   intersecting a reconstructed source.

 - `weaklensing`: if shear information is provided, a $\chi^2$-like fitness 
   measure will be provided for this. By default, the shear data will be
   interpreted as averaged ellipticities, but if desired (for testing
   purposes for example), real shear, or actual reduced shear can be
   specified as well (see the section about *'module parameters'*). For
   each data point, a weight can be specified.

 - `kappathreshold`: if it's certain that the convergence ($\kappa$) values
   at certain locations should not exceed a certain threshold, this fitness
   measure is added.

 - `bayesweaklensing`: TODO, bayesian, per galaxy

In general, when more than one fitness measure is used, at the end of the 
algorithm a set of mass maps is found in which no single one will be better
than another with respect to all fitness measures. While you can specify that
all of these solutions need to be saved to files, usually you'll just want
to keep a single solution. To choose this solution, the order of importance 
of the fitness measures can be left to the default or you can specify it 
yourself. See the section about *'module parameters'* for more information.

Input images data list
----------------------

In [`InversionWorkSpace.addImageDataToList](http://research.edm.uhasselt.be/~jori/grale2/grale_inversion.html#grale.inversion.InversionWorkSpace.addImageDataToList),
you specify the image type, as well as (optionally) extra parameters.
The images data set type name describes what kind of data
is contained in this entry, and can be one of the following:

 - `pointimages`: the images data set describes point images that originate
   from a single source (also a point of course). Alternatively, extended
   images may be provided as well, but in this case point groups specifying
   which points are really images of the same source plane point must be
   present. In general in a strong
   lensing scenario there will be more than one image, but single images are
   allowed as well. If this is the case, the image can't be used directly 
   because it originates from a single point source by definition, but it can
   still be useful for a null space constraint to inform the algorithm that
   this source should not lie in a region that produces multiple images.

    For the overlap fitness, two options are possible: if the same `groupname`
   is used for more multiple point images, the overlap of the points will be
   measured relative to the scale of only these back-projected points. By
   default, the scale of all backprojected points is used. If `usescalefrom`
   is present, then for those points the scale set by another set of points,
   identified by a specific `groupname`. 

    In case time delay information was added to the images data set, by 
   default a time delay fitness measure will be calculated as well. This can
   still be disabled by specifying a boolean extra parameter `timedelay` that
   is set to `False`.

    The fitness measures to which this type of data is relevant, are 
   `pointimageoverlap`, `pointimagenull` and `timedelay`.

 - `pointgroupimages`: an image data set with type is used in the
   `pointgroupoverlap` fitness measure, which calculates the RMS in the
   image plane. The input should either be point image data, or extended
   images in which point groups are used to indicate which points belong
   together.

    In case time delay information was added to the images data set, by 
   default a time delay fitness measure will be calculated as well. This can
   still be disabled by specifying a boolean extra parameter `timedelay` that
   is set to `False`.

    The fitness measures to which this type of data is relevant, are
   `pointgroupoverlap`, `pointimagenull` and `timedelay`.

 - `pointnullgrid`: the images data set should contain only one 'image', a
   triangulated grid, and is interpreted as a null space grid. It refers to 
   the previous images data set of `pointimages` type and the algorithm will 
   use this grid to estimate the number of extra images that are generated by 
   a trial solution. 
   
    One extra parameter can be specified: `weight`, which should be a real and 
   positive value. The estimated amount of images is multiplied by this
   number, and can be used to adjust the relative importance of some null 
   spaces. This can be useful to provide extra weight to singly imaged sources
   because they may otherwise provide only a small penalty.

    The only fitness measure to which this type of data is relevant, is 
   `pointimagenull`.

 - `extendedimages`: in this case the images data set describes the extended 
   images originating from a single source. Therefore, in principle each image
   should consist of at least three points. A set of point images may be
   provided as well, but since each separate point image does not provide
   a distance scale to measure overlap with, the average size of the other
   backprojected (extended images) is used as the distance scale instead.
   A single image
   containing exactly one point is allowed as well. In this case it will not
   be used as a constraint regarding the overlap of back-projected images
   in the source plane, but it can be used as a constraint in the null space,
   indicating that the corresponding source plane position should produce only
   one image.
 
    For extended images, the extra parameters `userectangles` and 
   `usepointgroups` (both taking a boolean value) can be useful. By default,
   the overlap of surrounding back-projected rectangles is always used, but
   can be disabled with the first option. If point groups (points in different
   images that correspond to each other) are available in the images data set,
   they are used by default as well. The second option can disable their use.

    In case time delay information was added to the images data set, by 
   default a time delay fitness measure will be calculated as well. This can
   still be disabled by specifying a boolean extra parameter `timedelay` that
   is set to `no`.

    The fitness measures to which this type of data is relevant, are 
   `extendedimageoverlap`, `extendedimagenull`, `pointimagenull` and `timedelay`.

 - `extendednullgrid`: similar to `pointnullgrid`, the images data set should 
   contain only one 'image', a triangulated grid, and is interpreted as a null 
   space grid. Here, typically the regions containing observed images or 
   regions that may harbor an as yet undetected image are removed from the 
   grid. The images data set refers to the previous images data set of 
   `extendedimages` type and the algorithm will use this grid to calculate the
   sizes of the regions in the image plane that contain additional images.
   
    One extra parameter can be specified: `weight`, which should be a real and 
   positive value. The total size of the additional images is multiplied by 
   this number, and can be used to adjust the relative importance of some null
   spaces. This can be useful to provide extra weight to singly imaged sources
   because they may otherwise provide only a small penalty.

    The only fitness measure to which this type of data is relevant, is 
   `extendedimagenull`.

 - `sheardata`: use this type to provide shear measurements to the
   algorithm. The images data set should contain only one 'image', a set of 
   points in the image plane for which shear components have been
   specified. One extra parameter (a real number) called `threshold` must be
   provided: this contains a threshold value for $|1-\kappa|$. Only when at
   a certain point the value for $|1-\kappa|$ exceeds the specified threshold,
   will it be included in the $\chi^2$-like calculation.

    The only fitness measure to which this type of data is relevant, is 
   `weaklensing`.

 - `bayesellipticities`: TODO, galaxy ellipticities, for `bayesweaklensing`

 - `kappathresholdpoints`: this data set should only contain a single 'image',
   a set of points at which the convergence $\kappa$ is calculated. The
   mandatory extra (real valued) parameter `threshold` specifies if a penalty
   is needed for a certain point: if its convergence is below the threshold,
   no penalty is added to the fitness measure, otherwise the amount by which
   the threshold is exceeded is added.

    The only fitness measure to which this type of data is relevant, is 
   `kappathreshold`.

 - `causticgrid`: when an images data set of this type is specified, it should
   contain a triangulated grid, which will be used to estimate the caustics.
   The algorithm will calculate the length of the caustics that intersect the 
   estimated source, and this source estimate is based on the previously
   encountered images data set marked as `extendedimages`.
 
    The only fitness measure to which this type of data is relevant, is 
   `causticpenalty`.

 - `singlyimagedpoints`: the images data set of this type should contain only
   one 'image' entry, which may consist of several points. Each point is assumed
   to have only a single image, i.e. it should not lie in a multiply imaged
   region. These data are therefore intended to be used together with a
   null space grid of type `pointnullgrid`. The points could also be specified 
   separately as `pointimages`, but this way is more convenient if several 
   points are needed for which the same null space grid can be used.

    The only fitness measure to which this type of data is relevant, is 
   `pointimagenull`.

Module parameters
-----------------

Extra parameters for this module can be set using the `fitnessObjectParameters`
argument in e.g. [`inversion.invert`](http://research.edm.uhasselt.be/~jori/grale2/grale_inversion.html#grale.inversion.InversionWorkSpace.invert).
The defaults can be obtained using the command
[`inversion.getDefaultModuleParameters`](http://research.edm.uhasselt.be/~jori/grale2/grale_inversion.html#grale.inversion.getDefaultModuleParameters), 
and are listed in the following table:

| Parameter name                      | Value                   |
|-------------------------------------|-------------------------|
| priority_causticpenalty             | 100                     |
| priority_pointimagenull             | 200                     |
| priority_extendedimagenull          | 200                     |
| priority_pointgroupoverlap          | 250                     |
| priority_pointimageoverlap          | 300                     |
| priority_extendedimageoverlap       | 300                     |
| priority_timedelay                  | 400                     |
| priority_weaklensing                | 500                     |
| priority_bayesweaklensing           | 500                     |
| priority_kappathreshold             | 600                     |
| scalepriority_pointimageoverlap     | 100                     |
| scalepriority_extendedimageoverlap  | 100                     |
| scalepriority_pointgroupoverlap     | 200                     |
| scalepriority_pointimagenull        | -1                      |
| scalepriority_extendedimagenull     | -1                      |
| scalepriority_weaklensing           | 300                     |
| scalepriority_timedelay             | -1                      |
| scalepriority_kappathreshold        | -1                      |
| scalepriority_causticpenalty        | -1                      |
| scalepriority_kappagradient         | -1                      |
| scalepriority_bayesweaklensing      | 300                     |
| fitness_pointgroupoverlap_rmstype   | 'AllBetas'              |
| fitness_pointimageoverlap_scaletype | 'MinMax'                |
| fitness_timedelay_type              | 'NoSrc'                 |
| fitness_weaklensing_type            | 'AveragedEllipticities' |
| fitness_bayesweaklensing_zdist_values          | None         |
| fitness_bayesweaklensing_zdist_range           | None         |
| fitness_bayesweaklensing_zdist_numsamples      | 16           |
| fitness_bayesweaklensing_b_over_a_distribution | None         |
| fitness_bayesweaklensing_sigmafactor           | 3.0          |
| fitness_bayesweaklensing_sigmasteps            | 7            |

In case input is provided with 'pointgroupimages' type, the 'pointgroupoverlap'
fitness calculation is used, which estimates the RMS in the image plane. It
does this by projecting the image points onto the source plane, determining
a source position based on these points, and using the magnification matrix
to convert differences in the source plane to difference in the image plane.
By default, each backprojected image point is used as a possible source
position, corresponding to the value 'AllBetas' of `fitness_pointgroupoverlap_rmstype`.
To use the averaged source position instead, you can set it to 'AverageBeta'.

If input is provided of 'pointimages' type, the 'pointimageoverlap' fitness
calculation will be activated. It projects the image points onto the source
plane, and uses the differences between the backprojected points to base the
fitness measure on. The distance scale with which these differences are measured
depends on all backprojected points and by default the size of the entire area
is used. This corresponds to the setting `fitness_pointimageoverlap_scaletype` 
to 'MinMax'. In case the [median of absolute deviations](https://en.wikipedia.org/wiki/Median_absolute_deviation)
should be used instead, this option can be set to 'MAD'.

The `fitness_timedelay_type` can also be 'Paper2009', to use the older one from
the [2009 article](https://ui.adsabs.harvard.edu/abs/2009MNRAS.397..341L/abstract).

For the weak lensing fitness, the data is by default interpreted as (averaged)
ellipticity measurements. This corresponds the default value of 'AveragedEllipticities'
for `fitness_weaklensing_type`. In case true shear is supplied, or the actual
reduced shear, the value can also be 'RealShear' or 'RealReducedShear' respectively.
This is mainly meant for testing purposes.

As the names suggest, the options that start with `priority_` describe priorities 
for fitness measures. These values do **not** have any effect on the way the 
genetic algorithm operates, only on the final solution that is chosen from the 
set of 'best' solutions.

Here, the lower the priority value, the more important it is considered to be,
so if for example both the `extendedimageoverlap` and `extendedimagenull`
fitness measures are used based on the provided images data sets, there may
be more than one 'best' solution found by the genetic algorithm. You may
decide to save all these solutions by setting `returnNds` to `True` in
the `invert` function, but typically you'll want to go on using a single solution. 
	
It is based on the priorities that were specified, that one solution will be 
chosen. Using the default settings in our example, first the solution(s) will 
be chosen with the lowest value for the `extendedimagenull` fitness value,
and then for `extendedimageoverlap`. You can force this order to be turned
around by overriding some of these values in the `fitnessObjectParameters`
argument.

As you can see, some priorities have the same value. In case two fitness
measures with the same priority are used, the number of images is typically
used as a tie breaker. For example, if you've specified both point images
and extended images as input, the fitness measures `pointimageoverlap` and
`extendedimageoverlap` will be used. When choosing a single final solution,
the fitness value that corresponds to most images will have the best
priority for the default settings. So if there are more point images than
extended images, the point image criterion will be considered first.

[1]: http://adsabs.harvard.edu/abs/2006MNRAS.367.1209L "A genetic algorithm 
for the non-parametric inversion of strong lensing systems"
[2]: http://adsabs.harvard.edu/abs/2010MNRAS.408.1916Z "Full lensing analysis 
of Abell 1703: comparison of independent lens-modelling techniques"
[3]: http://adsabs.harvard.edu/abs/2007MNRAS.380.1729L "Non-parametric 
inversion of gravitational lensing systems with few images using a 
multi-objective genetic algorithm"
[4]: http://adsabs.harvard.edu/abs/2009MNRAS.397..341L "Non-parametric strong 
lens inversion of SDSS J1004+4112"

)TEXT";
#else
	return "(Stub for Win32 - compiler says 'string too big')";
#endif // !WIN32
}

} // end namespace

