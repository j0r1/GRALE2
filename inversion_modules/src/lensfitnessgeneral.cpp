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
#include "fitnesscomponent.h"
#include <grale/imagesdataextended.h>
#include <grale/configurationparameters.h>
#include <limits>
#include <set>
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
	m_totalStoreIntens = false;
	m_totalStoreTimeDelay = false;
	m_totalStoreShear = false;
	m_totalDeflectionFlags.clear();
	m_totalDerivativeFlags.clear();
	m_totalPotentialFlags.clear();
	m_totalInverseFlags.clear();
	m_totalShearFlags.clear();
	m_totalConvergenceFlags.clear();

	m_shortInverse = false;
	m_shortShear = false; 
	m_shortConvergence = false;
	m_shortStoreIntens = false;
	m_shortStoreTimeDelay = false;
	m_shortStoreShear = false;
	m_shortDeflectionFlags.clear();
	m_shortDerivativeFlags.clear();
	m_shortPotentialFlags.clear();
	m_shortInverseFlags.clear();
	m_shortShearFlags.clear();
	m_shortConvergenceFlags.clear();

	delete m_pCache;
	m_pCache = 0;
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
#define COMPONENT_IDX_MAX						9

FitnessComponent *totalToShort(const vector<FitnessComponent *> &total, vector<int> &shortImageIndices)
{
	assert(total.size() == COMPONENT_IDX_MAX);

	FitnessComponent *pPointOverlap = total[COMPONENT_POINTOVERLAP_IDX];
	assert(!pPointOverlap || pPointOverlap->getObjectName() == "pointimageoverlap");

	FitnessComponent *pExtOverlap = total[COMPONENT_EXTENDEDOVERLAP_IDX];
	assert(!pExtOverlap || pExtOverlap->getObjectName() == "extendedimageoverlap");

	FitnessComponent *pPtGrpOverlap = total[COMPONENT_POINTGROUPOVERLAP_IDX];
	assert(!pPtGrpOverlap || pPtGrpOverlap->getObjectName() == "pointgroupoverlap");

	FitnessComponent *pWeak = total[COMPONENT_WEAK_IDX];
	assert(!pWeak || pWeak->getObjectName() == "weaklensing");

	// More than 1 image needs to be present in these to make them useful
	// (should also be checked in the components, but checking again for
	// safety)
	if (pPointOverlap && pPointOverlap->getNumberOfUsedImages() < 2)
		pPointOverlap = 0;
	if (pExtOverlap && pExtOverlap->getNumberOfUsedImages() < 2)
		pExtOverlap = 0;
	if (pPtGrpOverlap && pPtGrpOverlap->getNumberOfUsedImages() < 2)
		pPtGrpOverlap = 0;

	vector<FitnessComponent*> components { pPointOverlap, pExtOverlap, pPtGrpOverlap };
	sort(components.begin(), components.end(), [](FitnessComponent *c1, FitnessComponent *c2)
	{
		if (c1 == nullptr)
			return false;
		if (c2 == nullptr)
			return true;
		return c1->getNumberOfUsedImages() > c2->getNumberOfUsedImages();
	});

	if (components[0] != nullptr)
	{
		shortImageIndices = components[0]->getUsedImagesDataIndices();
		if (components[0] == pPointOverlap)
			return new FitnessComponent_PointImagesOverlap(nullptr);
		if (components[0] == pExtOverlap)
			return new FitnessComponent_ExtendedImagesOverlap(nullptr);

		assert(components[0] == pPtGrpOverlap);
		return new FitnessComponent_PointGroupOverlap(nullptr);
	}

	// As a final resort, use the weak lensing data to base the mass scale on
	if (pWeak)
	{
		shortImageIndices = pWeak->getUsedImagesDataIndices();
		return new FitnessComponent_WeakLensing(nullptr);
	}

	// Nothing useful could be found
	return nullptr;
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

inline string itos(int i) { return to_string(i); }

ConfigurationParameters *LensFitnessGeneral::getDefaultParametersInstance() const
{
	ConfigurationParameters *pParams = new ConfigurationParameters();

	pParams->setParameter("priority_pointimageoverlap", 300);
	pParams->setParameter("priority_extendedimageoverlap", 300);
	pParams->setParameter("priority_pointgroupoverlap", 250);
	pParams->setParameter("priority_pointimagenull", 200);
	pParams->setParameter("priority_extendedimagenull", 200);
	pParams->setParameter("priority_weaklensing", 500);
	pParams->setParameter("priority_timedelay", 400);
	pParams->setParameter("priority_kappathreshold", 600);
	pParams->setParameter("priority_causticpenalty", 100);

	pParams->setParameter("fitness_pointimageoverlap_scaletype", string("MinMax"));
	pParams->setParameter("fitness_pointgroupoverlap_rmstype", string("AllBetas"));
	pParams->setParameter("fitness_timedelay_type", string("Paper2009"));
	pParams->setParameter("fitness_timedelay_relative", false);
	pParams->setParameter("fitness_weaklensing_type", string("AveragedEllipticities"));
	return pParams;
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
	FitnessComponent_PointImagesOverlap *pPtComponent = new FitnessComponent_PointImagesOverlap(pCache);
	FitnessComponent_PointGroupOverlap *pPtGrpComponent = new FitnessComponent_PointGroupOverlap(pCache);
	FitnessComponent_TimeDelay *pTDComponent = new FitnessComponent_TimeDelay(pCache);
	FitnessComponent_WeakLensing *pWLComponent = new FitnessComponent_WeakLensing(pCache);

	/*vector<FitnessComponent *>*/ m_totalComponents = {
		pPtComponent,
		new FitnessComponent_ExtendedImagesOverlap(pCache),
		pPtGrpComponent,
		pWLComponent,
		new FitnessComponent_NullSpacePointImages(pCache),
		new FitnessComponent_NullSpaceExtendedImages(pCache),
		pTDComponent,
		new FitnessComponent_KappaThreshold(pCache),
		new FitnessComponent_CausticPenalty(pCache)
	};
	assert(m_totalComponents.size() == COMPONENT_IDX_MAX);

	vector<string> allKeys;
	pParams->getAllKeys(allKeys);

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
	}

	// Get the supported type names
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

	vector<bool> storeOrigIntensFlags, storeOrigTimeDelayFlags, storeOrigShearFlags;
	vector<bool> storeShortOrigIntensFlags, storeShortOrigTimeDelayFlags, storeShortOrigShearFlags;

	// Let the components process the images data
	int imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		string numStr(itos(imgDatCount+1));

		ImagesDataExtended *pImgDat = *it;
		assert(pImgDat);
	
		double Dds = pImgDat->getDds();
		double Ds = pImgDat->getDs();

		if (Dds <= 0 || Ds <= 0)
		{
			setErrorString("Source/lens and source/observer distances must be positive for images data set " + numStr);
			return false;
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

		// Process image in this component (component should ignore image if unknown type name)

		bool totNeedDefl = false, totNeedDeflDeriv = false, totNeedPotential = false,
			 totNeedInvMag = false, totNeedShear = false, totNeedConv = false,
			 totStoreOrigIntens = false, totStoreOrigTimeDelay = false, totStoreOrigShear = false;

		for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
		{
			FitnessComponent *pComp = m_totalComponents[i];
			assert(pComp);

			bool needDefl = false, needDeflDeriv = false, needPotential = false,
				 needInvMag = false, needShear = false, needConv = false,
				 storeOrigIntens = false, storeOrigTimeDelay = false, storeOrigShear = false;

			if (!pComp->inspectImagesData(imgDatCount, *pImgDat, needDefl, needDeflDeriv, needPotential,
						                  needInvMag, needShear, needConv, storeOrigIntens,
										  storeOrigTimeDelay, storeOrigShear))
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
			if (storeOrigIntens) totStoreOrigIntens = true;
			if (storeOrigTimeDelay) totStoreOrigTimeDelay = true;
			if (storeOrigShear) totStoreOrigShear = true;
		}

		m_totalDeflectionFlags.push_back(totNeedDefl);
		m_totalDerivativeFlags.push_back(totNeedDeflDeriv);
		m_totalPotentialFlags.push_back(totNeedPotential);
		m_totalInverseFlags.push_back(totNeedInvMag);
		m_totalShearFlags.push_back(totNeedShear);
		m_totalConvergenceFlags.push_back(totNeedConv);

		storeOrigIntensFlags.push_back(totStoreOrigIntens);
		storeOrigTimeDelayFlags.push_back(totStoreOrigTimeDelay);
		storeOrigShearFlags.push_back(totStoreOrigShear);
	}

	m_totalInverse = reduceFlags(m_totalInverseFlags);
	m_totalShear = reduceFlags(m_totalShearFlags);
	m_totalConvergence = reduceFlags(m_totalConvergenceFlags);
	m_totalStoreIntens = reduceFlags(storeOrigIntensFlags);
	m_totalStoreTimeDelay = reduceFlags(storeOrigTimeDelayFlags);
	m_totalStoreShear = reduceFlags(storeOrigShearFlags);

	// Clear the components that are not used, finalize others
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
			if (!pComp->finalize())
			{
				setErrorString("Unable to finalize component '" + pComp->getObjectName() + "':" + pComp->getErrorString());
				return false;
			}
		}
	}

	vector<int> shortIndices;
	set<int> shortIndicesSet;
	m_pShortComponent = totalToShort(m_totalComponents, shortIndices);
	
	// Store the relevant image indices in a set
	for (auto idx : shortIndices)
		shortIndicesSet.insert(idx);

	// Store the relevant images data instances in shortImages
	imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		if (shortIndicesSet.find(imgDatCount) != shortIndicesSet.end()) // This is a useful image
			shortImages.push_back(*it);
	}
	
	// Check that we have something to base the mass scale on
	if (m_pShortComponent == 0)
	{
		setErrorString("No images data type was present that we can base the mass scale on");
		return false;
	}

	// Run the shortImages through the newly created m_pShortComponent
	assert(m_pShortComponent);
	assert(shortImages.size() > 0);

	imgDatCount = 0;
	for (auto it = shortImages.begin() ; it != shortImages.end() ; ++it, imgDatCount++)
	{
		string numStr(itos(imgDatCount+1));

		ImagesDataExtended *pImgDat = *it;
		assert(pImgDat);

		bool needDefl = false, needDeflDeriv = false, needPotential = false,
			 needInvMag = false, needShear = false, needConv = false,
			 storeOrigIntens = false, storeOrigTimeDelay = false, storeOrigShear = false;

		if (!m_pShortComponent->inspectImagesData(imgDatCount, *pImgDat, needDefl, needDeflDeriv, needPotential,
									  needInvMag, needShear, needConv, storeOrigIntens,
									  storeOrigTimeDelay, storeOrigShear))
		{
			setErrorString("Unexpected: unable to process 'short' images data entry " + numStr + " by component " +
						   m_pShortComponent->getObjectName() + ": " + m_pShortComponent->getErrorString());
			return false;
		}
		
		// Sanity check
		if ((needInvMag || needShear || needConv) && !needDeflDeriv)
		{
			setErrorString("Unexpected: internal error processing 'short' images data " + numStr + " by component " +
						   m_pShortComponent->getObjectName() + ": didn't set need for deflection derivatives, "
						   "but inverse magnifications, shear or convergence was requested");
			return false;
		}

		m_shortDeflectionFlags.push_back(needDefl);
		m_shortDerivativeFlags.push_back(needDeflDeriv);
		m_shortPotentialFlags.push_back(needPotential);
		m_shortInverseFlags.push_back(needInvMag);
		m_shortShearFlags.push_back(needShear);
		m_shortConvergenceFlags.push_back(needConv);

		storeShortOrigIntensFlags.push_back(storeOrigIntens);
		storeShortOrigTimeDelayFlags.push_back(storeOrigTimeDelay);
		storeShortOrigShearFlags.push_back(storeOrigShear);
	}
	m_shortInverse = reduceFlags(m_shortInverseFlags);
	m_shortShear = reduceFlags(m_shortShearFlags);
	m_shortConvergence = reduceFlags(m_shortConvergenceFlags);
	m_shortStoreIntens = reduceFlags(storeShortOrigIntensFlags);
	m_shortStoreTimeDelay = reduceFlags(storeShortOrigTimeDelayFlags);
	m_shortStoreShear = reduceFlags(storeShortOrigShearFlags);

	//cerr << "Short finalize for " << m_pShortComponent->getObjectName() << endl;
	if (!m_pShortComponent->finalize())
	{
		setErrorString("Unable to finalize (short) component '" + m_pShortComponent->getObjectName() + "':" + m_pShortComponent->getErrorString());
		return false;
	}

	// Time delays are requested, make sure that z_d is set!
	if (m_totalComponents[COMPONENT_TIMEDELAY_IDX] != 0)
	{
		if (z_d <= 0)
		{
			setErrorString("Time delays are begin used, but lens redshift is not positive");
			return false;
		}
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
	m_calculationOrderComponents = m_totalComponents;

	sort(m_totalComponents.begin(), m_totalComponents.end(), componentSortFunction);
	for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
		m_totalComponents[i]->setPriorityOrder(i);

	// Check that each images data set is used
	vector<bool> imageCheckFlags(images.size(), false);
	for (size_t i = 0 ; i < m_totalComponents.size() ; i++)
	{
		FitnessComponent *pComp = m_totalComponents[i];
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
			setErrorString("Images data set " + itos(i+1) + " isn't used by any fitness component");
			return false;
		}
	}

	// Check for parameter names that are specified but not used
	imgDatCount = 0;
	for (auto it = images.begin() ; it != images.end() ; ++it, imgDatCount++)
	{
		string numStr(itos(imgDatCount+1));

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
| priority_kappathreshold             | 600                     |
| fitness_pointgroupoverlap_rmstype   | 'AllBetas'              |
| fitness_pointimageoverlap_scaletype | 'MinMax'                |
| fitness_timedelay_type              | 'Paper2009'             |
| fitness_timedelay_relative          | `False`                 |
| fitness_weaklensing_type            | 'AveragedEllipticities' |

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

The `fitness_timedelay_type` can also be 'ExperimentalI' or 'ExperimentalII',
but as you might guess, these are still experimental. The option
`fitness_timedelay_relative` is extremely experimental: don't enable this, it
does not work yet!

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
}

} // end namespace

