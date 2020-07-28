/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

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

/**
 * \file lensfitnessobject.h
 */

#ifndef GRALE_LENSFITNESSOBJECT_H

#define GRALE_LENSFITNESSOBJECT_H

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include <errut/errorbase.h>
#include <vector>
#include <list>

namespace grale
{
	
class ImagesDataExtended;
class ConfigurationParameters;

/** This class is a base class for fitness calculation code in the GA. */
class GRALE_IMPORTEXPORT LensFitnessObject : public errut::ErrorBase 
{
protected:
	LensFitnessObject()								{ }
public:
	virtual ~LensFitnessObject() 							{ }
	
	/** Initialize this object.
	 *  Initialize this object.
	 *  \param z_d The redshift of the lens.
	 *  \param images When calling this function, this list of images data basically
	 *                contains all the images data you specified as input, for example
	 *                using the GRALESHELL command \c 'imgdata/list/add'. You can modify
	 *                this list, and for all the images remaining in the list, the
	 *                points will be projected back onto the corresponding source plane
	 *                (or will be used to calculate other properties like shear or convergence).
	 *                So if you need images data as input for your module, but don't
	 *                actually need the points to be projected onto the source plane
	 *                (or need other properties that depend on a specific lens model),
	 *                you can remove these images data instances from the list.
	 *  \param massScaleImages As explained in \ref inversionmodule "writing inversion modules",
	 *                the genetic algorithm first determines a suitable mass scale by optimizing
	 *                the value returned by LensFitnessObject::calculateMassScaleFitness. When
	 *                this \c 'init' function is called, the \c massScaleImages list is empty, and you can
	 *                store in this list the entries from the \c images list which are needed in
	 *                the mass scale fitness calculation. If it remains empty, the same list as
	 *                \c images will be used in the mass scale fitness calculation.
	 */
	virtual bool init(double z_d, std::list<ImagesDataExtended *> &images, 
	                  std::list<ImagesDataExtended *> &massScaleImages,
					  const ConfigurationParameters *pParams) = 0;

	/** In case the lens fitness object uses additional parameters for initialization,
	 *  an instance containing the defaults can be created this way. */
	virtual ConfigurationParameters *getDefaultParametersInstance() const { return nullptr; }

	/** Obtain a usage description for this LensFitnessObject implementation. */
	virtual std::string getUsage() const = 0;

	/** Get information about the fitness values (e.g. "extendedimages extendednullspace") */
	virtual std::string getFitnessComponentsDescription() const = 0;

	/** Given the situation described by the ProjectedImagesInterface instance, this
	 *  function must calculate a single fitness value (lower is better) which will
	 *  be used to find an appropriate scale of the basis functions in the GA. */
	virtual bool calculateMassScaleFitness(const ProjectedImagesInterface &interface0, float &fitness) const = 0;

	/** This function should inspect the situation described by the ProjectedImagesInterface
	 *  instance, and calculate all corresponding fitness values. */ 
	virtual bool calculateOverallFitness(const ProjectedImagesInterface &interface0, float *pFitnessValues) const = 0;

	/** The number of fitness components your lens fitness object will calculate
	 *  (is 1 for a single-objective GA, 2 or more for a multi-objective GA). */
	virtual int getNumberOfFitnessComponents() const = 0;

	/** The length of the images data list in the LensFitnessObject::init function
	 *  must be a multiple of this number.
	 *  For more complex optimizations, typically several images data instances correspond
	 *  to a single source. For example an images data instance that corresponds to the actual
	 *  image points of a source could be followed by an instance specifying the null space grid.
	 *  In that case this function could return 2 and a check will be done that the list (after
	 *  calling LensFitnessObject::init) is a multiple of 2. It is only used for such a check
	 *  though, so you can also just return 1. */
	virtual int getImagesGroupSize() const = 0;
	
	// These are used in the factory initialization, to indicate what should be calculated and what original
	// information should be store
	/** @name Basic calculation flags
	 *
	 *  The following functions provide information to the BackProjectMatrixNew instances in the genetic
	 *  algorithm about what information should be calculated and stored.
	 **/
	///@{
	/** The LensFitnessObject instance should create these arrays of flags to indicate what
	 *  properties should be calculated for which source in the original (possibly modified)
	 *  images data list.
	 *  The LensFitnessObject instance should create these arrays of flags to indicate what
	 *  properties should be calculated for which source in the original (possibly modified)
	 *  images data list. You should resize each vector to match the size of the images
	 *  vector in the LensFitnessObject::init function, and for each source store what
	 *  properties should be calculated.
	 *  \param deflectionFlags For which sources should the deflection be calculated (you'll
	 *                         need this to access the ProjectedImagesInterface::getBetas
	 *                         functions from within the genetic algorithm)
	 *  \param derivativeFlags Flags to indicate for which sources the derivatives of the
	 *                         deflection angle should be calculated.
	 *  \param potentialFlags These indicate for which source indices the lens potential at
	 *                        the image points should be calculated.
	 *  */
	virtual void getTotalCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, std::vector<bool> &potentialFlags) const = 0;

	/** These flags indicate what information from the input images data should be stored in 
	 *  the ProjectedImagesInterface instance (you can save some memory this way, but if this
	 *  isn't an issue just set everything to \c true, it doensn't have any other effect on the
	 *  performance) */
	virtual void getTotalStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const = 0;

	/** Same as the LensFitnessObject::getTotalCalcFlags, but this is about the \c massScaleImages
	 *  list from the LensFitnessObject::init function. */
	virtual void getShortCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, std::vector<bool> &potentialFlags) const = 0;

	/** Same as the LensFitnessObject::getTotalStoreFlags, but this is about the \c massScaleImages
	 *  list from the LensFitnessObject::init function. */
	virtual void getShortStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const = 0;
	///@}

	// These are used in the genome fitness calculation. They all depend on the derivatives of the
	// deflection angle, so they can only be used for sources that have the derivative calculation
	// enabled in the ...CalcFlags functions below.
	// These functions then use these derivatives to calculate other properties
	/** @name Derived calculation flags
	 *  
	 *  The following functions provide information to the BackProjectMatrixNew instances in the genetic
	 *  algorithm about what information should be calculated based on the deflection angle derivatives.
	 *  You'll need to enable calculation of the derivatives of the deflection angles for a source to be able
	 *  to calculate these properties.
	 */
	///@{
	/** Should only return false if none of the sources need the inverse magnifications calculated. */
	virtual bool totalNeedInverseMagnifications() const = 0;

	/** Should only return false if none of the sources need the shear components (gamma1 and gamma2) to be calculated. */
	virtual bool totalNeedShearComponents() const = 0;

	/** Should only return false if none of the sources need the convergence to be calculated. */
	virtual bool totalNeedConvergence() const = 0;

	/** If LensFitnessObject::totalNeedInverseMagnifications returns true, this function will be called
	 *  and should return a list (of equal length as \c images from the LensFitnessObject::init function)
	 *  indicating for which sources the inverse magnification calculation should be done. */
	virtual const std::vector<bool> *getTotalInverseMagnificationFlags() const = 0;

	/** If LensFitnessObject::totalNeedShearComponents returns true, this function will be called
	 *  and should return a list (of equal length as \c images from the LensFitnessObject::init function)
	 *  indicating for which sources the shear components calculation should be done. */
	virtual const std::vector<bool> *getTotalShearComponentFlags() const = 0;

	/** If LensFitnessObject::totalNeedConvergence returns true, this function will be called
	 *  and should return a list (of equal length as \c images from the LensFitnessObject::init function)
	 *  indicating for which sources the convergence calculation should be done. */
	virtual const std::vector<bool> *getTotalConvergenceFlags() const = 0;

	/** Same as LensFitnessObject::totalNeedInverseMagnifications, but for the images data used in the
	 *  mass scale fitness calculation (the \c massScaleImages from the LensFitnessObject::init function). */
	virtual bool shortNeedInverseMagnifications() const = 0;

	/** Same as LensFitnessObject::totalNeedShearComponents, but for the images data used in the
	 *  mass scale fitness calculation (the \c massScaleImages from the LensFitnessObject::init function). */
	virtual bool shortNeedShearComponents() const = 0;

	/** Same as LensFitnessObject::totalNeedConvergence, but for the images data used in the
	 *  mass scale fitness calculation (the \c massScaleImages from the LensFitnessObject::init function). */
	virtual bool shortNeedConvergence() const = 0;

	/** Same as LensFitnessObject::getTotalInverseMagnificationFlags, but for the images data used in the
	 *  mass scale fitness calculation (the \c massScaleImages from the LensFitnessObject::init function). */
	virtual const std::vector<bool> *getShortInverseMagnificationFlags() const = 0;

	/** Same as LensFitnessObject::getTotalShearComponentFlags, but for the images data used in the
	 *  mass scale fitness calculation (the \c massScaleImages from the LensFitnessObject::init function). */
	virtual const std::vector<bool> *getShortShearComponentFlags() const = 0;

	/** Same as LensFitnessObject::getTotalConvergenceFlags, but for the images data used in the
	 *  mass scale fitness calculation (the \c massScaleImages from the LensFitnessObject::init function). */
	virtual const std::vector<bool> *getShortConvergenceFlags() const = 0;
	///@}
};

} // end namespace

#endif // GRALE_LENSFITNESSOBJECT_H

