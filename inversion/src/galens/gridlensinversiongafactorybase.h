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

#ifndef GRALE_GRIDLENSINVERSIONGAFACTORYBASE_H

#define GRALE_GRIDLENSINVERSIONGAFACTORYBASE_H

#include "graleconfig.h"
#include "randomnumbergenerator.h"
#include "backprojectmatrixnew.h"
#include "vector2d.h"
#include <mogal/gafactorydefaults.h>
#include <vector>

namespace grale
{

class GridLensInversionGAFactoryParams;
class GridLensInversionGenomeBase;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

// NOTE: the virtual inheritance is again very important!
class GRALE_IMPORTEXPORT GridLensInversionGAFactoryBase : public virtual mogal::GAFactory
{
public:
	GridLensInversionGAFactoryBase();
	~GridLensInversionGAFactoryBase();

	mogal::GAFactoryParams *createParamsInstance() const;
	const mogal::GAFactoryParams *getCurrentParameters() const;

	bool init(const mogal::GAFactoryParams *p);

	mogal::Genome *createNewGenome() const;

	size_t getMaximalFitnessSize() const							{ return sizeof(float)*(2+getNumberOfFitnessComponents()); } // one for scale factor, one for sheet scale factor some for fitness
	size_t getMaximalGenomeSize() const							{ return (m_numMasses+1)*sizeof(float); } // one extra for the mass sheet value

	bool writeGenome(serut::SerializationInterface &si, const mogal::Genome *g) const;
	bool readGenome(serut::SerializationInterface &si, mogal::Genome **g) const;
	bool writeGenomeFitness(serut::SerializationInterface &si, const mogal::Genome *g) const;
	bool readGenomeFitness(serut::SerializationInterface &si, mogal::Genome *g) const;
	bool writeCommonGenerationInfo(serut::SerializationInterface &si) const;
	bool readCommonGenerationInfo(serut::SerializationInterface &si);

	bool hasFloatingPointFitnessValues() const 						{ return true; }

	const mogal::RandomNumberGenerator *getRandomNumberGenerator() const			{ return &m_rndGen; }
	const std::vector<Vector2D<float> > &getPlummerPositions() const			{ return m_plummerPositions; }
	const std::vector<float> &getSquareSizes() const					{ return m_squareSizes; }
	float getGridSize() const								{ return m_gridSize; }
	Vector2D<float> getGridCenter() const							{ return m_gridCenter; }
	bool allowNegativeValues() const							{ return m_allowNegativeValues; }
	const float *getMassWeights() const							{ return &(m_massWeights[0]); }
	float getSheetScale() const								{ return m_sheetScale; }
	bool useLoopSheet() const								{ return m_useLoopSheet; }
	
	void getGenomeCalculationParameters(float &startfactor, float &stopfactor, int &numiterationsteps, int &numiterations, int &numiterationsteps2) const;
	void getGenomeSheetCalculationParameters(float &startfactor, float &stopfactor) const;
	bool useLogarithmicScaleSearch() const { return true; }

	virtual float getChanceMultiplier() = 0;
	virtual bool useAbsoluteMutation() = 0;
	virtual float getMutationAmplitude() = 0;

	GravitationalLens *createLens(const std::vector<float> &masses, float sheetValue, float scaleFactor, double *pTotalMass, std::string &errStr) const;

	DeflectionMatrix *getDeflectionMatrix()							{ return m_pDeflectionMatrix; }
	BackProjectMatrixNew *getShortBackProjectMatrix()					{ return m_pShortBPMatrix; }
	BackProjectMatrixNew *getTotalBackProjectMatrix()					{ return m_pTotalBPMatrix; }
	LensFitnessObject *getFitnessObject()							{ return m_pFitnessObject; }
protected:
	int getMaximumNumberOfGenerations() const						{ return m_maxGenerations; }

	virtual LensFitnessObject *createFitnessObject() = 0;

	// This should at least set the number of fitness components
	virtual bool subInit(LensFitnessObject *pFitnessObject) = 0;
#ifdef SHOWEVOLUTION
	void onSortedPopulation(const std::vector<mogal::GenomeWrapper> &population);
#endif // SHOWEVOLUTION

	void onCurrentBest(const std::list<mogal::Genome *> &bestGenomes) override;
private:
	bool localSubInit(double z_d, const std::vector<ImagesDataExtended *> &images, 
	                  const std::vector<std::pair<GravitationalLens *, Vector2D<double> > > &basisLenses,
					  const GravitationalLens *pBaseLens, bool useSheet, 
					  const ConfigurationParameters *pFitnessObjectParams);
	double getAngularScale() const								{ return m_pShortBPMatrix->getAngularScale(); }

	RandomNumberGenerator m_rndGen;
	GridLensInversionGAFactoryParams *m_pCurrentParams;
	int m_numMasses, m_maxGenerations;
	bool m_allowNegativeValues;
	bool m_useGenomeSheet;
	bool m_useLoopSheet;
	float m_sheetScale;
	bool m_wideSearch;

	std::vector<Vector2D<float> > m_plummerPositions;
	std::vector<float> m_squareSizes;
	std::vector<float> m_massWeights;
	float m_gridSize;
	Vector2D<float> m_gridCenter;

	std::vector<std::pair<GravitationalLens *, Vector2D<double> > > m_basisLenses;

	LensFitnessObject *m_pFitnessObject;
	DeflectionMatrix *m_pDeflectionMatrix;
	BackProjectMatrixNew *m_pShortBPMatrix, *m_pTotalBPMatrix;
};

} // end namespace

#endif // GRALE_GRIDLENSINVERSIONGAFACTORYBASE_H

