#include <grale/gridlensinversiongafactorybase.h>
#include <grale/gridlensinversiongafactoryparams.h>
#include <grale/gridlensinversiongenomebase.h>
#include <grale/imagesdataextended.h>
#include <grale/constants.h>
#include <grale/grid.h>
#include "lensfitnesssimplerectangles.h"
#include <mogal/gafactorymultiobjective.h>
#include <serut/vectorserializer.h>
#include <iostream>
#include <memory>
#include <stdlib.h>

using namespace grale;
using namespace std;
using namespace mogal;

class GridLensInversionGAFactory_Test : public GridLensInversionGAFactoryBase, public GAFactorySingleObjective
{
public:
	GridLensInversionGAFactory_Test() { }
	~GridLensInversionGAFactory_Test() { }
	LensFitnessObject *createFitnessObject() { return new LensFitnessSimpleRectangles(); }
	bool subInit(LensFitnessObject *pFitnessObject) { return true; }
	
	float getChanceMultiplier() { return 1.0f; }
	bool useAbsoluteMutation() { return true; }
	float getMutationAmplitude() { return 1.0; }
	
	void onGeneticAlgorithmStep(int generation, bool *pGenerationInfoChanged, bool *pStopAlgorithm)
	{
		*pStopAlgorithm = true;
	}

	Genome *selectPreferredGenome(const list<Genome *> &nonDominatedSet) const
	{
		return *(nonDominatedSet.begin());
	}
};

int main(int argc, char *argv[])
{
	setenv("GRALE_DEBUG_SEED", "12345", 1);

	double D_d = 1000*DIST_MPC;
	double D_s = 1600*DIST_MPC;
	double D_ds = 1200*DIST_MPC;
	shared_ptr<ImagesDataExtended> imgExt(new ImagesDataExtended(D_s, D_ds));
	
	if (!imgExt->create(3, false, false))
	{
		cerr << "Couldn't create images data: " << imgExt->getErrorString() << endl;
		return -1;
	}

	vector<vector<Vector2Dd>> points { { {-2.5,0}, {-1.5, 0}, {-2, 0.5}, {-2, -0.5} },
		{ {-20, 0}, {-22, 0}, {-20.5, 5}, {-20.5,-5} },
		{ {28, 0}, {27, 0}, {27, 7}, {27, -7} } };
	
	int idx = 0;
	for (const auto &img : points)
	{
		for (const auto &pt : img)
		{
			Vector2Dd realPt = pt;
			realPt *= ANGLE_ARCSEC;
			imgExt->addPoint(idx, realPt);
		}
		idx++;
	}

	vector<GridSquare> gridSquares; // TODO
	for (const auto &sq : vector<GridSquare> { { { -10.01, -10 }, 20 },
			                                   { { 10.01, -10 }, 20 },
											   { { 10.01, 10 }, 20 },
											   { { -10.01, 10 }, 20 },
											   { { 5.01, 0 }, 10 },
											   { { -5.01, 0 }, 10 },
											   })
	{
		double size = sq.getSize() * ANGLE_ARCSEC;
		auto ctr = sq.getCenter();
		ctr *= ANGLE_ARCSEC;
		gridSquares.push_back( { ctr, size } );
	}

	double massScale = 1e14*MASS_SOLAR;
	ConfigurationParameters *pFitnessObjectParams = nullptr; 

	for (auto useWeights : vector<bool> { false, true })
	{
		for (auto allowNeg : vector<bool> { false, true })
		{
			for (auto basisFunction : vector<GridLensInversionParameters::BasisFunctionType> { GridLensInversionParameters::PlummerBasis,
					GridLensInversionParameters::GaussBasis, GridLensInversionParameters::SquareBasis })
			{
				for (auto sheetSearchType : vector<GridLensInversionParameters::MassSheetSearchType> { GridLensInversionParameters::NoSheet,
						GridLensInversionParameters::Genome })
				{
					for (auto wideSearch : vector<bool> { false, true })
					{
						cout << endl;
						cout << "-------------------------------------------------------------------------------------------------" << endl;
						cout << "useWeights: " << useWeights << " allowNeg: " << allowNeg << " basisFunction: " << (int)basisFunction
							 << " sheetSearchType: " << (int)sheetSearchType << " wideSearch: " << wideSearch << endl;

						GridLensInversionGAFactoryParams origParams(1, // maxgenerations
							{ imgExt }, gridSquares, D_d, 0, // z_d
							massScale,
							useWeights, basisFunction, allowNeg, nullptr, // pBaseLens
							sheetSearchType, pFitnessObjectParams,
							ScaleSearchParameters(wideSearch));

						unique_ptr<GridLensInversionGAFactoryParams> pParamsCopy(origParams.createCopy());
						serut::VectorSerializer ser;

						GridLensInversionGAFactoryParams params;
						if (!pParamsCopy->write(ser) || !params.read(ser))
						{
							cerr << "Unable to serialize/deserialize inversion parameters" << endl;
							return -1;
						}

						GridLensInversionGAFactory_Test factory;
						if (!factory.init(&params))
						{
							cerr << "Couldn't init factory: " << factory.getErrorString() << endl;
							return -1;
						}

						auto pGenome = unique_ptr<GridLensInversionGenomeBase>(static_cast<GridLensInversionGenomeBase*>(factory.createNewGenome()));
						pGenome->calculateFitness();

						cout << "masses:";
						for (auto &m : pGenome->getBasisFunctionWeights())
							cout << " " << m;
						cout << endl;
						cout << "sheet:";
						for (auto &s : pGenome->getSheetValues())
							cout << " " << s;
						cout << endl;
						cout << "scale: " << pGenome->getScaleFactor() << endl;
					}
				}
			}
		}
	}
	return 0;
}
