#include "potentialgridlens.h"
#include "constants.h"
#include <eatk/vectorgenomefitness.h>
#include <eatk/jadeevolver.h>
#include <eatk/vectordifferentialevolution.h>
#include <eatk/mersennerandomnumbergenerator.h>
#include <eatk/evolutionaryalgorithm.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <eatk/multithreadedpopulationfitnesscalculation.h>
#include <eatk/stopcriterion.h>
#include <random>
#include <fstream>

using namespace grale;
using namespace eatk;
using namespace std;
using namespace errut;

class PotentialGridLensBase
{
public:
	PotentialGridLensBase(double Dd, Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY);
	~PotentialGridLensBase() { cleanup(); }

	const vector<double> &values() const { return m_values; }
	vector<double> &values() { return m_values; }

	// Must be called when values changed!
	errut::bool_t init();

	errut::bool_t getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	errut::bool_t getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	errut::bool_t getSurfaceMassDensity(Vector2D<double> theta, double &dens) const;
	errut::bool_t getProjectedPotential(Vector2D<double> theta, double *pPotentialValue) const;
private:
	void cleanup();
	
	const double m_Dd;
	const Vector2Dd m_bottomLeft, m_topRight;
	const int m_numX, m_numY;

	gsl_interp2d *m_pInterp;
	gsl_interp_accel *m_pXAccel, *m_pYAccel;
	std::vector<double> m_x, m_y;
	std::vector<double> m_values;
};

PotentialGridLensBase::PotentialGridLensBase(double Dd, Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY)
	: m_Dd(Dd), m_bottomLeft(bottomLeft), m_topRight(topRight), m_numX(numX), m_numY(numY)
{
	m_values.resize(numX*numY);
	m_pInterp = nullptr;
	m_pXAccel = nullptr;
	m_pYAccel = nullptr;

	for (int x = 0 ; x < numX ; x++)
		m_x.push_back((topRight.getX()-bottomLeft.getX())/(double)(numX-1) * x + bottomLeft.getX());

	for (int y = 0 ; y < numY ; y++)
		m_y.push_back((topRight.getY()-bottomLeft.getY())/(double)(numY-1) * y + bottomLeft.getY());
}

errut::bool_t PotentialGridLensBase::init()
{
	cleanup();

	if (m_bottomLeft.getX() >= m_topRight.getX() || m_bottomLeft.getY() >= m_topRight.getY())
		return "Corners need to be ordered correctly";

	if (m_numX*m_numY != m_values.size())
		return "Number of potential values doesn't match grid dimensions";

	if (m_numX < gsl_interp2d_type_min_size(gsl_interp2d_bicubic) ||
		m_numY < gsl_interp2d_type_min_size(gsl_interp2d_bicubic))
		return "Not enough points in grid in at least one dimension";

	m_pInterp = gsl_interp2d_alloc(gsl_interp2d_bicubic, m_numX, m_numY);
	if (!m_pInterp)
		return "Unable to allocate GSL interpolation workspace";

	m_pXAccel = gsl_interp_accel_alloc();
	m_pYAccel = gsl_interp_accel_alloc();

	int err = gsl_interp2d_init(m_pInterp, m_x.data(), m_y.data(), m_values.data(), m_numX, m_numY);
	if (err < 0)
		return "Unable to initialize GSL interpolation workspace: error " + to_string(err);

	return true;
}

void PotentialGridLensBase::cleanup()
{
	if (m_pXAccel)
		gsl_interp_accel_free(m_pXAccel);
	if (m_pYAccel)
		gsl_interp_accel_free(m_pYAccel);
	if (m_pInterp)
		gsl_interp2d_free(m_pInterp);

	m_pInterp = nullptr;
	m_pXAccel = nullptr;
	m_pYAccel = nullptr;	
}

bool_t PotentialGridLensBase::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	double ax = 0, ay = 0;
	auto getDeriv = [this, theta](
			int (*functionName)(const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d),
			double &result) -> bool_t
	{
		int err = functionName(m_pInterp, m_x.data(), m_y.data(), m_values.data(),
				            theta.getX(), theta.getY(), m_pXAccel, m_pYAccel, &result);
		if (err < 0)
			return "Error getting derivative of potential: GSL error code " + to_string(err);
		return true;
	};

	bool_t r = getDeriv(gsl_interp2d_eval_deriv_x_e, ax);
	if (!r)
		return r;
	r = getDeriv(gsl_interp2d_eval_deriv_y_e, ay);
	if (!r)
		return r;

	*pAlpha = Vector2Dd(ax, ay);
	return true;
}

bool_t PotentialGridLensBase::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	auto getDerivDeriv = [this, theta](
			int (*functionName)(const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d),
			double &result) -> bool_t
	{
		int err = functionName(m_pInterp, m_x.data(), m_y.data(), m_values.data(), 
				               theta.getX(), theta.getY(), m_pXAccel, m_pYAccel, &result);
		if (err < 0)
			return "Error getting second derivative of potential: GSL error code " + to_string(err);
		return true;
	};

	bool_t r = getDerivDeriv(gsl_interp2d_eval_deriv_xx_e, axx);
	if (!r)
		return r;

	r = getDerivDeriv(gsl_interp2d_eval_deriv_yy_e, ayy);
	if (!r)
		return r;

	r = getDerivDeriv(gsl_interp2d_eval_deriv_xy_e, axy);
	if (!r)
		return r;

	return true;
}

bool_t PotentialGridLensBase::getSurfaceMassDensity(Vector2D<double> theta, double &dens) const
{
	double axx, ayy, axy;

	bool_t r = getAlphaVectorDerivatives(theta, axx, ayy, axy);
	if (!r)
		return r;

	double kappa = 0.5*(axx + ayy);
	double sigmaCrit = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*m_Dd);
	dens = kappa*sigmaCrit;
	return true;
}

bool_t PotentialGridLensBase::getProjectedPotential(Vector2D<double> theta, double *pPotentialValue) const
{
	double result;
	int err = gsl_interp2d_eval_e(m_pInterp, m_x.data(), m_y.data(), m_values.data(),
	                              theta.getX(), theta.getY(), m_pXAccel, m_pYAccel, &result);
	if (err < 0)
		return "Unable to get projected potential at requested point: GSL error code " + to_string(err);

	*pPotentialValue = result;
	return true;
}

class PotentialGrid
{
public:
	PotentialGrid(const PotentialGridLens &l);

	void checkSameSettings(const PotentialGrid &g);
	void checkMaskSize(const vector<uint8_t> &mask);

	double getWidth() const { return m_topRight.getX() - m_bottomLeft.getX(); }
	double getHeight() const { return m_topRight.getY() - m_bottomLeft.getY(); }

	Vector2Dd topRight() const { return m_topRight; }
	Vector2Dd bottomLeft() const { return m_bottomLeft; }

	size_t getNumX() const { return m_numX; }
	size_t getNumY() const { return m_numY; }

	const vector<double> &values() const { return m_values; }
	double getMax() const { return m_maxValue; }

	double getDeltaPosSize() const { return m_deltaPosSize; }
private:
	vector<double> m_values;
	size_t m_numX, m_numY;
	Vector2Dd m_bottomLeft, m_topRight;
	double m_maxValue;

	double m_deltaPosSize;
};

PotentialGrid::PotentialGrid(const PotentialGridLens &l)
{
	auto pParams = l.getLensParameters();
	if (!pParams)
		throw runtime_error("Unexpected: lens parameters are null");

	auto &params = dynamic_cast<const PotentialGridLensParams &>(*pParams);

	m_values = params.getValues();
	m_numX = params.getNumX();
	m_numY = params.getNumY();
	m_topRight = params.getTopRight();
	m_bottomLeft = params.getBottomLeft();

	// Make the minimum value 0;
	double minValue = *min_element(m_values.begin(), m_values.end());
	for (double &x : m_values)
		x -= minValue;

	// get the max value, for a scale
	m_maxValue = *max_element(m_values.begin(), m_values.end());

	double pixelW = getWidth()/(getNumX()-1);
	double pixelH = getHeight()/(getNumY()-1);
	double pixelSize = sqrt(pixelW*pixelW+pixelH*pixelH);
	m_deltaPosSize = pixelSize/100.0;
}

void PotentialGrid::checkSameSettings(const PotentialGrid &g)
{
	if (!(m_bottomLeft == g.m_bottomLeft && m_topRight == g.m_topRight))
		throw "Grid corners differ";
	if (m_numX != g.m_numX || m_numY != g.m_numY)
		throw "Grid dimensions differ";
}

void PotentialGrid::checkMaskSize(const vector<uint8_t> &mask)
{
	if (mask.size() != m_numX*m_numY)
		throw runtime_error("Mask size is not compatible with potential grid size");
}

vector<uint8_t> createMask(Vector2Dd bottomLeft, Vector2Dd topRight, size_t numX, size_t numY, double r1, double r2)
{
	vector<uint8_t> data;
	
	double r1Sq = r1*r1;
	double r2Sq = r2*r2;
	for (size_t Y = 0 ; Y < numY ; Y++)
	{
		double y = ((double)Y/((double)numY-1.0)) * (topRight.getY() - bottomLeft.getY()) + bottomLeft.getY();
		for (size_t X = 0 ; X < numX ; X++)
		{
			double x = ((double)X/((double)numX-1.0)) * (topRight.getX() - bottomLeft.getX()) + bottomLeft.getX();
			double rSq = x*x + y*y;
			uint8_t m = 0;
			
			if (rSq < r1Sq)
				m = 0; // first lens file, inner region
			else if (rSq > r2Sq)
				m = 1; // second lens file, outer region
			else
				m = 2; // interpolation region			

			data.push_back(m);
		}
	}

	return data;
}

class PotentialBuilderShared
{
public:
	PotentialBuilderShared(const string &fileName1, const string &fileName2, double r1, double r2);

	size_t getGenomeSize() const { return m_genomeSize; }
	
	const PotentialGrid &getGrid(size_t i) { return *m_grids[i%2]; }
	const GravitationalLens &getLens(size_t i) { return *m_lenses[i%2]; }
	const vector<Vector2Dd> &getInterpolationPositions() const { return m_interpPos; }
	const vector<size_t> &getInterpolationIndices() const { return m_interpIndices; }
	const vector<uint8_t> &getMask() const { return m_mask; }

	double getLensDistance() const { return m_Dd; }
	double getR1() const { return m_r1; }
	double getR2() const { return m_r2; }
private:
	unique_ptr<GravitationalLens> m_lenses[2];
	unique_ptr<PotentialGrid> m_grids[2];
	vector<uint8_t> m_mask;
	size_t m_genomeSize;
	double m_Dd;
	double m_r1, m_r2;

	vector<Vector2Dd> m_interpPos;
	vector<size_t> m_interpIndices;
};

class PotentialBuilder
{
public:
	PotentialBuilder(const shared_ptr<PotentialBuilderShared> &sharedBuilderInfo);
	
	const GravitationalLens &getLens(size_t i) { return m_builder->getLens(i); }
	const vector<Vector2Dd> &getInterpolationPositions() const { return m_builder->getInterpolationPositions(); }

	const PotentialGridLensBase &buildLens(const vector<double> &genome);
	unique_ptr<GravitationalLens> buildLensFull(const vector<double> &genome);
private:
	shared_ptr<PotentialBuilderShared> m_builder;

	unique_ptr<PotentialGridLensBase> m_lensInner, m_lensOuter, m_lensInterp;
};

PotentialBuilder::PotentialBuilder(const shared_ptr<PotentialBuilderShared> &sharedBuilderInfo)
	: m_builder(sharedBuilderInfo)
{
	auto &grid = m_builder->getGrid(0);

	auto createLensBase = [this,&grid]()
	{
		return make_unique<PotentialGridLensBase>(m_builder->getLensDistance(), grid.bottomLeft(), grid.topRight(), grid.getNumX(), grid.getNumY());
	};

	m_lensInner = createLensBase();
	m_lensOuter = createLensBase();
	m_lensInterp = createLensBase();
}


PotentialBuilderShared::PotentialBuilderShared(const string &fileName1, const string &fileName2, double r1, double r2)
{
	m_r1 = r1;
	m_r2 = r2;

	string errStr;

	if (!GravitationalLens::load(fileName1, m_lenses[0], errStr) ||
		!GravitationalLens::load(fileName2, m_lenses[1], errStr))
		throw runtime_error("Can't load on of the lenses: " + errStr);

	if (m_lenses[0]->getLensDistance() != m_lenses[1]->getLensDistance())
		throw runtime_error("Lenses have different distances");

	m_Dd = m_lenses[0]->getLensDistance();

	auto &potLens1 = dynamic_cast<const PotentialGridLens &>(*m_lenses[0]);
	auto &potLens2 = dynamic_cast<const PotentialGridLens &>(*m_lenses[1]);

	m_grids[0] = make_unique<PotentialGrid>(potLens1);
	m_grids[1] = make_unique<PotentialGrid>(potLens2);

	Vector2Dd topRight = m_grids[0]->topRight();
	Vector2Dd bottomLeft = m_grids[0]->bottomLeft();
	size_t numX = m_grids[0]->getNumX();
	size_t numY = m_grids[0]->getNumY();

	m_mask = createMask(bottomLeft, topRight, numX, numY, r1, r2);

	m_grids[1]->checkSameSettings(*m_grids[0]);
	m_grids[1]->checkMaskSize(m_mask);

	// Get which non-zero mask entry maps to which actual position
	vector<size_t> maskPositions;
	for (size_t i = 0 ; i < m_mask.size() ; i++)
		if (m_mask[i] == 2)
			maskPositions.push_back(i);
	
	cout << "Found " << maskPositions.size() << " interpolated mask positions" << endl;
	if (!maskPositions.size())
		throw runtime_error("No positions found that should be interpolated");

	m_genomeSize = 5; // 1 for potential offset, 2 gradients for first map, 2 gradients for second and the interpolated positions

	// Calculate the positions for which the kappa gradients should be calculated

	// TODO: use more positions?
	for (size_t i = 0 ; i < m_mask.size() ; i++)
	{
		if (m_mask[i] == 2) // marked as interpolation
		{
			size_t X = i%m_grids[0]->getNumX();
			size_t Y = i/m_grids[0]->getNumY();
			assert(Y < numY);

			double Xfrac = (double)X/((double)numX-1);
			double Yfrac = (double)Y/((double)numY-1);

			double x = bottomLeft.getX() + (topRight.getX() - bottomLeft.getX())*Xfrac;
			double y = bottomLeft.getY() + (topRight.getY() - bottomLeft.getY())*Yfrac;

			m_interpPos.push_back({x,y});
			m_interpIndices.push_back(i);
		}
	}
}

double getSeriesValue(const vector<double> &coeffs, double x)
{
	double xPow = 1.0;
	double sum = 0;
	for (auto a : coeffs)
	{
		sum += xPow * a;
		xPow *= x;
	}
	return sum;
}

double getInterpolatedValue(double v1, double v2, double x)
{
	//static const vector<double> coeffs = { 1, 0, 0, 0, -35, 84, -70, 20 }; even smoother step
	static const vector<double> coeffs = { 1, 0, 0, -10, 15, -6 }; // smootherstep
	//static const vector<double> defaultCoeffs = { 1, 0, -3, 2 }; // smoothstep
	//static const vector<double> oeffs = { 1, -1 }; // linear
	double y = getSeriesValue(coeffs, x);
	return v1*y + (1.0-y)*v2;
}

const PotentialGridLensBase &PotentialBuilder::buildLens(const vector<double> &genome)
{
	assert(genome.size() == 5);
	double diffx1 = genome[0];
	double diffy1 = genome[1];
	double diffx2 = genome[2];
	double diffy2 = genome[3];
	double diffpot = genome[4];

	vector<double> *values[2] { &m_lensInner->values(), &m_lensOuter->values() };
	vector<double> &valuesNew = m_lensInterp->values();
	
	values[0]->clear();
	values[1]->clear();
	valuesNew.clear();

	const PotentialGrid &grid = m_builder->getGrid(0);
	const PotentialGrid &grid2 = m_builder->getGrid(1);

	Vector2Dd topRight = grid.topRight();
	Vector2Dd bottomLeft = grid.bottomLeft();
	size_t numX = grid.getNumX();
	size_t numY = grid.getNumY();

	auto fromGradient = [numX, numY, bottomLeft, topRight](size_t X, size_t Y, double derivx, double derivy)
	{
		double Xfrac = (double)X/((double)numX-1);
		double Yfrac = (double)Y/((double)numY-1);

		double dx = bottomLeft.getX() + (topRight.getX() - bottomLeft.getX())*Xfrac;
		double dy = bottomLeft.getY() + (topRight.getY() - bottomLeft.getY())*Yfrac;

		return dx*derivx + dy*derivy;
	};

	for (size_t Y = 0, i = 0 ; Y < numY ; Y++)
	{
		for (size_t X = 0 ; X < numX ; X++, i++)
		{
			values[0]->push_back(grid.values()[i] + fromGradient(X, Y, diffx1, diffy1));
			values[1]->push_back(grid2.values()[i] + diffpot + fromGradient(X, Y, diffx2, diffy2));
			auto maskValue = m_builder->getMask()[i];
			if (maskValue != 2) // not interpolated, store one of the values in the final grid
				valuesNew.push_back(values[maskValue]->back());
			else
				valuesNew.push_back(-12345); // Will be filled in later			
		}
	}

	// Process the new values for these lenses
	m_lensInner->init();
	m_lensOuter->init();

	const vector<Vector2Dd> &interpPos = m_builder->getInterpolationPositions();
	const vector<size_t> &interpIdx = m_builder->getInterpolationIndices();

	assert(interpPos.size() == interpIdx.size());

	auto getPotential = [](const PotentialGridLensBase &lens, Vector2Dd pos) -> double
	{
		double phi;
		bool_t r = lens.getProjectedPotential(pos,&phi);
		if (!r)
			throw runtime_error("Can't get potential: " + r.getErrorString());
		return phi;
	};

	double r1 = m_builder->getR1();
	double r2 = m_builder->getR2();

	for (size_t j = 0 ; j < interpIdx.size() ; j++)
	{
		auto pos = interpPos[j];
		auto idx = interpIdx[j];

		double r = pos.getLength();
		double x = (r - r1)/(r2 - r1);

		if (x < 0 || x > 1)
			throw runtime_error("Internal error: x = " + to_string(x) + " x - 1 = " + to_string(x-1));

		assert(idx >= 0 && idx < valuesNew.size());
		assert(valuesNew[idx] == -12345);
		
		valuesNew[idx] = getInterpolatedValue(getPotential(*m_lensInner, pos), getPotential(*m_lensOuter, pos), x);
	}

	m_lensInterp->init();
	return *m_lensInterp;
}

unique_ptr<GravitationalLens> PotentialBuilder::buildLensFull(const vector<double> &genome)
{
	buildLens(genome); // Fills in m_lensInterp

	auto lens = make_unique<PotentialGridLens>();
	auto &grid = m_builder->getGrid(0);

	PotentialGridLensParams params(grid.bottomLeft(), grid.topRight(), m_lensInterp->values(), grid.getNumX(), grid.getNumY());
	if (!lens->init(m_builder->getLensDistance(), &params))
		throw runtime_error("Can't init lens: " + lens->getErrorString());
	return lens;
}


class Creator : public IndividualCreation
{
public:
	Creator(const shared_ptr<RandomNumberGenerator> &rng, const shared_ptr<PotentialBuilderShared> &builder)
		: m_rng(rng), m_builder(builder) { }

	shared_ptr<Genome> createInitializedGenome() override;
	shared_ptr<Fitness> createEmptyFitness() override { return make_shared<ValueFitness<double>>(); }
private:
	shared_ptr<RandomNumberGenerator> m_rng;
	shared_ptr<PotentialBuilderShared> m_builder;
};

shared_ptr<Genome> Creator::createInitializedGenome()
{
	double potScale = std::max(m_builder->getGrid(0).getMax(),m_builder->getGrid(1).getMax());
	double sizeScale = std::max(std::max(m_builder->getGrid(0).getHeight(),m_builder->getGrid(1).getHeight()),
								std::max(m_builder->getGrid(0).getWidth(),m_builder->getGrid(1).getWidth()));

	auto genome = make_shared<VectorGenome<double>>(m_builder->getGenomeSize());
	vector<double> &values = genome->getValues();
	
	assert(values.size() == 5);

	double gradSize = (potScale/sizeScale)/10.0;
	double potSize = potScale/10.0;
	
	values[0] = m_rng->getRandomDouble(-gradSize, gradSize);
	values[1] = m_rng->getRandomDouble(-gradSize, gradSize);
	values[2] = m_rng->getRandomDouble(-gradSize, gradSize);
	values[3] = m_rng->getRandomDouble(-gradSize, gradSize);
	values[4] = m_rng->getRandomDouble(-potSize, potSize);
	
	return genome;
}

class FitnessCalc : public GenomeFitnessCalculation
{
public:
	FitnessCalc(const shared_ptr<PotentialBuilder> &builder) : m_builder(builder) { }

	double calculateFitness(const PotentialGridLensBase &lens)
	{
		bool_t r;
		double sum = 0;
		auto &lens1 = m_builder->getLens(0);
		auto &lens2 = m_builder->getLens(1);
		for (auto pos : m_builder->getInterpolationPositions())
		{
			double kappa;

			r = lens.getSurfaceMassDensity(pos, kappa);
			if (!r)
				throw runtime_error("Can't get surface mass density: " + r.getErrorString());
			
			double kappa1 = lens1.getSurfaceMassDensity(pos);
			double kappa2 = lens2.getSurfaceMassDensity(pos);
			double dk1 = kappa-kappa1;
			double dk2 = kappa-kappa2;

			sum += dk1*dk1 + dk2*dk2;
		}

		sum /= m_builder->getInterpolationPositions().size();

		return std::sqrt(sum);
	}

	bool_t calculate(const Genome &g, Fitness &f)
	{
		const VectorGenome<double> &genome = static_cast<const VectorGenome<double> &>(g);
		ValueFitness<double> &fitness = static_cast<ValueFitness<double>&>(f);

		auto &lens = m_builder->buildLens(genome.getValues());
		double value = calculateFitness(lens);

		fitness.setValue(value);
		fitness.setCalculated();
		return true;
	}
private:
	shared_ptr<PotentialBuilder> m_builder;
};

void save(const shared_ptr<Individual> &best, PotentialBuilder &builder, const string &fileName)
{
	const VectorGenome<double> &genome = static_cast<VectorGenome<double>&>(best->genomeRef());
	auto &values = genome.getValues();
	auto lens = builder.buildLensFull(values);

	if (!lens->save(fileName))
		throw runtime_error("Couldn't save solution: " + lens->getErrorString());

	cout << "Saved a solution to " << fileName << endl;
	cout << "Genome is: " << endl;
	cout << "  gx1: " << values[0]/ANGLE_ARCSEC << endl;
	cout << "  gy1: " << values[1]/ANGLE_ARCSEC << endl;
	cout << "  gx2: " << values[2]/ANGLE_ARCSEC << endl;
	cout << "  gy2: " << values[3]/ANGLE_ARCSEC << endl;
	cout << "  offset: " << values[4]/(ANGLE_ARCSEC*ANGLE_ARCSEC) << endl;
}

class Stop : public FixedGenerationsStopCriterion
{
public:
	Stop(size_t maxGen, const shared_ptr<PotentialBuilder> &builder, const string &tmpFilePrefix, size_t tmpFileInterval) 
		: FixedGenerationsStopCriterion(maxGen), m_builder(builder),
		  m_tmpFilePrefix(tmpFilePrefix), m_tmpFileInterval(tmpFileInterval) { }
	bool_t analyze(const std::vector<std::shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop) override;
private:
	shared_ptr<PotentialBuilder> m_builder;
	const string m_tmpFilePrefix;
	const size_t m_tmpFileInterval;
};

bool_t Stop::analyze(const std::vector<std::shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop)
{
	cout << "Best on generation " << generationNumber << ": " << currentBest[0]->fitnessRef().toString() << endl;
	if (generationNumber%m_tmpFileInterval == 1) // save intermediate
		save(currentBest[0], *m_builder, m_tmpFilePrefix + to_string(generationNumber) + ".lensdata");
	
	return FixedGenerationsStopCriterion::analyze(currentBest, generationNumber, shouldStop);
}

int mainCxx(const vector<string> &args)
{
	if (args.size() > 11)
		throw runtime_error("Too many arguments");

	string fileName1 = args.at(0);
	string fileName2 = args.at(1);
	double r1 = stod(args.at(2)) * ANGLE_ARCSEC;
	double r2 = stod(args.at(3)) * ANGLE_ARCSEC;

	if (r1 < 0 || r2 <= r1)
		throw runtime_error("Invalid r1 or r2");

	string seedStr = args.at(4);
	size_t numThreads = stoul(args.at(5));
	size_t popSize = stoul(args.at(6));
	size_t maxGen = stoul(args.at(7));
	string targetFileName = args.at(8);
	string tmpFileNamePrefix = args.at(9);
	size_t tmpFileInterval = stoul(args.at(10));
	
	if (numThreads > 128)
		throw runtime_error("Too many threads");
	if (numThreads == 0)
		throw runtime_error("No threads");

	shared_ptr<PotentialBuilderShared> builderShared = make_shared<PotentialBuilderShared>(fileName1, fileName2, r1, r2);

	vector<shared_ptr<PotentialBuilder>> builders;
	vector<shared_ptr<GenomeFitnessCalculation>> genomeCalculators;
	for (size_t i = 0 ; i < numThreads ; i++)
	{
		builders.push_back(make_shared<PotentialBuilder>(builderShared));
		genomeCalculators.push_back(make_shared<FitnessCalc>(builders.back()));
	}

	random_device rd;
	auto seed = rd();
	if (seedStr != "random")
		seed = stoi(seedStr);

	cout << "Using seed " << seed << endl;

	auto rng = make_shared<MersenneRandomNumberGenerator>(seed);
	auto mut = make_shared<VectorDifferentialEvolutionMutation<double>>();
	auto cross = make_shared<VectorDifferentialEvolutionCrossover<double>>(rng);
	auto comp = make_shared<ValueFitnessComparison<double>>();

	JADEEvolver evolver(rng, mut, cross, comp);
	EvolutionaryAlgorithm ea;

	shared_ptr<PopulationFitnessCalculation> popFitnessCalc;
	if (numThreads == 1)
		popFitnessCalc = make_shared<SingleThreadedPopulationFitnessCalculation>(genomeCalculators[0]);
	else
	{
		auto fc = make_shared<MultiThreadedPopulationFitnessCalculation>();
		popFitnessCalc = fc;
		bool_t r = fc->initThreadPool(genomeCalculators);
		if (!r)
			throw runtime_error("Can't init thread pool: " + r.getErrorString());
	}

	Creator creator(rng, builderShared);
	Stop stop(maxGen, builders[0], tmpFileNamePrefix, tmpFileInterval);
	
	bool_t r = ea.run(creator, evolver, *popFitnessCalc, stop, popSize, popSize, popSize*2);
	if (!r)
		throw runtime_error("Error running EA: " + r.getErrorString());

	cout << "EA Finished" << endl;
	save(evolver.getBestIndividuals()[0], *builders[0], targetFileName);

	cout << "Solution saved to " << targetFileName << endl;
	return 0;
}

int main(int argc, char *argv[])
{
	try
	{
		std::vector<std::string> args;
		for (int i = 1 ; i < argc ; i++)
			args.push_back(argv[i]);

		return mainCxx(args);
	}
	catch(const exception &e)
	{
		cerr << "Error: " << e.what() << endl;
		cerr << "Usage:\n";
		cerr << "  mergepotentials lens1inner lens2outer r1arcsec r2arcsec seed|random numthreads popsize maxgen savefile tmpprefix tmpinterval\n\n";
	}
	catch(...)
	{
		cerr << "Caught some unknown exception" << endl;
	}
	return -1;
}

