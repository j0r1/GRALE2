#include "graleconfig.h"
#include "lensgacalculatorregistry.h"
#include "lensinversiongafactorysingleplanecpu.h"
#include "lensinversiongafactorymultiplanegpu.h"
#include "lensfitnessobject.h"
#include "lensfitnessgeneral.h"
#include "lensgacalculatorregistry.h"
#include <serut/vectorserializer.h>

using namespace std;
using namespace errut;

namespace grale
{

template<class T, class P1>
class GAFactoryWrapperLensGAGenomeCalculator;

template<class T, class P1>
class GAFactoryHelper : public T
{
public:
	GAFactoryHelper(unique_ptr<LensFitnessObject> fitObj)
		: m_fitObj(move(fitObj)), m_pParent(nullptr) { }
	GAFactoryHelper() { }

	void setParent(GAFactoryWrapperLensGAGenomeCalculator<T,P1> *pParent) { m_pParent = pParent; }

	// Fr now, we're not going to create new one, just use the previously set one
	LensFitnessObject *createFitnessObject() override
	{
		// cout << "GAFactoryHelper: createFitnessObject " << (void*)m_fitObj.get() << endl;
		if (!m_fitObj.get())
		{
			T::setErrorString("LensFitnessObject was already retrieved");
			return nullptr;
		}
		auto obj = m_fitObj.release(); // The caller will take care of this
		return obj;
	}

	void sendMessage(const string &s) const override
	{
		if (!m_pParent)
			cerr << s << endl;
		else
			m_pParent->log(s);
	}
private:
	unique_ptr<LensFitnessObject> m_fitObj;
	GAFactoryWrapperLensGAGenomeCalculator<T,P1> *m_pParent;
};

template<class T, class P1>
class GAFactoryWrapperLensGAGenomeCalculator: public LensGAGenomeCalculator
{
public:
	GAFactoryWrapperLensGAGenomeCalculator(unique_ptr<GAFactoryHelper<T,P1>> helper)
		: m_helperFactory(move(helper))
	{
		m_helperFactory->setParent(this);
	}

	bool_t init(const LensInversionParametersBase &params) override
	{
		// TODO: check parameter type? (P1)
		// Convert the parameters to the right class
		serut::VectorSerializer ser;
		if (!params.write(ser))
			return "Error serializing parameters: " + params.getErrorString();
		
		if (!m_helperFactory->init(&params))
			return "Can't init helper factory: " + m_helperFactory->getErrorString();

		return true;
	}

	bool_t createLens(const LensGAGenome &genome, std::unique_ptr<GravitationalLens> &lens) const override
	{
		string errStr = "unknown error";

		lens = move(m_helperFactory->createLens(genome.m_weights, genome.m_sheets, genome.m_scaleFactor, errStr));
		if (!lens.get())
			return errStr;
		return true;
	}

	size_t getNumberOfObjectives() const override { return m_helperFactory->getNumberOfFitnessComponents(); }
	bool allowNegativeValues() const override { return m_helperFactory->allowNegativeValues(); }
	size_t getNumberOfBasisFunctions() const override { return m_helperFactory->getNumberOfBasisFunctions(); }
	size_t getNumberOfSheets() const override { return m_helperFactory->getNumberOfSheets(); }
	size_t getMaximumNumberOfGenerations() const override { return m_helperFactory->getMaximumNumberOfGenerations(); }
	const LensFitnessObject &getLensFitnessObject() const override { return m_helperFactory->getFitnessObject(); }

	bool_t calculate(const eatk::Genome &genome, eatk::Fitness &fitness) override
	{
		const LensGAGenome &g = static_cast<const LensGAGenome&>(genome);
		LensGAFitness &f = static_cast<LensGAFitness&>(fitness);
		if (!m_helperFactory->calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data()))
			return "Unable to calculate fitness: " + m_helperFactory->getErrorString();
		return true;		
	}
private:
	unique_ptr<GAFactoryHelper<T,P1>> m_helperFactory;
};

template<class T, class P1>
class GeneralFactory : public LensGACalculatorFactory
{
public:
	std::unique_ptr<LensInversionParametersBase> createParametersInstance() override
	{
		return make_unique<P1>();
	}

	std::unique_ptr<LensGAGenomeCalculator> createCalculatorInstance(unique_ptr<LensFitnessObject> fitObj) override
	{
		auto helper = make_unique<GAFactoryHelper<T,P1>>(move(fitObj));
		auto wrapper = make_unique<GAFactoryWrapperLensGAGenomeCalculator<T,P1>>(move(helper));
		return wrapper;
	}
};

void registerWrapperCalculators()
{
	LensGACalculatorRegistry::instance().registerCalculatorFactory("singleplanecpu",
		make_unique<GeneralFactory<LensInversionGAFactorySinglePlaneCPU,LensInversionParametersSinglePlaneCPU>>());

	LensGACalculatorRegistry::instance().registerCalculatorFactory("multiplanegpu",
		make_unique<GeneralFactory<LensInversionGAFactoryMultiPlaneGPU,LensInversionParametersMultiPlaneGPU>>());
}

}
