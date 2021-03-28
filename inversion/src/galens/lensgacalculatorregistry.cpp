#include "lensgacalculatorregistry.h"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace errut;

namespace grale
{

std::unique_ptr<LensGACalculatorRegistry> LensGACalculatorRegistry::s_instance;

LensGACalculatorRegistry &LensGACalculatorRegistry::instance()
{
    if (s_instance.get())
        return *s_instance;

    s_instance = std::move(unique_ptr<LensGACalculatorRegistry>(new LensGACalculatorRegistry()));

    // TODO: register defaults here?

    return *s_instance;
}

bool_t LensGACalculatorRegistry::registerCalculatorFactory(const std::string &name, unique_ptr<LensGACalculatorFactory> factory)
{
    // TODO
    return true;
}


LensGACalculatorFactory *LensGACalculatorRegistry::getFactory(const std::string &name)
{
    // TODO
    return nullptr;
}

LensGACalculatorRegistry::LensGACalculatorRegistry()
{
    if (s_instance.get()) // Shouldn't happen as constructor is private
    {
        cerr << "FATAL: one one instance of LensGACalculatorRegistry may exist!" << endl;
        exit(-1);
    }
}

LensGACalculatorRegistry::~LensGACalculatorRegistry()
{
}

}