#include "lensfitnessobject.h"

using namespace std;
using namespace errut;

namespace grale
{

unique_ptr<LensFitnessObjectRegistry> LensFitnessObjectRegistry::s_instance;

LensFitnessObjectRegistry &LensFitnessObjectRegistry::instance()
{
    if (s_instance.get())
        return *s_instance;
    
    s_instance = unique_ptr<LensFitnessObjectRegistry>(new LensFitnessObjectRegistry());

    // TODO: register defaults?

    return *s_instance;
}

LensFitnessObjectRegistry::LensFitnessObjectRegistry()
{
    if (s_instance.get()) // Shouldn't happen as constructor is private
    {
        cerr << "FATAL: one one instance of LensFitnessObjectRegistry may exist!" << endl;
        exit(-1);
    }
}

LensFitnessObjectRegistry::~LensFitnessObjectRegistry()
{
}

bool_t LensFitnessObjectRegistry::registerFitnessObjectFactory(const std::string &name, std::unique_ptr<LensFitnessObjectFactory> fitnessObjectFactory)
{
    auto it = m_registry.find(name);
    if (it != m_registry.end())
        return "Name already in use";
    m_registry[name] = move(fitnessObjectFactory);
    return true;
}

std::unique_ptr<LensFitnessObject> LensFitnessObjectRegistry::createFitnessObject(const std::string &name)
{
    auto it = m_registry.find(name);
    if (it == m_registry.end())
        return nullptr;

    return it->second->createFitnessObject();
}

}
