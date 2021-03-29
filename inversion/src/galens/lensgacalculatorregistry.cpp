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
    if (!factory.get())
        return "Attempting to register null as factory";

    auto it = m_registry.find(name);
    if (it != m_registry.end())
        return "Name already in use";
    m_registry[name] = move(factory);
    // cout << "Registered factory name " << name << " in " << (void*)this << endl;
    return true;
}


LensGACalculatorFactory *LensGACalculatorRegistry::getFactory(const std::string &name)
{
    // cout << "Querying " << (void*)this << " for name " << name << endl;
    // for (auto &it : m_registry)
    //     cout << "Known = " << it.first << endl;
    auto it = m_registry.find(name);
    if (it == m_registry.end())
        return nullptr;

    return it->second.get();
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