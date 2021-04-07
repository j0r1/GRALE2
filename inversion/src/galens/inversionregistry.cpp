#include "lensfitnessobject.h"
#include "lensgaoldfactorywrapper.h"

namespace grale
{

void registerDefaultInversionComponents()
{
	registerDefaultLensFitnessObjects();
	registerWrapperCalculators();
}

}

