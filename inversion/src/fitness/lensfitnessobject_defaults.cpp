#include "lensfitnessgeneral.h"

namespace grale
{

void registerDefaultLensFitnessObjects()
{
	LensFitnessObjectRegistry::instance().registerLensFitnessObject<LensFitnessGeneral>("general");
}

}
