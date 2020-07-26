#include <grale/lensinversiongafactorymultiplanegpu.h>
#include "general_template.h"

extern "C"
{
	GAMODS_EXPORT mogal::GAFactory *CreateFactoryInstance()
	{
		return new grale::LensInversionGAFactory_General<grale::LensInversionGAFactoryMultiPlaneGPU>();
	}

	GAMODS_EXPORT grale::LensFitnessObject *CreateFitnessObject()
	{
		return new grale::LensFitnessGeneral();
	}
}

