#pragma once

#include "graleconfig.h"
#include "cosmology.h"
#include "lensinversionbasislensinfo.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include "scalesearchparameters.h"
#include <errut/errorbase.h>
#include <serut/serializationinterface.h>

namespace grale
{

class GRALE_IMPORTEXPORT LensInversionParametersMultiPlaneGPU : public errut::ErrorBase
{
public:
    LensInversionParametersMultiPlaneGPU();
    LensInversionParametersMultiPlaneGPU(const Cosmology &cosmology,
        const std::vector<double> &lensRedshifts,
        const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &basisLenses,
        const std::vector<std::shared_ptr<ImagesDataExtended>> &sourceImages,
        double massEstimate,
        bool useMassSheets,
        const ConfigurationParameters *pFitnessObjectParameters,
        int maxGenerations,
        bool allowNegativeWeights,
        const ScaleSearchParameters &massScaleSearchParams);
    ~LensInversionParametersMultiPlaneGPU();

   	const Cosmology &getCosmology() const                                           { return m_cosmology; }
	const std::vector<double> &getLensRedshifts() const                             { return m_lensRedshifts; }
	const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &getBasisLenses() const { return m_basisLenses; }
   	const std::vector<std::shared_ptr<ImagesDataExtended>> &getSourceImages() const { return m_images; }
    double getMassEstimate() const                                                  { return m_massEstimate; }
    bool useMassSheetBasisFunctions() const                                         { return m_useSheets; }
    const ConfigurationParameters *getFitnessObjectParameters() const				{ return m_fitnessObjectParams.get(); }
    int getMaximumNumberOfGenerations() const                                       { return m_maxGenerations; }
    bool getAllowNegativeWeights() const                                            { return m_allowNeg; }
    const ScaleSearchParameters &getMassScaleSearchParameters() const               { return m_scaleSearchParams; }

 	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
    Cosmology m_cosmology;
    std::vector<double> m_lensRedshifts;
    std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> m_basisLenses;
    // should have a "z" property, Ds will be set later based on that, Dds will be set to 0
    std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
    double m_massEstimate;
    bool m_useSheets;
    std::shared_ptr<ConfigurationParameters> m_fitnessObjectParams;
    int m_maxGenerations;
    bool m_allowNeg;
    ScaleSearchParameters m_scaleSearchParams;
};

} // end namespace