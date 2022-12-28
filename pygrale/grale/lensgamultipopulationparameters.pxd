from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport grale.serut as serut

cdef extern from "grale/lensgamultipopulationparameters.h" namespace "grale":
    cdef cppclass LensGAMultiPopulationParameters:
        LensGAMultiPopulationParameters()

        void setNumberOfPopulations(size_t s)
        size_t getNumberOfPopulations() const

        void setNumberOfInitialGenerationsToSkip(size_t n)
        size_t getNumberOfInitialGenerationsToSkip() const

        void setMigrationGenerationFraction(double x)
        double getMigrationGenerationFraction() const

        void setNumberOfIndividualsToLeavePopulation(size_t n)
        size_t getNumberOfIndividualsToLeavePopulation() const

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si) const

        string getErrorString()

