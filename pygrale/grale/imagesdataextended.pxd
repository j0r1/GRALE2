from libcpp.string cimport string
from libcpp cimport bool as cbool
from libcpp.vector cimport vector

cimport grale.imagesdata as imagesdata
cimport grale.configurationparameters as configurationparameters
cimport grale.serut as serut

cdef extern from "grale/imagesdataextended.h" namespace "grale":
    cdef cppclass ImagesDataExtended(imagesdata.ImagesData):
        ImagesDataExtended()
        ImagesDataExtended(const ImagesDataExtended &dat)
        ImagesDataExtended(const imagesdata.ImagesData &dat)

        cbool read(serut.SerializationInterface &si)
        cbool write(serut.SerializationInterface &si) const

        void setDs(double v)
        void setDds(double v)
        
        size_t getNumberOfExtraParameters() const
        void getAllExtraParameterKeys(vector[string] &keys)

        void clearExtraParameters()
        void setExtraParameter(string key, cbool v)
        void setExtraParameter(string key, int v)
        void setExtraParameter(string key, double v)
        void setExtraParameter(string key, const string v)

        cbool getExtraParameter(string key, configurationparameters.TypedParameter &dst) const

