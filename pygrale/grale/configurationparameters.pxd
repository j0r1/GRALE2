from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cimport grale.serut as serut
cimport grale.errut as errut

cdef extern from "grale/configurationparameters.h" namespace "grale::TypedParameter":
    cdef enum Type:
        pass

cdef extern from "grale/configurationparameters.h" namespace "grale":
    cdef cppclass TypedParameter(errut.ErrorBase):
        TypedParameter(const TypedParameter &src)
        TypedParameter()
        TypedParameter(bool v)
        TypedParameter(int v)
        TypedParameter(double v)
        TypedParameter(string &v)

        Type getType()
        bool isEmpty()
        bool isBoolean()
        bool isInteger()
        bool isReal()
        bool isString()

        bool getBooleanValue()
        int getIntegerValue()
        double getRealValue()
        string getStringValue()

        const TypedParameter &operator=(const TypedParameter &src)

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si)

    cdef cppclass ConfigurationParameters(errut.ErrorBase):
        ConfigurationParameters()
        ConfigurationParameters(const ConfigurationParameters &cfg)

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si)

        const ConfigurationParameters &operator=(const ConfigurationParameters &cfg)

        size_t getNumberOfParameters()
        void getAllKeys(vector[string] &keys)

        void clearParameters()
        void setParameter(const string &key, bool v)
        void setParameter(const string &key, int v)
        void setParameter(const string &key, double v)
        void setParameter(const string &key, const string &v)

        bool hasParameter(const string &key)
        bool getParameter(const string &key, TypedParameter &dst)
        bool getParameter(const string &key, bool &v)
        bool getParameter(const string &key, int &v)
        bool getParameter(const string &key, double &v)
        bool getParameter(const string &key, string &v)

        void clearRetrievalMarkers()
        void getUnretrievedKeys(vector[string] &keys)

