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
        TypedParameter(bool v, int n=1)
        TypedParameter(int v, int n=1)
        TypedParameter(double v, int n=1)
        TypedParameter(string &v, int n=1)
        TypedParameter(const vector[bool] &v)
        TypedParameter(const vector[int] &v)
        TypedParameter(const vector[double] &v)
        TypedParameter(const vector[string] &v)

        void dump()

        Type getType()
        bool isArray()

        bool isEmpty()
        bool isBoolean()
        bool isInteger()
        bool isReal()
        bool isString()
        int getNumberOfEntries()

        bool getBooleanValue()
        int getIntegerValue()
        double getRealValue()
        string getStringValue()

        const vector[bool] &getBooleanValues()
        const vector[int] &getIntegerValues()
        const vector[double] &getRealValues()
        const vector[string] &getStringValues()

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
        void setParameterEmpty(const string &key)
        void setParameter(const string &key, bool v)
        void setParameter(const string &key, int v)
        void setParameter(const string &key, double v)
        void setParameter(const string &key, const string &v)
        void setParameter(const string &key, const vector[bool] &v)
        void setParameter(const string &key, const vector[int] &v)
        void setParameter(const string &key, const vector[double] &v)
        void setParameter(const string &key, const vector[string] &v)

        bool hasParameter(const string &key)
        const TypedParameter *getParameter(const string &key)
        bool getParameter(const string &key, TypedParameter &dst)

        bool getParameter(const string &key, bool &v)
        bool getParameter(const string &key, int &v)
        bool getParameter(const string &key, double &v)
        bool getParameter(const string &key, string &v)
        bool getParameter(const string &key, vector[bool] &v)
        bool getParameter(const string &key, vector[int] &v)
        bool getParameter(const string &key, vector[double] &v)
        bool getParameter(const string &key, vector[string] &v)

        void clearRetrievalMarkers()
        void getUnretrievedKeys(vector[string] &keys)

