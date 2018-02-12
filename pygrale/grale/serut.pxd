from libcpp.vector cimport vector

cdef extern from "serut/serializationinterface.h" namespace "serut":

    cdef cppclass SerializationInterface:
        pass

cdef extern from "serut/dummyserializer.h" namespace "serut":

    cdef cppclass DummySerializer(SerializationInterface):
        size_t getBytesWritten()

cdef extern from "serut/memoryserializer.h" namespace "serut":

    cdef cppclass MemorySerializer(SerializationInterface):
        MemorySerializer(const void *pReadBuffer, size_t readSize, void *pWriteBuffer, size_t writeSize)

cdef extern from "serut/vectorserializer.h" namespace "serut":
    cdef cppclass VectorSerializer(SerializationInterface):
        VectorSerializer()
        const vector[unsigned char] &getBuffer()
        unsigned char *getBufferPointer()
        int getBufferSize()

