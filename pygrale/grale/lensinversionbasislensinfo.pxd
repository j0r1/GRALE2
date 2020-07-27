from libcpp.memory cimport shared_ptr

cimport grale.vector2d as vector2d
cimport grale.gravitationallens as gravitationallens

cdef extern from "grale/lensinversionbasislensinfo.h" namespace "grale":
    cdef cppclass LensInversionBasisLensInfo:
        LensInversionBasisLensInfo(shared_ptr[gravitationallens.GravitationalLens] &lens,
                      vector2d.Vector2Dd center,
                      double relevantLensingMass)
