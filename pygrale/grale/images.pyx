"""This module contains several classes that are related to the images
in a gravitational lensing system.

The :class:`ImagesData` class can be used to store properties of point-like
or extended images, and is used to pass data to a lens inversion algorithm.

The :class:`LensPlane` class calculates the deflection angles in the lens plane
of a specific gravitational lens. This mapping can subsequently be used to
calculate an image plane, represented by the :class:`ImagePlane` class. Such
an image plane is defined for a specific source, and hence requires you to
specify the angular diameter distance of the source. 

Once calculated, it can be used to obtain the critical lines and caustics for
the lens and source distance involved. The image plane instance can also be
used to render what the images of a specific source shape look like. The
source shapes are derived from the class :class:`SourceImage` and can be

 - :class:`CircularSource`: describes a source that has the shape of a disc
 - :class:`EllipticalSource`: describes a source that has an elliptical shape
 - :class:`PolygonSource`: describes a source which can be drawn as a polygon
 - :class:`DiscreteSource`: describes a source that consists of a grid of pixels
 - :class:`PolygonSource`: describes a point source
"""

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cython.operator cimport dereference as deref
import cython
import random
import numpy as np
cimport numpy as np
from cpython.array cimport array,clone
from . import lenses
from . import privutil

include "stringwrappers.pyx"

cimport grale.imagesdata as imagesdata
cimport grale.imagesdataextended as imagesdataextended
cimport grale.vector2d as vector2d
cimport grale.triangleindices as triangleindices
cimport grale.lensplane as lensplane
cimport grale.imageplane as imageplane
cimport grale.sourceimage as sourceimage
cimport grale.polygon2d as polygon2d
cimport grale.serut as serut
cimport grale.configurationparameters as configurationparameters

ctypedef triangleindices.TriangleIndices TriangleIndices
ctypedef vector2d.Vector2Dd Vector2Dd

class ImagesDataException(Exception):
    """This exception is raised if something goes wrong in the :class:`ImagesData` class."""
    pass

_imagesDataRndId = -int(random.random()*0xFFFFFF)

cdef class ImagesData:
    """This class is used as a container for the images of a single source. It's
    possible to specify a number of images, several points within each image, 
    and a triangulation of the points of an image. You can also add intensity
    information, time delay information and data about the shear."""
    
    cdef imagesdata.ImagesData *m_pImgData

    @staticmethod
    cdef imagesdata.ImagesData * _getImagesData(ImagesData img):
        return img.m_pImgData
    
    def __cinit__(self):
        self.m_pImgData = new imagesdata.ImagesData()

    def __init__(self, int numImages, cbool intensities = False, cbool shearInfo = False):
        """
        __init__(numImages, intensities = False, shearInfo = False)

        Parameters:

         - ``numImages``: the number of images this instance will contain. Can still be increased
           using the :func:`addImage` member function.
         - ``intensities``: flag indicating if intensity information will be stored.
         - ``shearInfo``: flag indicating if shear info will be stored.
        """
        
        if numImages == _imagesDataRndId: # use this special marked to create an uninitialized instance
            return

        # Not a special empty instance, initialize it
        if not self.m_pImgData.create(numImages, intensities, shearInfo):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

    def __dealloc__(self):
        del self.m_pImgData

    def addImage(self):
        """addImage()

        Add an image, and return the index that needs to be used to refer to this image."""

        cdef imageNum = self.m_pImgData.addImage()
        if imageNum < 0:
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))
        return imageNum

    def addPoint(self, int imageNum, point, intensity = None, shear = None):
        """addPoint(imageNum, point, intensity = None, shear = None)

        Add a point to an image.

        Parameters:

         - ``imageNum``: the index of the image to which the point should be added
         - ``point``: the 2D angular coordinates of the point
         - ``intensity``: if intensity information was specified in the constructor,
           then this intensity will be specified for the point
         - ``shear``: if shear info was specified in the constructor, then this 2D
           shear info will be stored for the point.
        """

        cdef pointNum = 0
        cdef vector2d.Vector2Dd p = vector2d.Vector2Dd(point[0], point[1])
        if intensity is None:
            if shear is None:
                pointNum = self.m_pImgData.addPoint(imageNum, p)
            else:
                pointNum = self.m_pImgData.addPoint(imageNum, p, shear[0], shear[1])
        else:
            if shear is None:
                pointNum = self.m_pImgData.addPoint(imageNum, p, intensity)
            else:
                pointNum = self.m_pImgData.addPoint(imageNum, p, intensity, shear[0], shear[1])

        if pointNum < 0:
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))
        return pointNum

    def addGroup(self):
        """addGroup()

        To specify which image points in different images correspond to each other, you
        can use an image group. This function creats a new group, and returns the identifier
        for the new group.
        """
        cdef groupNum = self.m_pImgData.addGroup()
        if groupNum < 0:
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))
        return groupNum

    def addGroupPoint(self, int groupNumber, int imageIndex, int pointIndex):
        """addGroupPoint(self, groupNumber, imageIndex, pointIndex)

        Add the specified point to a certain point group.

        Parameters:

         - ``groupNumber``: the ID of the group to add the point to
         - ``imageIndex``: specifies the ID of the image the point refers to
         - ``pointIndex``: the ID of the point itself
        """

        if not self.m_pImgData.addGroupPoint(groupNumber, imageIndex, pointIndex):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

    def addTimeDelayInfo(self, int imageIndex, int pointIndex, double timeDelay):
        """addTimeDelayInfo(imageIndex, pointIndex, timeDelay)

        Specifies that the time delay ``timeDelay`` should be associated with the
        point with image ID ``imageIndex`` and point ID ``pointIndex``.
        """
        if not self.m_pImgData.addTimeDelayInfo(imageIndex, pointIndex, timeDelay):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

    def addTriangle(self, int imageNumber, int index1, int index2, int index3):
        """addTriangle(imageNumber, index1, index2, index3)

        This function allows you to customize the triangulation of an image.
        It says that a triangle should be defined within the image with ID 
        ``imageNumber``, and that its three points are specified by
        point ID ``index1``, ``index2`` and  ``index3``.
        """
        if not self.m_pImgData.addTriangle(imageNumber, index1, index2, index3):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

    @staticmethod
    def load(fileName):
        """load(fileName)

        This function attempts to interpret the file with name ``fileName`` as
        an images data set, and returns the loaded instance if successful.
        """

        img = ImagesData(_imagesDataRndId) # Make sure an empty instance is allocated

        # Try to load the data
        if not img.m_pImgData.load(B(fileName)):
            raise ImagesDataException(S(img.m_pImgData.getErrorString()))

        return img

    def save(self, fileName):
        """save(fileName)

        Saves the current images data set to the file with name ``fileName``.
        """
        if not self.m_pImgData.save(B(fileName)):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

    def hasIntensities(self):
        """hasIntensities()

        Returns a boolean indicating if intensity information is present.
        """
        return self.m_pImgData.hasIntensities()

    def hasShearInfo(self):
        """hasShearInfo()

        Returns a boolean indicating if shear information is present.
        """
        return self.m_pImgData.hasShearInfo()

    def hasTimeDelays(self):
        """hasTimeDelays()

        Returns a boolean indicating if time delay information is present.
        """
        return self.m_pImgData.hasTimeDelays()
    
    def getNumberOfImages(self):
        """getNumberOfImages()

        Returns the number of images that are contained in this images data set.
        If this function returns ``Ni``, valid image IDs will range from 0 to ``Ni``-1.
        """
        return self.m_pImgData.getNumberOfImages()

    cdef _checkImageNumber(self, int i):
        if i < 0 or i >= self.m_pImgData.getNumberOfImages():
            raise ImagesDataException("Invalid image number")

    cdef _checkImagePointNumber(self, int image, int point):
        self._checkImageNumber(image)
        if point < 0 or point >= self.m_pImgData.getNumberOfImagePoints(image):
            raise ImagesDataException("Invalid image point number")

    cdef _checkIntensities(self):
        if not self.hasIntensities():
            raise ImagesDataException("No intensities are stored in this data set");

    cdef _checkShear(self):
        if not self.hasShearInfo():
            raise ImagesDataException("No shear data is stored in this data set");

    cdef _checkGroupNumber(self, int group):
        if group < 0 or group >= self.m_pImgData.getNumberOfGroups():
            raise ImagesDataException("Invalid group number")

    cdef _checkGroupPoint(self, int group, int pointnr):
        cdef int np = 0
        self._checkGroupNumber(group)
        np = self.m_pImgData.getNumberOfGroupPoints(group)
        if pointnr < 0 or pointnr >= np:
            raise ImagesDataException("Invalid group point number")

    cdef _checkTimeDelayIndex(self, int index):
        if index < 0 or index >= self.m_pImgData.getNumberOfTimeDelays():
            raise ImagesDataException("Invalid time delay point index")

    def getNumberOfImagePoints(self, i):
        """getNumberOfImagePoints(i)

        Returns the number of image points are stored for the image with ID ``i``.
        If this function returns ``Np``, then valid point IDs will range from 0
        to ``Np``-1.
        """
        self._checkImageNumber(i)
        return self.m_pImgData.getNumberOfImagePoints(i)

    def getImagePointPosition(self, image, point):
        """getImagePointPosition(image, point)

        Returns the position of the point with ID ``point`` in image with ID ``image``.
        """
        self._checkImagePointNumber(image, point)
        pos = self.m_pImgData.getImagePointPosition(image, point)
        return np.array([pos.getX(), pos.getY()])

    def getAllImagePoints(self):
        """getAllImagePoints()

        A convenience function that calls calls :func:`getImagePoints <grale.images.ImagesData.getImagePoints>`
        for each image, returning a list of the lists of dictionaries described in
        that function call. This means that the main list returned will be indexed
        by image number and each entry of the main list will itself be a list of
        points in an image. 
        """
        return [ list(self.getImagePoints(i)) for i in range(self.getNumberOfImages()) ]

    def getImagePoints(self, image):
        """getImagePoints(image)

        Returns an iterator that can be used to obtain all points inside
        the image with ID ``image``. For each point, an dictionary will be
        provided with the following entries:

         - ``position``: contains the 2D position of the point
         - ``intensity`` (optional): contains the intensity information associated
           to the point
         - ``shear`` (optional): contains the two shear components associated to
           the point
        """
        numPoints = self.getNumberOfImagePoints(image)
        intens = self.hasIntensities()
        shear = self.hasShearInfo()
        for i in range(numPoints):
            obj = { }
            obj["position"] = self.getImagePointPosition(image, i)
            if intens:
                obj["intensity"] = self.getImagePointIntensity(image, i)
            if shear:
                obj["shear"] = self.getShearComponents(image, i)

            yield obj

    def getImagePointIntensity(self, image, point):
        """getImagePointIntensity(image, point)

        Returns the intensity information stored for the point with ID
        ``point`` inside the image with ID ``image``.
        """
        self._checkImagePointNumber(image, point)
        self._checkIntensities()
        return self.m_pImgData.getImagePointIntensity(image, point)

    def getShearComponent1(self, image, point):
        """getShearComponent1(image, point)

        Returns the first shear component stored for the point with ID
        ``point`` inside the image with ID ``image``.
        """
        self._checkImagePointNumber(image, point)
        self._checkShear()
        return self.m_pImgData.getShearComponent1(image, point)

    def getShearComponent2(self, image, point):
        """getShearComponent2(image, point)

        Returns the second shear component stored for the point with ID
        ``point`` inside the image with ID ``image``.
        """
        self._checkImagePointNumber(image, point)
        self._checkShear()
        return self.m_pImgData.getShearComponent2(image, point)

    def getShearComponents(self, image, point):
        """getShearComponents(image, point)

        Returns the both shear component stored for the point with ID
        ``point`` inside the image with ID ``image``.
        """
        self._checkImagePointNumber(image, point)
        self._checkShear()
        return np.array([self.m_pImgData.getShearComponent1(image, point), self.m_pImgData.getShearComponent2(image, point)])

    def getNumberOfGroups(self):
        """getNumberOfGroups()

        Returns the number of point groups stored in this object. A point group
        is a set of corresponding points in several images. If ``Ng`` is returned,
        valid point group IDs range from 0 to ``Ng``-1.
        """
        return self.m_pImgData.getNumberOfGroups()

    def getNumberOfGroupPoints(self, group):
        """getNumberOfGroupPoints(group)

        Returnes the number of points in the group with ID ``group``. If ``Ngp``
        is returned, valid group point indices range from 0 to ``Ngp``-1.
        """
        self._checkGroupNumber(group)
        return self.m_pImgData.getNumberOfGroupPoints(group)

    def getGroupPointIndices(self, group, pointnr):
        """getGroupPointIndices(group, pointnr)

        For group with ID ``group`` and point ID ``pointnr`` within this group,
        return the ``(img, point)`` tuple containing the image ID ``img`` and
        point ID ``point`` that the group point refers to.
        """
        cdef int img = 0
        cdef int point = 0
        self._checkGroupPoint(group, pointnr)
        self.m_pImgData.getGroupPointIndices(group, pointnr, cython.address(img), cython.address(point))
        return (img, point)

    def getNumberOfTimeDelays(self):
        """getNumberOfTimeDelays()

        Returns the number of time delays stored in this instance. If ``Nt`` is the value
        returned, valid time delay IDs range from 0 to ``Nt``-1.
        """
        return self.m_pImgData.getNumberOfTimeDelays()

    def getTimeDelay(self, index):
        """getTimeDelay(index)

        For the time delay with ID ``index``, return the ``(img, pointnr, delay)`` tuple
        that describes the point and the time delay associated with the point. The point
        is specified by giving the image ID ``img`` and point ID ``pointnr`` within that
        image.
        """
        cdef int img = 0
        cdef int pointnr = 0
        cdef double delay = 0
        self._checkTimeDelayIndex(index)
        self.m_pImgData.getTimeDelay(index, cython.address(img), cython.address(pointnr), cython.address(delay))
        return (img, pointnr, delay)

    def hasTriangulation(self):
        """hasTriangulation()

        Returns a boolean indicating if this instance contains a triangulation of points
        within one or more images.
        """
        return self.m_pImgData.hasTriangulation()

    def getTriangles(self, int image):
        """getTriangles(image)

        Returns the triangles for the triangulation that's stored for the image with ID ``image``,
        as a list of tuples, each containing three point indices within this image.
        """
        cdef vector[TriangleIndices] triangles

        if not self.m_pImgData.getTriangles(image, triangles):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

        l = [ ]
        for t in triangles:
            l.append( (t.getIndex(0), t.getIndex(1), t.getIndex(2)) )

        return l

    def clearTriangulation(self):
        """clearTriangulation()

        Clears all triangulations in this images data instance.
        """
        self.m_pImgData.clearTriangulation()

    def getTopRightCorner(self):
        """getTopRightCorner()

        Finds out what the rectangular area surrounded by all image points is, and returns
        the top right corner.
        """
        cdef vector2d.Vector2Dd corner = self.m_pImgData.getTopRightCorner()
        return np.array([corner.getX(), corner.getY()])

    def getBottomLeftCorner(self):
        """getBottomLeftCorner()

        Finds out what the rectangular area surrounded by all image points is, and returns
        the bottom left corner.
        """
        cdef vector2d.Vector2Dd corner = self.m_pImgData.getBottomLeftCorner()
        return np.array([corner.getX(), corner.getY()])

    def centerOnPosition(self, double ra, double dec):
        """centerOnPosition(ra, dec)

        Assuming that all image point coordinates were added using right ascention as
        the first coordinate and declination as the second, this function recalculates
        all coordinates so that they are now specified relative to ``ra`` and ``dec``.
        """
        self.m_pImgData.centerOnPosition(ra, dec)

    def uncenterOnPosition(self, double ra, double dec):
        """uncenterOnPosition(ra, dec)

        Performs the inverse operation of :func:`centerOnPosition <grale.images.ImagesData.centerOnPosition>`
        """
        self.m_pImgData.uncenterOnPosition(ra, dec)

    def subtractIntensity(self, double v):
        """subtractIntensity()

        For all intensities stored in this images data instance, the value ``v`` will
        be subtracted.
        """
        self.m_pImgData.subtractIntensity(v)

# For internal use
cdef class ImagesDataExtended:
    cdef imagesdataextended.ImagesDataExtended *m_pImgDataExt

    def __cinit__(self):
        self.m_pImgDataExt = NULL

    cdef _check(self):
        if self.m_pImgDataExt == NULL:
            raise ImagesDataException("No internal extended images data instance has been set")

    def __init__(self, ImagesData img):
        cdef imagesdata.ImagesData *pImgDat
        pImgDat = ImagesData._getImagesData(img)
        self.m_pImgDataExt = new imagesdataextended.ImagesDataExtended(deref(pImgDat))

    def __dealloc__(self):
        del self.m_pImgDataExt

    def setDs(self, Ds):
        self._check()
        self.m_pImgDataExt.setDs(Ds)

    def setDds(self, Dds):
        self._check()
        self.m_pImgDataExt.setDds(Dds)

    def setExtraParameter(self, key, value):
        cdef string bKey
        cdef cbool bvalue
        cdef int ivalue
        cdef double dvalue
        cdef string svalue
        
        self._check()

        bKey = B(key)
        bvalue = False
        ivalue = 0
        dvalue = 0.0
        svalue = ""

        if type(value) == bool:
            bvalue = value
            self.m_pImgDataExt.setExtraParameter(bKey, bvalue)
        elif type(value) == int:
            ivalue = value
            self.m_pImgDataExt.setExtraParameter(bKey, ivalue)
        elif type(value) == float:
            dvalue = value
            self.m_pImgDataExt.setExtraParameter(bKey, dvalue)
        else: # assume it's a string
            svalue = B(value)
            self.m_pImgDataExt.setExtraParameter(bKey, svalue)

    def getExtraParameters(self):
        cdef vector[string] keys;
        cdef configurationparameters.TypedParameter param
        
        retObj = { }

        self._check()
        self.m_pImgDataExt.getAllExtraParameterKeys(keys)
        for key in keys:
            if not self.m_pImgDataExt.getExtraParameter(key, param):
                raise ImagesDataException(S(self.m_pImgDataExt.getErrorString()))

            strKey = S(key)
            if param.isBoolean():
                retObj[strKey] = param.getBooleanValue()
            elif param.isInteger():
                retObj[strKey] = param.getIntegerValue()
            elif param.isReal():
                retObj[strKey] = param.getRealValue()
            elif param.isString():
                retObj[strKey] = S(param.getStringValue())
            else:
                retObj[strKey] = None

        return retObj
    
    def clearExtraParameters(self):
        self._check()
        self.m_pImgDataExt.clearExtraParameters()

    @staticmethod
    def fromBytes(bytes b):
        cdef array[char] buf = chararrayfrombytes(b)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(b), NULL, 0)
        cdef imagesdataextended.ImagesDataExtended *pImgDat = NULL
        cdef string errorString
        
        try:
            pImgDat = new imagesdataextended.ImagesDataExtended()
            if not pImgDat.read(deref(m)):
                errorString = pImgDat.getErrorString()
                del pImgDat
                raise ImagesDataException(S(errorString))
        finally:
            del m

        imgDat = ImagesDataExtended(ImagesData(1))
        del imgDat.m_pImgDataExt
        imgDat.m_pImgDataExt = pImgDat
        return imgDat

    def toBytes(self):
        cdef serut.VectorSerializer vSer

        self._check()
        if not self.m_pImgDataExt.write(vSer):
            raise ImagesDataException(S(self.m_pImgDataExt.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

def centerOnPosition(position, centerRaDec):
    """centerOnPosition(position, centerRaDec)
    
    Recalculates ``position`` as relative to ``centerRaDec``"""
    # For now, we'll just abuse the method in ImagesData

    img = ImagesData(1)
    img.addPoint(0, position)
    img.centerOnPosition(centerRaDec[0], centerRaDec[1])
    return img.getImagePointPosition(0,0)

def uncenterOnPosition(position, centerRaDec):
    """uncenterOnPosition(position, centerRaDec)

    Performs the opposite calculation of :func:`centerOnPosition <grale.images.centerOnPosition>`."""
    # For now, we'll just abuse the method in ImagesData

    img = ImagesData(1)
    img.addPoint(0, position)
    img.uncenterOnPosition(centerRaDec[0], centerRaDec[1])
    return img.getImagePointPosition(0,0)

class LensPlaneException(Exception):
    """This exception is raised if something goes wrong in the :class:`LensPlane` class."""
    pass

cdef class LensPlane:
    cdef lensplane.LensPlane *m_pLensPlane
    cdef object feedback

    cdef _check(self):
        if self.m_pLensPlane == NULL:
            raise LensPlaneException("No internal lens plane has been set")

    def __cinit__(self):
        self.m_pLensPlane = NULL

    def __dealloc__(self):
        del self.m_pLensPlane

    def __init__(self, lens, bottomLeft, topRight, int numX, int numY, renderer = None, feedbackObject = "default"):
        """__init__(lens, bottomLeft, topRight, numX, numY, renderer = None, feedbackObject = "default")

        This creates a LensPlane instance that covers the area specified by the ``bottomLeft`` and
        ``topRight`` corners. Such a LensPlane instance calculates and stores the deflection angles
        for the :class:`grale.lenses.GravitationalLens` instance ``lens`` on a grid of ``numX`` points wide by
        ``numY`` points high.

        If ``renderer`` is ``None``, the mapping is calculated single threaded, within this
        Python process. Other renderers can be specified as well, for example to calculate the
        mapping faster using multiple cores with the MPI renderer. See the :mod:`grale.renderers`
        module for more information. 
        
        Feedback while rendering can be provided by specifying a ``feedbackObject`` parameter. See
        the :mod:`grale.feedback` module for more information about allowed values.
        """
        if lens is None and bottomLeft is None and topRight is None:
            return

        self.m_pLensPlane = new lensplane.PyLensPlane(self)

        cdef serut.MemorySerializer *m = NULL
        cdef serut.MemorySerializer *m2 = NULL
        cdef array[char] buf
        cdef array[char] buf2

        renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "LENSPLANE")
        self.feedback = feedbackObject

        if renderer and renderer.renderType != "LENSPLANE":
            raise LensPlaneException("Specified renderer does not appear to be for the LENSPLANE type")

        try:
            lensData = lens.toBytes()
            buf = chararrayfrombytes(lensData)
            m = new serut.MemorySerializer(buf.data.as_voidptr, len(lensData), NULL, 0)

            if renderer is None:
                if not self.m_pLensPlane.init(deref(m), Vector2Dd(bottomLeft[0], bottomLeft[1]), Vector2Dd(topRight[0], topRight[1]),
                                              numX, numY):
                    raise LensPlaneException("Unable to initialize lens plane: " + S(self.m_pLensPlane.getErrorString()))
            else:
                b = renderer.render(lensData, bottomLeft, topRight, numX, numY)
                buf2 = chararrayfrombytes(b)
                m2 = new serut.MemorySerializer(buf2.data.as_voidptr, len(b), NULL, 0)
                if not self.m_pLensPlane.init(deref(m), Vector2Dd(bottomLeft[0], bottomLeft[1]), Vector2Dd(topRight[0], topRight[1]),
                                              numX, numY, deref(m2)):
                    raise LensPlaneException("Unable to initialize lens plane: " + S(self.m_pLensPlane.getErrorString()))
        finally:
            del m
            del m2

    def getLens(self):
        """getLens()

        This returns a copy of the :class:`grale.lenses.GravitationalLens` instance that was
        used to create the deflections.
        """
        cdef vector[unsigned char] buf
        cdef string errStr

        self._check()

        if not lensplane.PyLensPlane.getLensBytes(self.m_pLensPlane, buf, errStr):
            raise LensPlaneException(S(errStr))

        return lenses.GravitationalLens.fromBytes(<bytes>(&buf[0])[:len(buf)])

    def createDeflectionGridLens(self):
        """createDeflectionGridLens()

        This creates a :class:`grale.lenses.GravitationalLens` instance that uses the calculated
        deflections. In between the grid points, the values are interpolated. Outside the specified
        region, the lens effect will not be correct.
        """
        cdef vector[unsigned char] buf
        cdef string errStr

        self._check()

        if not lensplane.PyLensPlane.createDeflectionGridLens(self.m_pLensPlane, buf, errStr):
            raise LensPlaneException(S(errStr))

        return lenses.GravitationalLens.fromBytes(<bytes>(&buf[0])[:len(buf)])
    
    def getRenderInfo(self):
        """getRenderInfo()

        Returns a dictionary with the following entries:

         - ``bottomleft``: the bottom-left corner that was specified in the constructor
           of this instance
         - ``topright``: the top-right corner that was specified in the constructor of
           this instance
         - ``xpoints``: the number of points in the x-direction, between the left and right
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled
         - ``ypoints``: the number of points in the y-direction, between the bottom and top
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled
        """
        cdef Vector2Dd bl, tr
        cdef int xpoints, ypoints, xpixels, ypixels

        self._check()

        bl = self.m_pLensPlane.getBottomLeft()
        tr = self.m_pLensPlane.getTopRight()
        xpoints = self.m_pLensPlane.getNumXPoints()
        ypoints = self.m_pLensPlane.getNumYPoints()

        obj = { }
        obj["bottomleft"] = np.array([bl.getX(), bl.getY()])
        obj["topright"] = np.array([tr.getX(), tr.getY()])
        obj["xpoints"] = xpoints
        obj["ypoints"] = ypoints
        return obj

    def getAlphas(self):
        """TODO:"""
        cdef np.ndarray[double,ndim=2] alphax, alphay
        cdef int xpoints, ypoints, x, y;
        cdef Vector2Dd alpha;

        self._check()

        xpoints = self.m_pLensPlane.getNumXPoints()
        ypoints = self.m_pLensPlane.getNumYPoints()
        
        alphax = np.zeros([ypoints, xpoints], dtype = np.double)
        alphay = np.zeros([ypoints, xpoints], dtype = np.double)
        for y in range(ypoints):
            for x in range(xpoints):
                alpha = self.m_pLensPlane.getAlpha(x, y)
                alphax[y,x] = alpha.getX()
                alphay[y,x] = alpha.getY()

        return { "alpha_x": alphax, "alpha_y": alphay }

    def getAlphaVectorDerivatives(self):
        """TODO:"""
        cdef np.ndarray[double,ndim=2] alphaxx, alphayy, alphaxy
        cdef int xpoints, ypoints, x, y;
        cdef double axx, ayy, axy

        self._check()

        axx = 0
        ayy = 0
        axy = 0
        xpoints = self.m_pLensPlane.getNumXPoints()
        ypoints = self.m_pLensPlane.getNumYPoints()
        alphaxx = np.zeros([ypoints, xpoints], dtype = np.double)
        alphayy = np.zeros([ypoints, xpoints], dtype = np.double)
        alphaxy = np.zeros([ypoints, xpoints], dtype = np.double)
        for y in range(ypoints):
            for x in range(xpoints):
                self.m_pLensPlane.getAlphaDerivatives(x, y, axx, ayy, axy)
                alphaxx[y,x] = axx
                alphayy[y,x] = ayy
                alphaxy[y,x] = axy

        return {  "alpha_xx": alphaxx, "alpha_yy": alphayy, "alpha_xy": alphaxy }

    def save(self, fileName):
        """save(fileName)

        Store this instance in a file called ``fileName``.
        """
        self._check()

        if not self.m_pLensPlane.save(B(fileName)):
            raise LensPlaneException(S(self.m_pLensPlane.getErrorString()))

    @staticmethod
    def load(fileName):
        """load(fileName)

        Attempts to load a LensPlane instance from the file called ``fileName``.
        If successful, this returns the newly loaded instance.
        """
        cdef string errorString
        cdef lensplane.LensPlane *pLensPlane = NULL

        if not lensplane.LensPlane.load(B(fileName), cython.address(pLensPlane), errorString):
            raise LensPlaneException(S(errorString))

        lensPlane = LensPlane(None, None, None, 0, 0)
        lensPlane.m_pLensPlane = pLensPlane

        return lensPlane

    def _onFeedbackStatus(self, s):
        try:
            if self.feedback:
                self.feedback.onStatus(s)
        except Exception as e:
            print("Warning: ignoring exception ({}) in onStatus".format(e))

    def _onFeedbackPercentage(self, x):
        try:
            if self.feedback:
                self.feedback.onProgress(x)
        except Exception as e:
            print("Warning: ignoring exception ({}) in onProgress".format(e))

class ImagePlaneException(Exception):
    """This exception is raised if something goes wrong in the :class:`ImagePlane` class."""
    pass

cdef class ImagePlane:
    cdef imageplane.ImagePlane *m_pImgPlane

    def __cinit__(self):
        self.m_pImgPlane = new imageplane.ImagePlane()

    def __dealloc__(self):
        del self.m_pImgPlane

    def __init__(self, LensPlane lensplane, double Ds, double Dds):
        """__init__(lensplane, Ds, Dds)

        Based on the :class:`LensPlane` instance in ``lensplane``, which contains the deflection
        field at a number of grid points, an ImagePlane instance is created for a certain source
        plane. This source plane has angular diameter distance ``Ds`` relative to the observer,
        and ``Dds`` relative to the lens itself.
        """
        if not self.m_pImgPlane.init(<imageplane.LensPlane*>lensplane.m_pLensPlane, Ds, Dds):
            raise ImagePlaneException("Unable to initialize image plane: " + S(self.m_pImgPlane.getErrorString()))

    def getDs(self):
        """getDs()

        Returns the ``Ds`` parameter that was specified in the constructor.
        """
        return self.m_pImgPlane.getDs()

    def getDds(self):
        """getDds()

        Returns the ``Dds`` parameter that was specified in the constructor.
        """
        return self.m_pImgPlane.getDds()

    def getCriticalLines(self):
        """getCriticalLines()

        This returns a list describing the critical lines associated with this
        image plane. Each entry in the list is itself a list of 2D points, describing
        a connected part of a critical line.
        """
        cdef vector[vector[imageplane.Vector2Dd]] segments = self.m_pImgPlane.getCriticalLineSegments()
        
        return ImagePlane._segsToList(segments)

    @staticmethod
    cdef _segsToList(vector[vector[imageplane.Vector2Dd]] &segments):
        segs = [ ]
        for line in segments:
            l = [ ]
            for point in line:
                l.append(np.array([ point.getX(), point.getY() ]))
            segs.append(l)

        return segs

    def getCaustics(self):
        """getCaustics()

        This returns a list describing the caustics associated with this
        image plane. Each entry in the list is itself a list of 2D points, describing
        a connected part of a caustic.
        """
        cdef vector[vector[imageplane.Vector2Dd]] segments = self.m_pImgPlane.getCausticSegments()
        
        return ImagePlane._segsToList(segments)

    def traceBeta(self, beta):
        """traceBeta(beta)

        Estimates the image plane positions to which the source plane position ``beta``
        corresponds. Returns a list of 2D points.
        """
        cdef vector[imageplane.Vector2Dd] thetas

        self.m_pImgPlane.traceBeta(imageplane.Vector2Dd(beta[0], beta[1]), thetas)
        return [ np.array([t.getX(), t.getY()]) for t in thetas ]

    def getRenderInfo(self):
        """getRenderInfo()

        Returns a dictionary with the following entries:

         - ``bottomleft``: the bottom-left corner that is relevant for this instance. This
           is taken from the LensPlane instance specified in the constructor.
         - ``topright``: the top-right corner that is relevant for this instance. This
           is taken from the LensPlane instance specified in the constructor.
         - ``xpoints``: the number of points in the x-direction, between the left and right
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled. This is taken from the LensPlane
           instance specified in the constructor.
         - ``ypoints``: the number of points in the y-direction, between the bottom and top
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled. This is taken from the LensPlane
           instance specified in the constructor.
         - ``xpixels``: based on the image plane to source plane mappings that are known at
           ``xpoints`` * ``ypoints`` grid points, a number of pixels can be defined
           that can contain light from the source plane, and these pixels will
           be used when rendering the image plane or the source plane with
           :func:`renderImages` or :func:`renderSources`. This ``xpixels`` value
           specifies the number of pixels in the x-direction, and is one less
           than ``xpoints``.
         - ``ypixels``: similar to ``xpixels``, but for the y-direction.
        """
        cdef Vector2Dd bl, tr
        cdef int xpoints, ypoints, xpixels, ypixels

        bl = self.m_pImgPlane.getBottomLeft()
        tr = self.m_pImgPlane.getTopRight()
        xpoints = self.m_pImgPlane.getNumXPoints()
        ypoints = self.m_pImgPlane.getNumYPoints()
        xpixels = self.m_pImgPlane.getNumXPixels()
        ypixels = self.m_pImgPlane.getNumYPixels()

        obj = { }
        obj["bottomleft"] = np.array([bl.getX(), bl.getY()])
        obj["topright"] = np.array([tr.getX(), tr.getY()])
        obj["xpoints"] = xpoints
        obj["ypoints"] = ypoints
        obj["xpixels"] = xpixels
        obj["ypixels"] = ypixels
        return obj

    def renderSources(self, sourceList, np.ndarray[double, ndim=2] plane = None, int subSamples = 9):
        """renderSources(sourceList, plane = None, subSamples = 9)

        For the list of :class:`SourceImage` derived classes in ``sourceList``, this function
        calculates what the sources look like based on the dimensions and number of pixels for
        this ImagePlane instance. This is what the image plane would look like if the
        gravitational lens effect could be turned off.

        The function returns a 2D NumPy array containing ``ypixels`` rows, each of ``xpixels``
        pixels wide (see also :func:`getRenderInfo`). If ``plane`` is specified, the results are
        stored in that 2D NumPy instance, which must have the same dimensions.

        Each pixel is sub-sampled ``sqrt(subSamples)`` times in x- and y- direction, to be able to
        roughly approximate the integration that's needed over the surface area of a pixel.

        Note that a call to only
        `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_.
        will plot the 0,0 value in ``plane`` (the bottom-left value) as the top-left corner
        causing the result to appear mirrored in the y-direction (the y-axis will point down).
        A subsequent call to `invert_yaxis <https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.invert_yaxis.html>`_
        might be useful.
        """
        cdef int numX = self.m_pImgPlane.getNumXPixels()
        cdef int numY = self.m_pImgPlane.getNumYPixels()
        cdef int x, y, storeY
        cdef vector[sourceimage.SourceImage *] sources
        cdef sourceimage.SourceImage *pSrc = NULL

        if plane is None:
            plane = np.zeros([numY, numX], dtype = np.double)

        if plane.shape[0] != numY or plane.shape[1] != numX:
            raise ImagePlaneException("Render plane must have {} rows and {} columns".format(numY, numX))

        for s in sourceList:
            pSrc = SourceImage._getSourceImage(s)
            if pSrc != NULL:
                sources.push_back(pSrc)

        for y in range(numY):
            for x in range(numX):
                plane[y,x] = self.m_pImgPlane.getSourceIntensityAccurate(sources, x, y, subSamples)

        return plane
    
    def renderImages(self, sourceList, np.ndarray[double, ndim=2] plane = None, int subSamples = 9):
        """renderImages(sourceList, plane = None, subSamples = 9)

        For the list of :class:`SourceImage` derived classes in ``sourceList``, this function
        calculates what the images look like based on the dimensions and number of pixels for
        this ImagePlane instance. 

        The function returns a 2D NumPy array containing ``ypixels`` rows, each of ``xpixels``
        pixels wide (see also :func:`getRenderInfo`). If ``plane`` is specified, the results are
        stored in that 2D NumPy instance, which must have the same dimensions.

        Each pixel is sub-sampled ``sqrt(subSamples)`` times in x- and y- direction, to be able to
        roughly approximate the integration that's needed over the surface area of a pixel.

        Note that a call to only
        `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_.
        will plot the 0,0 value in ``plane`` (the bottom-left value) as the top-left corner
        causing the result to appear mirrored in the y-direction (the y-axis will point down).
        A subsequent call to `invert_yaxis <https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.invert_yaxis.html>`_
        might be useful.
        """
        cdef int numX = self.m_pImgPlane.getNumXPixels()
        cdef int numY = self.m_pImgPlane.getNumYPixels()
        cdef int x, y, storeY
        cdef vector[sourceimage.SourceImage *] sources
        cdef sourceimage.SourceImage *pSrc = NULL

        if plane is None:
            plane = np.zeros([numY, numX], dtype = np.double)

        if plane.shape[0] != numY or plane.shape[1] != numX:
            raise ImagePlaneException("Render plane must have {} rows and {} columns".format(numY, numX))

        for s in sourceList:
            pSrc = SourceImage._getSourceImage(s)
            if pSrc != NULL:
                sources.push_back(pSrc)

        for y in range(numY):
            for x in range(numX):
                plane[y,x] = self.m_pImgPlane.getImageIntensityAccurate(sources, x, y, subSamples)

        return plane

    @staticmethod
    cdef _findSegmentAt(int x, int y, np.ndarray[np.uint8_t, ndim=2] markers):
        cdef int numX = markers.shape[1]
        cdef int numY = markers.shape[0]

        if markers[y, x] == 0:
            return [ ]

        coords = [ ]
        neighbours = [ (x,y) ]
        while len(neighbours) > 0:
            x, y = neighbours.pop()
            markers[y, x] = 0
            coords.append((x,y))

            if y > 0 and markers[y-1,x] > 0:
                neighbours.append( (x, y-1) )
            if y < numY-1 and markers[y+1,x] > 0:
                neighbours.append( (x, y+1) )
            if x > 0 and markers[y,x-1] > 0:
                neighbours.append( (x-1, y) )
            if x < numX-1 and markers[y,x+1] > 0:
                neighbours.append( (x+1, y) )

        return coords

    @staticmethod
    def static_segment(np.ndarray[double, ndim=2] plane, bottomLeft, topRight, threshold = 0.0):
        """static_segment(plane, bottomLeft, topRight, threshold = 0.0)

        Similar tp :func:`segment`, but with the bottom left and rop right corners
        set manually.
        """
        cdef int numX
        cdef int numY
        cdef int x, y
        cdef double pixw, pixh, X, Y
        cdef Vector2Dd bl, tr
        cdef np.ndarray[np.uint8_t,ndim=2] markers

        bl = Vector2Dd(bottomLeft[0], bottomLeft[1])
        tr = Vector2Dd(topRight[0], topRight[1])

        if plane is None:
            raise ImagePlaneException("No plane was specified")

        numX = plane.shape[1]
        numY = plane.shape[0]

        yis, xis = np.where(plane > threshold)
        markers = np.zeros((numY, numX), dtype=np.uint8)

        for i in range(len(yis)):
            xi = xis[i]
            yi = yis[i]
            markers[yi,xi] = 1
        
        segments = [ ]

        for y in range(numY):
            for x in range(numX):
                if markers[y,x]:
                    segment = ImagePlane._findSegmentAt(x, y, markers)
                    if len(segment) > 0:
                        segments.append(segment)

        coordSegments = [ ]
        pixw = (tr.getX()-bl.getX())/numX
        pixh = (tr.getY()-bl.getY())/numY
        for segment in segments:
            cs = [ ]
            for pixelpos in segment:
                x = pixelpos[0]
                y = pixelpos[1]
                X = bl.getX() + pixw*(0.5+x)
                Y = bl.getY() + pixh*(0.5+y)
                cs.append(np.array([X, Y], dtype = np.double))

            coordSegments.append(cs)

        return coordSegments

    def segment(self, np.ndarray[double, ndim=2] plane, threshold = 0.0):
        """segment(plane, threshold = 0.0)

        For the image plane ``plane`` that was rendered using :func:`renderImages`,
        this function looks at all the pixels that have a value larger than ``threshold``.
        These pixels are divided into regions that are coherent, and a list of these
        regions is returned. Each region is itself a list of 2D coordinates describing
        the centers of the pixels.
        """
        numX = self.m_pImgPlane.getNumXPixels()
        numY = self.m_pImgPlane.getNumYPixels()

        ri = self.getRenderInfo()
        bl = ri["bottomleft"]
        tr = ri["topright"]

        if plane is None:
            raise ImagePlaneException("No plane was specified")

        if plane.shape[0] != numY or plane.shape[1] != numX:
            raise ImagePlaneException("Specified plane must have {} rows and {} columns".format(numY, numX))

        return ImagePlane.static_segment(plane, bl, tr, threshold)

    @staticmethod
    def static_traceBetaApproximately(beta, np.ndarray[double, ndim=3] betaMap, bottomLeft, topRight):
        cdef vector[Vector2Dd] betas
        cdef vector[Vector2Dd] thetaPoints
        cdef string errorString
        cdef int numX, numY, x, y

        if betaMap is None:
            raise ImagePlaneException("No beta map was specified")

        if betaMap.shape[2] != 2:
            raise ImagePlaneException("Third dimension of beta map must be 2 (for x and y components of beta vector)")

        numX, numY = betaMap.shape[1], betaMap.shape[0]
        betas.resize(numX*numY)

        for y in range(numY):
            for x in range(numX):
                betas[x+y*numX] = Vector2Dd(betaMap[y,x,0], betaMap[y,x,1])

        if not imageplane.ImagePlane.staticTraceBetaApproximately(Vector2Dd(beta[0], beta[1]), thetaPoints, betas,
                                                           Vector2Dd(bottomLeft[0], bottomLeft[1]),
                                                           Vector2Dd(topRight[0], topRight[1]), numX, numY,
                                                           errorString):
            raise ImagePlaneException(S(errorString))

        points = []
        for i in range(int(thetaPoints.size())):
            points.append([thetaPoints[i].getX(), thetaPoints[i].getY()])

        return points

cdef public api void cy_call_void_string_function(object self, const char *method, char *s):
    try:
        f = getattr(self, S(method))
        f(S(s))
    except Exception as e:
        print("Warning: exception in " + S(method) + ": " + str(e))

cdef public api void cy_call_void_int_function(object self, const char *method, int x):
    try:
        f = getattr(self, S(method))
        f(x)
    except Exception as e:
        print("Warning: exception in " + S(method) + ": " + str(e))

class SourceImageException(Exception):
    """This exception is raised if something goes wrong in one of the :class:`SourceImage` 
    derived classes."""
    pass

cdef class SourceImage:
    """This is a base class for different shapes of sources.

    It should not be instantiated directly, but only through one of the subclasses.
    """
    cdef sourceimage.SourceImage *m_pSrc

    @staticmethod
    cdef sourceimage.SourceImage * _getSourceImage(SourceImage s):
        return s.m_pSrc

    def __cinit__(self):
        self.m_pSrc = NULL

    def __dealloc__(self):
        del self.m_pSrc

    cdef _check(self):
        if self.m_pSrc == NULL:
            raise SourceImageException("No internal source image has been set")

    def getAngularPosition(self):
        """getAngularPosition()

        Returns the 2D position of this source in the source plane.
        """
        cdef Vector2Dd pos
        
        self._check()
        pos = self.m_pSrc.getAngularPosition()
        return np.array([pos.getX(), pos.getY()])

    def addToAngularPosition(self, p):
        """addToAngularPosition(p)

        Adds a 2D vector to the position of this source in the source plane.
        """
        self._check()
        self.m_pSrc.addToAngularPosition(Vector2Dd(p[0], p[1]))

    def setAngularPosition(self, p):
        """setAngularPosition(p)

        Sets the 2D position of the source in the source plane to some value.
        """
        self._check()
        self.m_pSrc.setAngularPosition(Vector2Dd(p[0], p[1]))

    def addToAngle(self, ang):
        """addToAngle(ang)

        Adds the value ``ang`` to the rotation angle of this source. The value
        must be specified in degrees.
        """
        self._check()
        self.m_pSrc.addToAngle(float(ang))

    def setAngle(self, ang):
        """setAngle(ang)

        Sets the rotation angle of this source to ``ang``, which must be specified
        in degrees.
        """
        self._check()
        self.m_pSrc.setAngle(float(ang))

    def getAngle(self):
        """getAngle()

        Returns the rotation angle for the source shape, specified in degrees.
        """
        self._check()
        return self.m_pSrc.getAngle()

    def getMaximumRadius(self):
        """getMaximumRadius()

        Returns the radius outside of which the source does not produce any light."
        """
        self._check()
        return self.m_pSrc.getMaxRadius()

    def getIntensity(self, p):
        """getIntensity(p)

        Returns the intensity of the source at the specified angular position (in the source plane).
        """
        p = np.array(p)
        return self._reshapeAndCall1D(lambda x: self._getIntensity1(x), p, 2, 1)

    cdef _getIntensity1(self, np.ndarray[double, ndim=1] positions):
        
        cdef np.ndarray[double,ndim=1] intens
        cdef int l,i,j

        if positions.shape[0] % 2 == 0:
            l = positions.shape[0]
            intens = np.zeros([(l//2)], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                intens[j] = self.m_pSrc.getIntensity(Vector2Dd(positions[i], positions[i+1]))
                j += 1
        else:
            raise SourceImageException("Bad 1D array dimensions, must be a multiple of two")

        return intens

    # For now, just copied the code from 'lenses'. Can we import some common code here?
    cdef _reshapeAndCall1D(self, functionName, thetas, int coreNumIn, int coreNumOut):
        cdef int totalElements, l, i

        self._check()

        l = len(thetas.shape)
        if l == 1:
            return functionName(thetas)

        if l > 1 and thetas.shape[l-1] == coreNumIn:
            totalElements = 1
            for i in range(l):
                totalElements *= thetas.shape[i]

            outShape = thetas.shape[:-1] 
            if coreNumOut > 1:
                outShape += (coreNumOut,)
            return np.reshape(functionName(np.reshape(thetas,[totalElements])), outShape)

        raise SourceImageException("Bad array dimensions")

cdef class CircularSource(SourceImage):
    def __init__(self, position, angularRadius, brightnessScale = 1.0, fade = False):
        """__init__(position, angularRadius, brightnessScale = 1.0, fade = False)

        Creates a circular source shape with center at 2D location ``position``,
        and with radius ``angularRadius``. The central brightness of the source is
        specified by ``brightnessScale``, and depending on ``fade`` the brightness
        will either stay constant within the circular region or will fade to the border.
        """
        cdef Vector2Dd pos = Vector2Dd(position[0], position[1])
        cdef double r = float(angularRadius)
        cdef double s = float(brightnessScale)
        cdef cbool f = fade
        cdef sourceimage.CircularSource *pSrc = new sourceimage.CircularSource(pos, r, s)

        pSrc.setFade(f)
        self.m_pSrc = pSrc

    cdef sourceimage.CircularSource * _src(self):
        self._check()
        if self.m_pSrc.getSourceType() != sourceimage.Circle:
            raise SourceImageException("Internal error: source image is not of type CircularSource")
        return <sourceimage.CircularSource *>self.m_pSrc

    def getAngularRadius(self):
        """getAngularRadius()

        Returns the radius of the source in the source plane, as specified in the
        constructor.
        """
        return self._src().getAngularRadius()

    def setAngularRadius(self, a):
        """setAngularRadius(a)

        Sets the radius of the source to ``a``.
        """
        cdef double r = float(a)
        self._src().setAngularRadius(r)

    def setFade(self, cbool f):
        """setFade(f)

        Adjusts the fade parameter (see constructor) to ``f``.
        """
        self._src().setFade(f)

    def getFade(self):
        """getFade()

        Returns a boolean indicating if the source brightness fades to zero towards the
        boundary of the circular region, or if it stays constant.
        """
        return self._src().getFade()

cdef class EllipticalSource(SourceImage):
    def __init__(self, position, halfAxis, eccentricity, angle = 0.0, brightnessScale = 1.0, fade = False):
        """__init__(position, halfAxis, eccentricity, angle = 0.0, brightnessScale = 1.0, fade = False)

        Creates an elliptical source shape at location `position`` in the source plane,
        with half long axis ``halfAxis`` and eccentricity described by ``eccentricity``. The
        shape is rotated counter clockwise over an angle ``angle`` (in degrees). The central 
        brightness of the source is specified by ``brightnessScale``, and depending on ``fade`` 
        the brightness will either stay constant within the circular region or will fade to 
        the border.
        """
        cdef Vector2Dd pos = Vector2Dd(position[0], position[1])
        cdef double axis = float(halfAxis)
        cdef double ecc = float(eccentricity)
        cdef double s = float(brightnessScale)
        cdef double ang = float(angle)
        cdef cbool f = fade
        cdef sourceimage.EllipticalSource *pSrc = new sourceimage.EllipticalSource(pos, axis, ecc, ang, s)

        pSrc.setFade(f)
        self.m_pSrc = pSrc

    cdef sourceimage.EllipticalSource * _src(self):
        self._check()
        if self.m_pSrc.getSourceType() != sourceimage.Ellipse:
            raise SourceImageException("Internal error: source image is not of type EllipticalSource")
        return <sourceimage.EllipticalSource *>self.m_pSrc

    def setFade(self, cbool f):
        """setFade(f)

        Adjusts the fade parameter (see constructor) to ``f``.
        """
        self._src().setFade(f)

    def getFade(self):
        """getFade()

        Returns a boolean indicating if the source brightness fades to zero towards the
        boundary of the circular region, or if it stays constant.
        """
        return self._src().getFade()

cdef class PolygonSource(SourceImage):
    def __init__(self, position, polygonPoints, cbool calcHull, angle = 0, brightnessScale = 1.0):
        """__init__(position, polygonPoints, cbool calcHull, angle = 0, brightnessScale = 1.0)

        Uses a polygon as a source shape, and places the first point at location ``position``.
        The polygon can be specified by specifying a list of points in ``polygonPoints``.
        In this case, the polygon is closed automatically, so the last point should not be
        the same as the first. For this usage the ``calcHull`` parameter should be ``False``.

        Alternatively, the ``polygonPoints`` can also list a number of points of which the
        convex hull needs to be calculated and used as the source shape. In this case, the
        ``calcHull`` parameter should be set to ``True``.

        The polygon shape is rotated counter clockwise over angle ``angle`` (in degrees),
        and a brightness scale of ``brightnessScale`` is used when rendering the source.
        """
        cdef Vector2Dd pos = Vector2Dd(position[0], position[1])
        cdef double ang = float(angle)
        cdef double s = float(brightnessScale)
        cdef vector[Vector2Dd] polyPoints
        cdef double x, y
        cdef polygon2d.Polygon2Dd polygon

        if len(polygonPoints) < 3:
            raise SourceImageException("At least three points are required for a polygon source")
        for p in polygonPoints:
            if len(p) != 2:
                raise SourceImageException("Each point in the polygon source should have two coordinates")
            
            x = p[0]
            y = p[1]
            polyPoints.push_back(Vector2Dd(x,y))

        polygon.init(polyPoints, calcHull)

        cdef sourceimage.PolygonSource *pSrc = new sourceimage.PolygonSource(pos, polygon, ang, s)
        self.m_pSrc = pSrc

cdef class DiscreteSource(SourceImage):
    def __init__(self, np.ndarray[double,ndim=2] data, angularWidth, angularHeight,
                 position, angle = 0.0, brightnessScale = 1.0):
        """__init__(data, angularWidth, angularHeight, position, angle = 0.0, brightnessScale = 1.0)

        Uses the 2D NumPy array ``data`` as a pixellated source. This image has width
        and height described by ``angularWidth`` and ``angularHeight``, and the
        center of the source image is placed at location ``position``. It is rotated
        counter clockwise over the angle in ``angle`` (in degrees), and the pixel
        values will be scaled by the specified brightness scale ``brightnessScale``.
        """

        cdef Vector2Dd pos = Vector2Dd(position[0], position[1])
        cdef double w = float(angularWidth)
        cdef double h = float(angularHeight)
        cdef double a = float(angle)
        cdef double s = float(brightnessScale)
        cdef vector[double] dataVector
        cdef int numX, numY, x, y, p
        cdef sourceimage.DiscreteSource *pSrc = new sourceimage.DiscreteSource(pos, a, s)
    
        try:
            if data is None:
                raise SourceImageException("Grid data for source may not be None")

            numY = data.shape[0]
            numX = data.shape[1]
            dataVector.resize(numX*numY)
            p = 0
            for y in range(numY):
                for x in range(numX):
                    dataVector[p] = data[y,x]
                    p += 1

            if not pSrc.setData(dataVector, numX, numY, w, h):
                raise LensPlaneException(S(pSrc.getErrorString()))

            self.m_pSrc = pSrc
        except:
            del pSrc
            raise

cdef class PointSource(SourceImage):
    def __init__(self, position, brightnessScale = 1.0):
        """__init__(position, brightnessScale = 1.0)

        Places a point source at location ``position``, with the specified brightness scale.
        """
        cdef Vector2Dd pos = Vector2Dd(position[0], position[1])
        cdef double s = float(brightnessScale)

        self.m_pSrc = new sourceimage.PointSource(pos, s)

def hoursMinutesSecondsToDegrees(s):
    """hoursMinutesSecondsToDegrees(s)
    
    Converts a string or array of three parts, specifying hours, minutes and seconds,
    into a single floating point number that corresponds to a number of degrees.

    As an example, passing ``"01:23:45"``, ``"01 23 45"``, ``["01","23","45"]`` and
    ``[1,23,45]`` will all produce the same output of 20.9375.
    """

    if type(s) != list:
        if ":" in s:
            s = s.split(":")
        else: # Perhaps a space is used?
            s = s.split(" ")
            
    if len(s) != 3:
        raise Exception("Expected three parts, but detected {}".format(len(s)))
        
    h, m, s = float(s[0]), float(s[1]), float(s[2])
    
    if h < 0 or h >= 24:
        raise Exception("'Hours' must lie between 0 and 24, but is {}".format(h))
    if m < 0 or m >= 60:
        raise Exception("'Minutes must lie between 0 and 60, but is {}".format(m))
    if s < 0 or s >= 60:
        raise Exception("'Seconds must lie between 0 and 60, but is {}".format(s))
    
    return (((s / 60.0) + m)/60.0 + h)/24.0 * 360.0

def degreesMinutesSecondsToDegrees(s):
    """degreesMinutesSecondsToDegrees(s)
    
    Converts a string or array of three parts, specifying degrees, minutes and seconds,
    into a single floating point number that corresponds to a number of degrees.

    As an example, passing ``"-1:23:45"``, ``"-1 23 45"``, ``["-1","23","45"]`` and
    ``[-1,23,45]`` will all produce the same output of -1.39583.
    """
    if type(s) != list:
        if ":" in s:
            s = s.split(":")
        else: # Perhaps a space is used
            s = s.split(" ")
            
    if len(s) != 3:
        raise Exception("Expected three parts, but detected {}".format(len(s)))
        
    d, m, s = float(s[0]), float(s[1]), float(s[2])
    
    sign = 1.0
    if d < 0:
        sign = -1.0
        d = -d

    if m < 0 or m >= 60:
        raise Exception("'Minutes must lie between 0 and 60, but is {}".format(m))
    if s < 0 or s >= 60:
        raise Exception("'Seconds must lie between 0 and 60, but is {}".format(s))
    
    return (((s / 60.0) + m)/60.0 + d) * sign

# For internal use
def _processImagesDataDict(imgs):
    if not imgs:
        raise Exception("Empty images set")
    
    imgDat = ImagesData(len(imgs))
    groupInfo = { }
    
    imgIdx = 0
    for imgId in imgs:
        imgPoints = imgs[imgId]
        if len(imgPoints) == 0:
            raise Exception("Empty image")
            
        for ptInfo in imgPoints:
            ptIdx = imgDat.addPoint(imgIdx, [ ptInfo["x"], ptInfo["y"] ])
            if "group" in ptInfo:
                grpId = ptInfo["group"]
                if not grpId in groupInfo:
                    groupInfo[grpId] = []
                groupInfo[grpId].append((imgIdx, ptIdx))
                
            if "timedelay" in ptInfo:
                imgDat.addTimeDelayInfo(imgIdx, ptIdx, ptInfo["timedelay"])
        
        imgIdx += 1
        
    # Process the stored groups
    for grpKey in groupInfo:
        grpIdx = imgDat.addGroup()
        for imgIdx, ptIdx in groupInfo[grpKey]:
            imgDat.addGroupPoint(grpIdx, imgIdx, ptIdx)
    
    return imgDat

def readInputImagesFile(inputData, lineAnalyzer, isPointImagesFile, centerOn = [0, 0]):
    """readInputImagesFile(inputData, lineAnalyzer, isPointImagesFile, centerOn = [0, 0])
    
    This function can process a text file (or previously read text data) into
    one or more :class:`ImagesData` instances.

    If the function you pass in ``lineAnalyzer`` returns a source identifier
    and image identifier, these will be used to group images per source. Otherwise,
    the behaviour is different depending on the value of ``isPointImagesFile``. If
    ``True``, then each single line is interpreted as a single point image, and a
    blank line groups the point images per source. If ``isPointImagesFile`` is ``False``,
    then a single blank line separates the points that belong to different images,
    and a double blank link separates the images from different sources.

    Arguments:
     - `inputData`: this can be an open file object, a filename or just
       the text data that's already read from a file. Lines that start
       with a ``#`` sign are considered to be comments, and are ignored
       completely.

     - `lineAnalyzer`: here you should pass a function that interprets a single
       (non-empty) line in the input data. It should return a dictionary with
       these entries:
       
        - ``x``: the x-coordinate of the point, converted to radians (the
          :ref:`pre-defined constants <constants>` can be useful here).
          The helper functions :func:`degreesMinutesSecondsToDegrees` and
          :func:`hoursMinutesSecondsToDegrees` can also be of assistance.
        - ``y``: the y-coordinate of the point, converted to radians (the
          :ref:`pre-defined constants <constants>` can be useful here)
          The helper function :func:`degreesMinutesSecondsToDegrees` can also be of assistance.
        - (optionally) ``srcnr``: the source identifier of this point.
        - (optionally) ``imgnr``: the image identifier of this point.
        - (optionally) ``z``: the redshift for the source. Note that if this
          is specified for more than one image in a source, they should all
          have exactly the same redshift.
        - (optionally): ``group``: if certain points in different images
          belong together, this can be specified by having the same group
          identifier.
        - (optionally): ``timedelay``: if timedelay information is known for
          a point, it can be passed along this way.

     - `isPointImagesFile`: as explained above, this flag indicates how input
       is treated when no source and image identifiers are specified in the
       `lineAnalyzer` function.

     - `centerOn`: if present, all coordinates will be recentered on this
       point, using the :func:`ImagesData.centerOnPosition` function.

    The return value of this function is a list of dictionaries with entries:

     - ``imgdata``: the :class:`ImagesData` object for a specific source
     - (optionally) ``z``: if it was present in the input, the redshift for
       the source is stored.
     - (optionally) ``srcnr``: if a specific source identifier was specified
       in the file, and passed along in the `lineAnalyzer` function, it is
       stored here.

    Examples: suppose that each line contains entries like the one
    below (from `A free-form lensing grid solution for A1689 with new multiple images <http://adsabs.harvard.edu/abs/2015MNRAS.446..683D>`_) ::

        # i  ID   B05  REF  RAJ2000(h:m:s)  DECJ2000(d:m:s)  z     Delta_beta
        1    1.1  1.1  B05  13:11:26.257    -1:19:58.753     3.04  1.03
        2    1.2  1.2  B05  13:11:26.088    -1:20:02.261     3.04  0.73
        3    1.3  1.3  B05  13:11:29.584    -1:21:09.475     3.04  2.50
        ...

    Then, you could use the following function as the `lineAnalyzer` parameter:

    .. code-block:: python

        def la(l):
            p = l.split()
            if len(p) != 8:
                raise Exception("Expecting 8 input fields in line '{}'".format(l))
            
            info = { }
            info["srcnr"], info["imgnr"] = map(int, p[1].split("."))
            info["z"] = p[6]
            info["x"] = hoursMinutesSecondsToDegrees(p[4])*ANGLE_DEGREE
            info["y"] = degreesMinutesSecondsToDegrees(p[5])*ANGLE_DEGREE
            
            return info

    If we'd like to center the coordinates on that first point, we could calculate

    .. code-block:: python

        centerRa = hoursMinutesSecondsToDegrees("13:11:26.257")*ANGLE_DEGREE
        centerDec = degreesMinutesSecondsToDegrees("-1:19:58.753")*ANGLE_DEGREE

    and pass ``[centerRa, centerDec]`` as the `centerOn` parameter.

    As another example, in `Dark matter dynamics in Abell 3827: new data consistent with standard Cold Dark Matter <http://adsabs.harvard.edu/abs/2017arXiv170804245M>`_
    you'll find coordinates like the ones below::

        # Name  RA         Dec
        Ao1     330.47479  -59.94358
        Ao2     330.46649  −59.94665
        ...


    The 'A' is a label for the source, and as all points refer to the same image, is present
    for each point. The second character, 'o' in this case, refers to a feature in an image.
    The third part, '1' and '2' above, indicates the image the point is a part of. We can
    treat this input in two ways: the image number goes from 1 to 7, so we can process the
    input in the following way to create seven images, where each point is part of a specific
    group:

    .. code-block:: python

        def la(l):
            p = l.split()
            if len(p) != 3:
                raise Exception("Expecting 3 input fields in '{}'".format(l))
            
            info = { }
            info["srcnr"] = 'A'
            info["imgnr"] = p[0][2]
            info["group"] = p[0][1]
            info["x"] = float(p[1])*ANGLE_DEGREE
            info["y"] = float(p[2])*ANGLE_DEGREE
            
            return info

    Alternatively, we can treat each feature, each set of corresponding points (marked with 'o'
    for example), as a different point source. This is wat the following `lineAnalyzer`
    function would do:

    .. code-block:: python

        def la(l):
            p = l.split()
            if len(p) != 3:
                raise Exception("Expecting 3 input fields in '{}'".format(l))
            
            info = { }
            info["srcnr"] = p[0][1]
            info["imgnr"] = p[0][2]
            info["x"] = float(p[1])*ANGLE_DEGREE
            info["y"] = float(p[2])*ANGLE_DEGREE
            
            return info
    """

    try:
        # Check if it's a file that's already open
        lines = inputData.readlines()
    except:
        try:
            with open(inputData, "rt") as f:
                lines = f.readlines()
        except:
            # Check if it's a string that can be split in lines
            lines = inputData.splitlines()

    internalIdPrefix = "justsomelongandunlikelyprefix_"
    internalSourceId = 0
    internalImageId = 0
    sourceCount = 0 # To be able to keep the same ordering of the sources as in the file
    
    def idStr(x):
        return "{}_{:05d}".format(internalIdPrefix, x)
    
    emptyLineCount = 0
    sourceData = { }
    
    for l in lines:
        l = l.strip()
        if l.startswith("#"): # a commented line, ignore completely
            continue
        
        if not l:
            emptyLineCount += 1
        else:
            emptyLineCount = 0
          
        if emptyLineCount == 0: # New info
            if isPointImagesFile:
                internalImageId += 1
            
            info = lineAnalyzer(l)
            imagesId = idStr(internalImageId) if not "imgnr" in info else info["imgnr"]
            sourceId = idStr(internalSourceId) if not "srcnr" in info else info["srcnr"]
            
            if not sourceId in sourceData:
                imagesData = { "points": { } }
                sourceData[sourceId] = imagesData
                sourceData[sourceId]["sourcecount"] = sourceCount
                sourceCount += 1
            else:
                imagesData = sourceData[sourceId]
                
            if not imagesId in imagesData["points"]:
                imagePoints = [ ]
                imagesData["points"][imagesId] = imagePoints
            else:
                imagePoints = imagesData["points"][imagesId]
            
            if "z" in info:
                if "z" in imagesData:
                    if info["z"] != imagesData["z"]:
                        raise Exception("Error in line '{}': redshift '{}' is not the same as previous recorded redshift '{}'".format(info["z"],imagesData["z"]))
                
                imagesData["z"] = info["z"]
            
            ptInfo = {
                "x": float(info["x"]),
                "y": float(info["y"])
            }
            if "group" in info:
                ptInfo["group"] = info["group"]
            if "timedelay" in info:
                ptInfo["timedelay"] = info["timedelay"]
            
            imagePoints.append(ptInfo)
    
        elif emptyLineCount == 1 and not isPointImagesFile:
            internalImageId += 1
        elif emptyLineCount == 2 or (emptyLineCount == 1 and isPointImagesFile):
            internalSourceId += 1
            internalImageId = 0
    
    srcListTmp = sorted([ (sourceData[s]["sourcecount"], s) for s in sourceData])
        
    srcList = [ ]
    for tmpCount, srcId in srcListTmp:
        info = { }
        if "z" in sourceData[srcId]:
            info["z"] = sourceData[srcId]["z"]
        if not str(srcId).startswith(internalIdPrefix):
            info["srcnr"] = srcId
        
        info["imgdata"] = _processImagesDataDict(sourceData[srcId]["points"])
        info["imgdata"].centerOnPosition(centerOn[0], centerOn[1])
        srcList.append(info)
    
    # TODO: if point images file, check that every images data set is indeed a point image set
    #       if not a point image set, _optionally_ check that they're not point images
    
    #return sourceData
    return srcList
