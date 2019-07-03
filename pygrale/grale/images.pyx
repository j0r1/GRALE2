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
 - :class:`PointSource`: describes a point source
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

    @staticmethod
    cdef void _swapImagesData(ImagesData i1, ImagesData i2):
        cdef imagesdata.ImagesData *pTmp = i1.m_pImgData
        i1.m_pImgData = i2.m_pImgData
        i2.m_pImgData = pTmp
    
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

    def addPoint(self, int imageNum, point, intensity = None, shear = None, shearWeight = 1):
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
        cdef vector2d.Vector2Dd sh

        if intensity is None:
            if shear is None:
                pointNum = self.m_pImgData.addPoint(imageNum, p)
            else:
                sh = vector2d.Vector2Dd(shear[0], shear[1])
                pointNum = self.m_pImgData.addPoint(imageNum, p, sh, shearWeight)
        else:
            if shear is None:
                pointNum = self.m_pImgData.addPoint(imageNum, p, intensity)
            else:
                sh = vector2d.Vector2Dd(shear[0], shear[1])
                pointNum = self.m_pImgData.addPoint(imageNum, p, intensity, sh, shearWeight)

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

    def setImagePointPosition(self, image, point, position):
        """setImagePointPosition(image, point, position)

        Changes the stored position for the point with ID `point` and image with ID
        `image`, to the coordinates in `pos`.
        """
        cdef vector2d.Vector2Dd p = vector2d.Vector2Dd(position[0], position[1])

        self._checkImagePointNumber(image, point)
        self.m_pImgData.setImagePointPosition(image, point, p)

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
                obj["shearweight"] = self.getShearWeight(image, i)

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

    def getShearWeight(self, image, point):
        """getShearComponent2(image, point)

        Returns the second shear component stored for the point with ID
        ``point`` inside the image with ID ``image``.
        """
        self._checkImagePointNumber(image, point)
        self._checkShear()
        return self.m_pImgData.getShearWeight(image, point)

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

    def getConvexHull(self, image, asIndices = False):
        """getConvexHull(image, asIndices = False)

        Returns the convex hull of the points for image ID `image`. By default, a list of coordinates is returned,
        but if `asIndices` is ``True``, then the point indices are returned instead.
        If successful, the first and last point in the list will be the same.
        """
        from scipy.spatial import ConvexHull

        num = self.getNumberOfImagePoints(image)
        points = [ self.getImagePointPosition(image, pt) for pt in range(num) ]
        if len(points) < 3:
            raise ImagesDataException("Need at least three points for a convex hull")

        hull = ConvexHull(points)
        if asIndices:
            h = [ pt for pt in hull.vertices ]
        else:
            h = [ self.getImagePointPosition(image, pt) for pt in hull.vertices ]
        return h + h[:1]

    def getBorder(self, image, asIndices = False):
        """getBorder(image, asIndices = False)

        For the image with ID `image`, return the border, based on the triangulation
        that's stored for this image. By default, a list of coordinates is returned,
        but if `asIndices` is ``True``, then the point indices are returned instead.
        If successful, the first and last point in the list will be the same.
        """
        tList = self.getTriangles(image)
        edgeCount = { }

        def getEdges(i1, i2, i3):
            return [ tuple(sorted([a,b])) for a,b in [(i1,i2), (i2,i3), (i1,i3)]]

        for i1,i2,i3 in tList:
            edges = getEdges(i1, i2, i3)
            for e in edges:
                if not e in edgeCount:
                    edgeCount[e] = 1
                else:
                    edgeCount[e] += 1
                    if edgeCount[e] > 2:
                        raise ImagesDataException("Unexpected amount of shared edges")

        border = sorted([ e for e in edgeCount if edgeCount[e] == 1 ])
        if not border:
            raise ImagesDataException("No border found")

        #pprint.pprint(border)
        
        # TODO: check that each vertex is used two times

        end1, end2 = border[0]
        borderPoints = [ end1, end2 ]
        border = border[1:]
        if end1 == end2:
            raise ImagesDataException("Error: edge between identical points")

        while True:
            for idx in range(len(border)):
                i1, i2 = border[idx]
                if i1 == end1:
                    end1 = i2
                    borderPoints = [i2] + borderPoints
                    del border[idx]
                    break
                elif i1 == end2:
                    end2 = i2
                    borderPoints = borderPoints + [i2]
                    del border[idx]
                    break
                elif i2 == end1:
                    end1 = i1
                    borderPoints = [i1] + borderPoints
                    del border[idx]
                    break
                elif i2 == end2:
                    end2 = i1
                    borderPoints = borderPoints + [i1]
                    del border[idx]
                    break
            else:
                raise ImagesDataException("No next edge found")

            if end1 == end2: # We've closed the border
                break

        if border:
            raise ImagesDataException("Closed border, but still have edges left")

        if borderPoints[0] != borderPoints[-1]:
            raise ImagesDataException("Internal error: border not closed as expected")

        #borderPoints.pop()
        #pprint.pprint(borderPoints)

        if asIndices:
            borderPoints = [ ptIdx for ptIdx in borderPoints ]
        else:
            borderPoints = [ self.getImagePointPosition(image, ptIdx) for ptIdx in borderPoints ]
        return borderPoints

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(b)

        This function attempts to interpret the bytes ``b``
        an images data set, and returns the new instance if successful.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(b), NULL, 0)
        cdef imagesdata.ImagesData2 *pImgDat

        img = ImagesData(_imagesDataRndId) # Make sure an empty instance is allocated

        # Try to load the data
        pImgDat2 = <imagesdata.ImagesData2 *>img.m_pImgData
        if not pImgDat2.read2(deref(m)):
            del m
            raise ImagesDataException(S(img.m_pImgData.getErrorString()))
        
        del m
        return img

    def toBytes(self):
        """toBytes()

        Returns a binary representation of the current images data set."
        """
        cdef serut.VectorSerializer vSer
        cdef imagesdata.ImagesData2 *pImgDat2 = <imagesdata.ImagesData2 *>self.m_pImgData

        if not pImgDat2.write2(vSer):
            raise ImagesDataException(S(self.m_pImgData.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

    def __getstate__(self):
        return self.toBytes()

    def __setstate__(self, state):
        img = ImagesData.fromBytes(state)
        ImagesData._swapImagesData(self, img)

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

    cdef void _swapLensPlanes(LensPlane l1, LensPlane l2):
        cdef lensplane.LensPlane *pTmp = l1.m_pLensPlane
        cdef object tmpFeedback = l1.feedback

        l1.m_pLensPlane = l2.m_pLensPlane
        l2.m_pLensPlane = pTmp

        l1.feedback = l2.feedback
        l2.feedback = tmpFeedback

    def __cinit__(self):
        self.m_pLensPlane = NULL

    def __dealloc__(self):
        del self.m_pLensPlane

    def __init__(self, lens, bottomLeft, topRight, int numX, int numY, renderer = "default", feedbackObject = "default"):
        """__init__(lens, bottomLeft, topRight, numX, numY, renderer = "default", feedbackObject = "default")

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

        if renderer and renderer.renderType != "LENSPLANE":
            raise LensPlaneException("Specified renderer does not appear to be for the LENSPLANE type")

        self.feedback = feedbackObject

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
            self.feedback = None

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
        """Returns a dictionary containing entries ``alpha_x`` and
        ``alpha_y``, containing the deflections stored in the lensplane.
        Using the names returned by :func:`getRenderInfo`, each of these
        entries is a NumPy grid of shape ``(ypoints, xpoints)``, containing
        the X and Y components of the deflections. The [0,0] entry of a
        grid corresponds to the value at ``bottomleft``, the 
        ``[ypoints-1, xpoints-1]`` entry corresponds to ``topright``.
        """

        cdef np.ndarray[double,ndim=2] alphax, alphay
        cdef int xpoints, ypoints, x, y;
        cdef Vector2Dd alpha;

        self._check()

        xpoints = self.m_pLensPlane.getNumXPoints()
        ypoints = self.m_pLensPlane.getNumYPoints()
        
        alphax = np.empty([ypoints, xpoints], dtype = np.double)
        alphay = np.empty([ypoints, xpoints], dtype = np.double)
        for y in range(ypoints):
            for x in range(xpoints):
                alpha = self.m_pLensPlane.getAlpha(x, y)
                alphax[y,x] = alpha.getX()
                alphay[y,x] = alpha.getY()

        return { "alpha_x": alphax, "alpha_y": alphay }

    def getAlphaVectorDerivatives(self):
        """Returns a dictionary similar to :func:`getAlphas` but now with
        names ``alpha_xx``, ``alpha_yy`` and ``alpha_xy``, describing 
        respectively the derivative of the X-component of the deflection 
        angle in the X-direction, the derivative of the Y-component of 
        the deflection angle in the Y direction and the derivative of the
        X-component of the deflection angle in the Y-direction.
        """
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

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(b)

        This function attempts to interpret the bytes ``b``
        a lens plane instance, and returns the new instance if successful.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(b), NULL, 0)
        cdef string errorString
        cdef lensplane.LensPlane *pLensPlane = NULL

        if not lensplane.LensPlane.read(deref(m), cython.address(pLensPlane), errorString):
            del m
            raise LensPlaneException(S(errorString))

        del m

        lensPlane = LensPlane(None, None, None, 0, 0)
        lensPlane.m_pLensPlane = pLensPlane

        return lensPlane

    def toBytes(self):
        """toBytes()

        Returns a binary representation of the current images data set."
        """
        cdef serut.VectorSerializer vSer

        self._check()

        if not self.m_pLensPlane.write(vSer):
            raise ImagesDataException(S(self.m_pLensPlane.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

    def __getstate__(self):
        return self.toBytes()

    def __setstate__(self, state):
        lp = LensPlane.fromBytes(state)
        LensPlane._swapLensPlanes(self, lp)

class ImagePlaneException(Exception):
    """This exception is raised if something goes wrong in the :class:`ImagePlane` class."""
    pass

cdef class ImagePlane:
    cdef imageplane.ImagePlane *m_pImgPlane
    # Store the original lens plane (or better, a copy) and the distances to
    # be able to pickle
    cdef object m_lensPlane 
    cdef double m_Ds, m_Dds

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
        self._internalConstructor(lensplane, Ds, Dds)

    def _internalConstructor(self, LensPlane lensplane, double Ds, double Dds):
        if not self.m_pImgPlane.init(<imageplane.LensPlane*>lensplane.m_pLensPlane, Ds, Dds):
            raise ImagePlaneException("Unable to initialize image plane: " + S(self.m_pImgPlane.getErrorString()))

        # Store the constructor values
        import copy
        self.m_lensPlane = copy.copy(lensplane)
        self.m_Ds = Ds
        self.m_Dds = Dds

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

    def traceThetaApproximately(self, thetas):
        """Use the already calculated theta/beta mapping (image plane position to
        source plane positions), to estimate the mapping for theta vectors that
        have not been calculated exactly."""

        thetas = np.array(thetas)
        reshapedThetas = thetas.reshape((-1,2))
        reshapedBetas = self._traceThetaApprox(reshapedThetas)
        return reshapedBetas.reshape(thetas.shape)

    cdef _traceThetaApprox(self, np.ndarray[double, ndim=2] reshapedThetas):
        cdef Vector2Dd beta, theta, bottomLeft, topRight
        cdef np.ndarray[double, ndim=2] betas = np.zeros((reshapedThetas.shape[0], reshapedThetas.shape[1]))
        cdef int num, i

        bottomLeft = self.m_pImgPlane.getBottomLeft()
        topRight = self.m_pImgPlane.getTopRight()

        num = betas.shape[0]
        for i in range(num):
            theta = Vector2Dd(reshapedThetas[i,0], reshapedThetas[i,1])
            if theta.getX() < bottomLeft.getX() or theta.getY() < bottomLeft.getY() or theta.getX() > topRight.getX() or theta.getY() > topRight.getY():
                raise ImagePlaneException("Not all points lie withing the image plane boundaries")

            if not self.m_pImgPlane.traceThetaApproximately(theta, cython.address(beta)):
                raise ImagePlaneException(S(self.m_pImgPlane.getErrorString()))

            betas[i,0] = beta.getX()
            betas[i,1] = beta.getY()

        return betas

    def __getstate__(self):
        return [ self.m_lensPlane.toBytes(), self.m_Ds, self.m_Dds ]

    def __setstate__(self, state):
        self._internalConstructor(LensPlane.fromBytes(state[0]), state[1], state[2])

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

# import some things from a pure python module, to avoid recompilation if something
# changes there.
from privimages import getDefaultLineAnalyzer, setDefaultLineAnalyzer, readInputImagesFile, hoursMinutesSecondsToDegrees, degreesMinutesSecondsToDegrees, createGridTriangles, enlargePolygon, createSourceFromImagesData, createPointImagesData

