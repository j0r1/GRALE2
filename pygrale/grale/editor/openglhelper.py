from PyQt5 import QtWidgets, QtOpenGL, QtGui, QtCore
import numpy as np
import sys

def needcontext(f):

    def wrapper(*args, **kwargs):
        
        isOpenGLObject = False
        if len(args) > 0:
            if issubclass(type(args[0]), OpenGLObject):
                isOpenGLObject = True
                self = args[0]

        if isOpenGLObject:
            if not self.ctx.makeCurrent(self.surface):
                raise Exception("Unable to make context current")

            return f(*args, **kwargs)

        else:
            if not QtCore.QCoreApplication.instance():
                raise Exception("No Qt application instance has been created yet")
            t = QtCore.QThread.currentThread()
            if not t:
                raise Exception("Couldn't obtain current thread")
            d = t.eventDispatcher()
            if not d:
                raise Exception("Event dispatcher not found")

            surface = QtGui.QOffscreenSurface()
            fmt = surface.format()
            fmt.setVersion(2, 0)
            surface.setFormat(fmt)

            surface.create()
            if not surface.isValid():
                raise Exception("Couldn't create valid surface")

            ctx = QtGui.QOpenGLContext()
            if not ctx.create():
                raise Exception("Coudln't create OpenGL context")

            if not ctx.makeCurrent(surface):
                raise Exception("Unable to make context current")

            p = QtGui.QOpenGLVersionProfile()
            p.setVersion(2, 0)
            gl = ctx.versionFunctions(p)
            if not gl:
                raise Exception("Couldn't get access to OpenGL function")

            return f(*args, **kwargs, gl=gl)

    return wrapper

class OpenGLObject(QtCore.QObject):
    def __init__(self, parent = None):
        super(OpenGLObject, self).__init__(parent)

        if not QtCore.QCoreApplication.instance():
            raise Exception("No Qt application instance has been created yet")
        t = QtCore.QThread.currentThread()
        if not t:
            raise Exception("Couldn't obtain current thread")
        d = t.eventDispatcher()
        if not d:
            raise Exception("Event dispatcher not found")

        self.surface = QtGui.QOffscreenSurface()
        fmt = self.surface.format()
        fmt.setVersion(2, 0)
        self.surface.setFormat(fmt)

        self.surface.create()
        if not self.surface.isValid():
            raise Exception("Couldn't create valid surface")

        self.ctx = QtGui.QOpenGLContext()
        if not self.ctx.create():
            raise Exception("Coudln't create OpenGL context")

        self._getFunctions()

    @needcontext
    def _getFunctions(self):
        p = QtGui.QOpenGLVersionProfile()
        p.setVersion(2, 0)
        self.gl = self.ctx.versionFunctions(p)
        if not self.gl:
            raise Exception("Couldn't get access to OpenGL function")

class BackProjectObject(OpenGLObject):
    def __init__(self, parent = None):
        super(BackProjectObject, self).__init__(parent)

    @needcontext
    def backProject(self, imgPlane, inputImage, center, sizes, outputDimensions):
        return _backProjectInternal(imgPlane, inputImage, center, sizes, outputDimensions, self.gl)

@needcontext
def trace(imgPlane, srcImage, center, sizes, outputDimensions, subSample, **kwargs):
    gl = kwargs["gl"]
    return _trace(imgPlane, srcImage, center, sizes, outputDimensions, subSample, gl)

def _trace(imgPlane, srcImage, center, sizes, outputDimensions, subSample, gl):

    print("_trace", center, sizes)
    from grale.constants import ANGLE_ARCSEC

    ri = imgPlane.getRenderInfo()
    bl, tr = np.array(ri["bottomleft"]), np.array(ri["topright"])
    mapNumX, mapNumY = ri["xpoints"], ri["ypoints"]

    fb = QtGui.QOpenGLFramebufferObject(*outputDimensions)
    if not fb.bind():
        raise Exception("Couldn't bind framebuffer")


    #imgPlaneRect = _getImagePlaneRect(imgPlane)

    eps = 0
    mapTexData = np.empty([mapNumY, mapNumX, 2], dtype=np.double)
    mapTexData[:,:,0], mapTexData[:,:,1] = np.meshgrid(np.linspace(bl[0]+eps, tr[0]-eps, mapNumX),
                                                       np.linspace(bl[1]+eps, tr[1]-eps, mapNumY))

    mapTexData = (imgPlane.traceThetaApproximately(mapTexData)/ANGLE_ARCSEC).astype(np.float32).reshape((-1,))

    GL_TEXTURE0 = 0x84C0
    GL_TEXTURE_2D = 0x0DE1
    GL_TEXTURE_MAG_FILTER = 0x2800
    GL_TEXTURE_MIN_FILTER = 0x2801
    GL_LINEAR = 0x2601
    GL_TEXTURE_WRAP_S = 0x2802
    GL_TEXTURE_WRAP_T = 0x2803
    GL_CLAMP_TO_EDGE = 0x812F
    GL_RG32F = 0x8230
    GL_RG = 0x8227
    GL_FLOAT = 0x1406

    mapTexture = gl.glGenTextures(1)
    gl.glActiveTexture(GL_TEXTURE0)
    gl.glBindTexture(GL_TEXTURE_2D, mapTexture);
    gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    b = mapTexData.tobytes()
    import array
    a = array.array("f")
    a.frombytes(b)
    gl.glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, mapNumX, mapNumY, 0, GL_RG, GL_FLOAT, a)

    srcTexture = QtGui.QOpenGLTexture(QtGui.QOpenGLTexture.Target2D)
    if not srcTexture.create():
        raise Exception("Unable to create texture")

    srcTexture.bind(1)
    srcTexture.setWrapMode(QtGui.QOpenGLTexture.ClampToEdge)
    srcTexture.setData(srcImage)

    print("1 glGetError =", gl.glGetError())

    prog = QtGui.QOpenGLShaderProgram()
    prog.addShaderFromSourceCode(QtGui.QOpenGLShader.Vertex,"""
		precision highp float;

		attribute vec2 a_position;
		varying vec2 v_tpos;
        uniform vec2 u_mapPixSize;

		void main()
		{
            vec2 xy = (a_position+1.0)/2.0;
            // We need to rescale slightly because the map values 
            // are at the center of the pixels
			v_tpos = (1.0-u_mapPixSize)*xy + 0.5*u_mapPixSize;
			gl_Position = vec4(a_position.xy, 0, 1);
		}
        """)

    prog.addShaderFromSourceCode(QtGui.QOpenGLShader.Fragment,"""
		precision highp float;

		uniform sampler2D u_mapTexture;
		uniform sampler2D u_srcTexture;
		uniform vec2 u_srcBottomLeft;
		uniform vec2 u_srcTopRight;
        uniform vec2 u_dstart;
        uniform vec2 u_substep;
		varying vec2 v_tpos;

		void main()
		{{
            vec4 c = vec4(0,0,0,0);
            float y = v_tpos.y+u_dstart.y;
            
            for (int i = 0 ; i < {subSample} ; i++)
            {{
                float x = v_tpos.x+u_dstart.x;

                y += u_substep.y;
                
                for (int j = 0 ; j < {subSample} ; j++)
                {{
                    x += u_substep.x;

                    vec2 srcCoord = texture2D(u_mapTexture, vec2(x,y)).xy;
                    srcCoord -= u_srcBottomLeft;
                    srcCoord /= (u_srcTopRight-u_srcBottomLeft);

                    if (all(greaterThanEqual(srcCoord, vec2(0.0,0.0))) && all(lessThanEqual(srcCoord, vec2(1.0,1.0))))
                        c += texture2D(u_srcTexture, srcCoord);
                }}
            }}
            gl_FragColor = c/({subSample}*{subSample});
		}}
    """.format(subSample=subSample))
    if not prog.link():
        raise Exception("Unable to link program")
    prog.bind()

    print("1 glGetError =", gl.glGetError())

    coords = np.array([-1, -1, 
                       -1, 1,
                       1, -1,
                       1, 1], dtype=np.float32)

    buf = QtGui.QOpenGLBuffer(QtGui.QOpenGLBuffer.VertexBuffer)
    if not buf.create():
        raise Exception("Couldn't create buffer")

    buf.bind()
    buf.setUsagePattern(QtGui.QOpenGLBuffer.StaticDraw)
    buf.allocate(coords, coords.shape[0]*4)

    posLoc = prog.attributeLocation("a_position")
    prog.setAttributeBuffer(posLoc, gl.GL_FLOAT, 0, 2)
    prog.enableAttributeArray(posLoc)

    subStepWidth, subStepHeight = [ (1.0/subSample)/outputDimensions[i] for i in range(2) ]
    dStart = [ (0.5/subSample-0.5)/outputDimensions[i] for i in range(2) ] 

    prog.setUniformValue("u_mapPixSize", QtCore.QPointF(1.0/mapNumX, 1.0/mapNumY))
    prog.setUniformValue("u_mapTexture", 0)
    prog.setUniformValue("u_srcTexture", 1)
    prog.setUniformValue("u_srcBottomLeft", QtCore.QPointF(center[0]-sizes[0]/2, center[1]-sizes[1]/2))
    prog.setUniformValue("u_srcTopRight", QtCore.QPointF(center[0]+sizes[0]/2, center[1]+sizes[1]/2))
    prog.setUniformValue("u_dstart", QtCore.QPointF(dStart[0], dStart[1]))
    prog.setUniformValue("u_substep", QtCore.QPointF(subStepWidth, subStepHeight))

    gl.glViewport(0, 0, fb.width(), fb.height())
    gl.glClearColor(0, 0, 0, 0)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT)
    gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, 0, coords.shape[0]//2)

    gl.glFinish()
    gl.glDeleteTextures(1, [mapTexture])
    print("glGetError =", gl.glGetError())

    img = fb.toImage()
    if img.isNull():
        raise Exception("Unable to get image from framebuffer")

    gl.glFinish()
    print("2 glGetError =", gl.glGetError())

    return img

def _getCoords(w, h, x0, x1, y0, y1, mapFunction):
    x = np.linspace(x0, x1, w+1)
    y = np.linspace(y0, y1, h+1)
    xy = np.empty([h+1, w+1, 2], dtype=np.double)
    xy[:,:,0], xy[:,:,1] = np.meshgrid(x, y)
    xy = mapFunction(xy)
    xy = xy.astype(np.float32)

    c0 = xy[:-1,:-1,:]
    c1 = xy[1:,:-1,:]
    c2 = xy[:-1,1:,:]
    c3 = xy[1:,1:,:]

    coords = np.empty([h, w, 12], dtype=np.float32)
    coords[:,:,0:2] = c0
    coords[:,:,2:4] = c1
    coords[:,:,4:6] = c2
    coords[:,:,6:8] = c3
    coords[:,:,8:10] = c2
    coords[:,:,10:12] = c1
    coords = coords.reshape((-1,))

    buf = QtGui.QOpenGLBuffer(QtGui.QOpenGLBuffer.VertexBuffer)
    if not buf.create():
        raise Exception("Couldn't create buffer")

    buf.bind()
    buf.setUsagePattern(QtGui.QOpenGLBuffer.StaticDraw)
    buf.allocate(coords, coords.shape[0]*4)
    return coords, buf

@needcontext
def backProject(imgPlane, inputImage, center, sizes, outputDimensions, **kwargs):
    gl = kwargs["gl"]
    return _backProjectInternal(imgPlane, inputImage, center, sizes, outputDimensions, gl) # gl is set by decorator

def _backProjectInternal(imgPlane, inputImage, center, sizes, outputDimensions, gl):

    if inputImage.isNull():
        raise Exception("Input image is not valid")

    lrbt = [ center[0]-sizes[0]/2.0, center[0]+sizes[0]/2.0, center[1]-sizes[1]/2.0, center[1]+sizes[1]/2.0 ]
    outputWH = [ outputDimensions[0], outputDimensions[1] ]

    from grale.constants import ANGLE_ARCSEC

    texCoords, texCoordBuf = _getCoords(inputImage.width(), inputImage.height(), 0, 1, 0, 1, lambda x: x)
    bpCoords, bpCoordBuf = _getCoords(inputImage.width(), inputImage.height(), *lrbt, 
                                           lambda x: imgPlane.traceThetaApproximately(x*ANGLE_ARCSEC)/ANGLE_ARCSEC)

    fb = QtGui.QOpenGLFramebufferObject(*outputWH)
    if not fb.bind():
        raise Exception("Coudln't bind framebuffer")

    texture = QtGui.QOpenGLTexture(QtGui.QOpenGLTexture.Target2D)
    if not texture.create():
        raise Exception("Unable to create texture")

    texture.bind()
    texture.setData(inputImage)

    prog = QtGui.QOpenGLShaderProgram()
    prog.addShaderFromSourceCode(QtGui.QOpenGLShader.Vertex,"""
        attribute vec2 a_texcoord;
        attribute vec2 a_bpcoord;
        varying vec2 v_tpos;
        uniform vec2 u_xy0, u_xy1;

        void main()
        {
            v_tpos = a_texcoord;

            vec2 diff = u_xy1-u_xy0;
            gl_Position = vec4(2.0*(a_bpcoord-u_xy0)/diff-1.0, 0, 1);
        }
        """)

    prog.addShaderFromSourceCode(QtGui.QOpenGLShader.Fragment,"""
        uniform sampler2D u_tex;
        varying vec2 v_tpos;

        void main()
        {
            vec4 c = texture2D(u_tex, v_tpos);
            gl_FragColor = c;
        }
    """)
    if not prog.link():
        raise Exception("Unable to link program")
    prog.bind()

    bpCoords = bpCoords.reshape((-1,2))
    minX, minY, maxX, maxY = bpCoords[:,0].min(), bpCoords[:,1].min(), bpCoords[:,0].max(), bpCoords[:,1].max()
    print("minMax", minX, minY, maxX, maxY)

    prog.setUniformValue("u_xy0", QtCore.QPointF(minX, minY))
    prog.setUniformValue("u_xy1", QtCore.QPointF(maxX, maxY))

    bpLoc = prog.attributeLocation("a_bpcoord")
    texLoc = prog.attributeLocation("a_texcoord")

    texCoordBuf.bind()
    prog.setAttributeBuffer(texLoc, gl.GL_FLOAT, 0, 2)
    prog.enableAttributeArray(texLoc)

    bpCoordBuf.bind()
    prog.setAttributeBuffer(bpLoc, gl.GL_FLOAT, 0, 2)
    prog.enableAttributeArray(bpLoc)

    gl.glViewport(0, 0, fb.width(), fb.height())
    gl.glClearColor(0, 0, 0, 0)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT)
    gl.glDrawArrays(gl.GL_TRIANGLES, 0, texCoords.shape[0]//2)

    gl.glFinish()
    print("glGetError =", gl.glGetError())

    img = fb.toImage()
    if img.isNull():
        raise Exception("Unable to get image from framebuffer")

    return img, [minX, minY], [maxX, maxY]

def main():

    import grale.images as images
    import grale.lenses as lenses
    import grale.plotutil as plotutil
    from grale.constants import ANGLE_ARCSEC, DIST_MPC, MASS_SUN

    app = QtWidgets.QApplication(sys.argv)

    img = QtGui.QImage(1,1, QtGui.QImage.Format_RGB32)
    img.load("/tmp/040518_LG_dark-matter-galaxy_feat.jpg")

    l = lenses.PlummerLens(1000*DIST_MPC, { "mass": 1e14*MASS_SUN, "width": 5*ANGLE_ARCSEC })
    li = plotutil.LensInfo(l, size=60*ANGLE_ARCSEC,numxy=8)
    ip = images.ImagePlane(li.getLensPlane(), 1800*DIST_MPC, 800*DIST_MPC)

    area = li.getArea()
    center = (area["bottomleft"]+area["topright"])/(2.0*ANGLE_ARCSEC)
    sizes =  (area["topright"]-area["bottomleft"])/ANGLE_ARCSEC

    #bp = BackProjectObject()
    #img = bp.backProject(ip, img, center, sizes, [800, 600])
    #del bp

    img, bl, tr = backProject(ip, img, center, sizes, [800, 600])
    img.save("testimg.png")

    img = trace(ip, img, [ (bl[0]+tr[0])/2, (bl[1]+tr[1])/2 ], [(-bl[0]+tr[0]), (-bl[1]+tr[1]) ], [1024, 768], subSample=1)
    img.save("testimg2.png")

if __name__ == "__main__":
    main()
