from grale.all import *

def getWeakData(partial = False):
    ellImg = images.ImagesData(1, shear=True, shearsigma=True, redshift=True, redshiftsigma=True)
    for l in open("./shearcatalog_redshifs.txt"):
        if l.startswith("#"):
            continue

        x, y, e1, e2, zs = map(float, l.split())
        add = True
        if partial:
            lim = 200
            if (x < -lim and y > lim) or (x > lim and y < -lim):
                add = False
        if add:
            ellImg.addPoint(0, V(x, y)*ANGLE_ARCSEC, shear=V(e1,e2), redshift=zs, redshiftsigma=0, shearsigma=[0.0,0.0])
    return ellImg

def main():
    img = getWeakData(partial=True)
    plotutil.plotImagesData(img)
    plt.gca().set_aspect("equal")
    plt.show()

if __name__ == "__main__":
    main()
