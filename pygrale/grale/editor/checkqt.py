from PyQt5 import QtWidgets, QtCore

_qtAvailFlag = [ None ]

def _qtCheckProcess(q):
    app = QtWidgets.QApplication(["appname"])
    w = QtWidgets.QMainWindow()
    t = QtCore.QTimer(w)
    t.setInterval(100)
    t.setSingleShot(True)
    t.timeout.connect(app.quit)
    t.start()
    app.exec_()
    if q:
        q.put(True)

# Test availability of Qt using another process; this prevents
# this process getting killed if e.g. DISPLAY is not set (for X11)
def isQtAvailable():

    if _qtAvailFlag[0] is None:

        import multiprocessing
        from multiprocessing import Process, Queue

        q = Queue()
        p = Process(target=_qtCheckProcess, args=(q,))
        p.start()
        p.join()
        try:
            d = False
            d = q.get(False, 1)
        except Exception as e:
            pass
        _qtAvailFlag[0] = d

    return _qtAvailFlag[0]

def checkQtAvailable():
    if isQtAvailable():
        return
    raise Exception("Interactive use with Qt does not seem to be available")

def main():
    checkQtAvailable()

if __name__ == "__main__":
    main()
