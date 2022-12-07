#!/usr/bin/env python3

from grale.all import *
from grale.privutil import PipePair
import subprocess
import socket
from threading import Thread
import fcntl
import termios
import array
import select
import os
import sys

def getBytesAvailable(x, sockSize):
    fcntl.ioctl(x, termios.FIONREAD, sockSize)
    return sockSize[0]

def writeAll(desc, data):
    bytesWritten = 0
    while bytesWritten < len(data):
        num = os.write(desc, data[bytesWritten:])
        if num < 0:
            raise Exception("Error writing!")
        bytesWritten += num

def inputToPipe(s, desc, localStopFlag, globalStopFlag):
    try:
        pollObj = select.poll()
        pollObj.register(s, select.POLLIN|select.POLLERR|select.POLLHUP)
        sockSize = array.array('i', [0])
        while not localStopFlag[0] and not globalStopFlag[0]:
            r = pollObj.poll(500)
            if not r:
                continue
            num = getBytesAvailable(s, sockSize)
            if num == 0:
                break

            data = s.recv(num)
            writeAll(desc, data)
    except Exception as e:
        print("ERROR:", e)

    print("inputToPipe thread done")

def pipeToOutput(s, desc, localStopFlag, globalStopFlag):
    try:
        pollObj = select.poll()
        pollObj.register(desc, select.POLLIN|select.POLLERR|select.POLLHUP)
        sockSize = array.array('i', [0])
        while not localStopFlag[0] and not globalStopFlag[0]:
            r = pollObj.poll(500)
            if not r:
                continue

            num = getBytesAvailable(desc, sockSize)
            if num == 0:
                break

            data = os.read(desc, num)
            if len(data) != num:
                raise Exception("Incomplete read")
            s.sendall(data)
    except Exception as e:
        print("ERROR:", e)

    print("pipeToOutput thread done")

def main():
    portNumber = int(sys.argv[1])
    mpirunCommandWithoutPipes = sys.argv[2:]

    globalStopFlag = [ False ]

    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("127.0.0.1", portNumber))
            s.listen()
            while True:
                print("Waiting for connection...")
                conn, addr = s.accept()
                localStopFlag = [ False ]
                with conn:
                    print("Got connection!")
                    pp = PipePair() # Keep it for the lifetime of the object
                    print("Using pipes", pp.rdFileName, pp.wrFileName)

                    t1 = Thread(target=inputToPipe, args=[conn, pp.wrPipeDesc, localStopFlag, globalStopFlag])
                    t2 = Thread(target=pipeToOutput, args=[conn, pp.rdPipeDesc, localStopFlag, globalStopFlag])

                    print("Starting MPI process and threads")
                    mpiProc = subprocess.Popen(mpirunCommandWithoutPipes + [ pp.wrFileName, pp.rdFileName ])
                    t1.start()
                    t2.start()

                    print("Waiting for MPI process to finish")
                    mpiProc.wait()
                    print("MPI process finished")

                    localStopFlag[0] = True
                    print("Waiting for threads to stop")
                    t1.join()
                    t2.join()
    except KeyboardInterrupt:
        print("Got keyboard interrupt")
        globalStopFlag[0] = True # Signal stop to threads
    
    except Exception as e:
        print("ERROR:", e)
        globalStopFlag[0] = True # Signal stop to threads

    print("Done")

if __name__ == "__main__":
    main()

