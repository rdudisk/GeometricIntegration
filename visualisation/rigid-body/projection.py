import sys
from sys import __stdout__ as stdout
from mpmath import *
import csv

def hat(w):
    return matrix([[0.0,-w[2],w[1]],[w[2],0.0,-w[0]],[-w[1],w[0],0.0]])

# parameters used for the rotation of referential and projection in 2D.
# see cahier II p. 52--54 for details
a1 = 7*pi/32
a2 = -pi/6
R1 = matrix([[cos(a1),-sin(a1),0.0],[sin(a1),cos(a1),0.0],[0.0,0.0,1.0]])
uhat = hat(R1*matrix([0,1.0,0]))
R2 = diag([1.0,1.0,1.0]) + sin(a2)*uhat + (1-cos(a2))*uhat*uhat
_P = matrix([[0,1,0],[0,0,1]])*R1.T*R2.T
_w = R2*R1*matrix([1,0,0])

# project w onto 2D plane
def proj(w):
    global _P
    return _P*w

def main(csvFilename):
    csvFile = open(csvFilename,'r')
    csvData = csv.reader(csvFile)

    fo = open('output.csv','w')
    fo.write("t,x,y,z,px,py\n")
    for e in csvData:
        # p = diag(I)*w[i]
        # el = proj(p)
        el = proj(matrix([float(e[1]),float(e[2]),float(e[3])]))
        fo.write("{0},{1},{2},{3},{4},{5}\n".format(e[0],e[1],e[2],e[3],el[0],el[1]))
    fo.close()

if __name__=="__main__":
    main(sys.argv[1])
