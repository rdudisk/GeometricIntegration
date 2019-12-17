import sys
from sympy import *
from sympy.diffgeom import *

h,q,w = symbols('h q w')
m,I = symbols('m I')

class lag(Function):
    @classmethod
    def eval(cls,x):
        return 0.5*x

r = lag(w)
print(r)
print(lag(w).diff(w))

class hat(Function):
    @classmethod
    def eval(cls,w):
        return Matrix(([0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]))

class expso3(Function):
    @classmethod
    def eval(cls,w):
        theta = w.norm()
        K = hat(w/theta)
        return eye(3)+sin(theta)*K+(1.0-cos(theta))*K*K

class cayso3(Function):
    @classmethod
    def eval(cls,w):
        return eye(3)+(4.0/4.0+w.norm())*(hat(w)+hat(w)*hat(w)/2.0)

class cayso3alginv(Function):
    @classmethod
    def eval(cls,w):
        return (eye(3)-hat(w)/2.0).inv()*(eye(3)+hat(w)/2.0)

print(cayso3(Matrix([1,0,0])))

M = Manifold('M',3)
p = Patch('P',M)
rect = CoordSystem('rect',p)
x = BaseScalarField(rect,0)
y = BaseScalarField(rect,1)
z = BaseScalarField(rect,2)
e_x = BaseVectorField(rect,0)
e_y = BaseVectorField(rect,1)
e_z = BaseVectorField(rect,2)

field = cayso3(Matrix([x,y,z]))
dcay = Differential(field)

el = cayso3alginv(Matrix([x,y,z]))
pprint(el)

class dRcay(Differential):
    @classmethod
    def __call__(self,*v):
        return dcay(v)*el

help(dcay)
help(dRcay)
pprint(dRcay)
pprint(dRcay.__call__(e_x))
# pprint(dcay)
# pprint(dcay(e_x))
