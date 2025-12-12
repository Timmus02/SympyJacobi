from JacobiBibliothek import *
import numpy as np
from sympy import pprint, Matrix, simplify

_0T1 = Cdh(np.array(["q1", 0]), 2000, 0, 90)
_1T2 = Cdh(np.array(["q2", 90]), 0, 0, 90)
_2T3 = Cdh(np.array([0, 0]), "d3", 0, 0)

#_0T1 = Cdh(np.array([0, 0]), 2000, 0, 90)
#_1T2 = Cdh(np.array([0, 90]), 0, 0, 90)
#_2T3 = Cdh(np.array([0, 0]), 1000, 0, 0)

_0T2 = _0T1.getTrans() * _1T2.getTrans()
_0T3 = _0T2 * _2T3.getTrans()
print("######0T1#######")
pprint(_0T1.getTrans())
print("######0T2#######")
pprint(_0T2)
print("######0T3#######")
pprint(_0T3)

trans1, Or1 = calcJacobiRot(_0T3, "", Matrix([0,0,0,1]), 1)
trans2, Or2 = calcJacobiRot(_0T3, _0T1.getTrans(), Matrix([0,0,0,1]), 0)
trans3, Or3 = calcJacobiTrans(_0T2, 0)

print("######Jacobi#######")
Jacobi = Matrix([ [trans1, trans2, trans3],
                  [Or1, Or2, Or3]])
pprint(simplify(Jacobi))