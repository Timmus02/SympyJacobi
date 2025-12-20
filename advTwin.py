from JacobiBibliothek import *
import numpy as np
from sympy import pprint, Matrix, simplify

q1, q2, q3, q4, l2 = symbols("q1 q2 q3 q4 l2")
l2 = "l2"

l1 = 100
l3 = 600
l4 = 600
l5 = 600
l6 = 300

_0T1 = EulerABC(0, 0, 0, Matrix([0,0, l2]))
_1T2 = EulerABC("q1", 0, 0, Matrix([0, l1+l3, 0]))
_2T3 = EulerABC("q2", 0, 0, Matrix([0, l4, 0]))
_3T4 = EulerABC("q3", 90, 0, Matrix([-l5, 0, 0]))
_4T5 = EulerABC("q4", 0, 0, Matrix([0, 0, -l6]))

_0T5_ABC = _0T1.getTrans() * _1T2.getTrans() * _2T3.getTrans() * _3T4.getTrans() * _4T5.getTrans()
print("#####_0T5 nach ABC##########")
pprint(simplify(_0T5_ABC))
_0T5_ABC = _0T5_ABC.subs([(l2,0), (q1,0), (q2, convertToRad(-90)), (q3,0),(q4,0)])

_0T1 = Cdh(np.array([0, 0]), l2, l1+l3, 0)
_1T2 = Cdh(np.array(["-q1", 0]), 0, l4, 0)
_2T3 = Cdh(np.array(["-q2", 0]), 0, l5, 0)
_3T4 = Cdh(np.array(["-q3", 90]), 0, 0, 90)
_4T5 = Cdh(np.array(["q4", 0]), l6, 0, 0)
_0T5 = _0T1.getTrans() * _1T2.getTrans() * _2T3.getTrans() * _3T4.getTrans() * _4T5.getTrans()
print("#####_0T5 nach DH###########")
pprint(simplify(_0T5))
_0T5 = _0T5.subs([(l2,0), (q1,0), (q2,0), (q3,0),(q4,0)])

print("#########_0T5 ABC in Nullage DH############")
pprint(simplify(_0T5_ABC))
print("#####_0T5 nach DH in Nullage###########")
pprint(simplify(_0T5))