from JacobiBibliothek import *
import numpy as np
from sympy import pprint, Matrix, simplify, sqrt, solve
q1, q2, q3, l1, l2, l3 = symbols("q1 q2 q3 l1 l2 l3")

_0T1 = Cdh(np.array(["q1", -90]), l1, 0, -90)
_1T2 = Cdh(np.array(["q2", -90]), l2, 0, -90)
_2T3 = Cdh(np.array(["q3", 0]), 0, l3, 0)

print("##OT1##")
pprint(_0T1.getTrans())
print("##1T2##")
pprint(_1T2.getTrans())
print("##2T3##")
pprint(_2T3.getTrans())

_0T3 = _0T1.getTrans() * _1T2.getTrans() * _2T3.getTrans()
_0T3 = simplify(_0T3)
print("##0T3##")
pprint(_0T3)

print("################ Aufgabe 3.3 ###############")

print("##OT1##")
pprint(_0T1.getTrans())
_0T2 = _0T1.getTrans() * _1T2.getTrans()
_0T2 = simplify(_0T2)
print("##0T2##")
pprint(_0T2)
print("##0T3##")
pprint(_0T3)

trans1, Or1 = calcJacobiRot(_0T3, "", Matrix([0,0,0,1]), 1)
trans2, Or2 = calcJacobiRot(_0T3, _0T1.getTrans(), Matrix([0,0,0,1]), 0)
trans3, Or3 = calcJacobiRot(_0T3, _0T2, Matrix([0,0,0,1]), 0)
print("##Jacobi##")
Jacobi = Matrix([ [trans1, trans2, trans3],
                  [Or1, Or2, Or3]])
pprint(simplify(Jacobi))
JacobiDet = Jacobi.copy()
JacobiDet.row_del(5)
JacobiDet.row_del(4)
JacobiDet.row_del(3)
JacobiDet = simplify(JacobiDet)
JacobiDet = JacobiDet.det()
print(simplify(JacobiDet))
JacobiDet = JacobiDet.subs([(l1, 2), (l2, 1), (l3, 2)])
print(solve(JacobiDet, [q1, q2, q3]))

print("####Jacobi in Nullage###")
Jacobi = Jacobi.subs([(l1, 2), (l2, 1), (l3, 2), (q1, 0), (q2, 0), (q3, 0)])
pprint(Jacobi)

print("############# Aufgabe 3.4 ##########")

_0T1 = Cdh(np.array(["q1", -90]), l1, 0, -90)
_1T2 = Cdh(np.array(["q2", -90]), l2, 0, -90)
_2T3 = Cdh(np.array(["q3", 0]), 0, l3, -90)
_3T4 = Cdh(np.array(["q4", 0]), 0, l3, 0)

print("##OT1##")
pprint(_0T1.getTrans())
print("##1T2##")
pprint(_1T2.getTrans())
print("##2T3##")
pprint(_2T3.getTrans())

_0T4 = _0T1.getTrans() * _1T2.getTrans() * _2T3.getTrans() * _3T4.getTrans()
_0T4 = simplify(_0T4)
print("##0T4##")
pprint(_0T4)

_0T2 = _0T1.getTrans() * _1T2.getTrans()
_0T2 = simplify(_0T2)
_0T3 = simplify(_0T2*_2T3.getTrans())

trans1, Or1 = calcJacobiRot(_0T4, "", Matrix([0,0,0,1]), 1)
trans2, Or2 = calcJacobiRot(_0T4, _0T1.getTrans(), Matrix([0,0,0,1]), 0)
trans3, Or3 = calcJacobiRot(_0T4, _0T2, Matrix([0,0,0,1]), 0)
trans4, Or4 = calcJacobiRot(_0T4, _0T3, Matrix([0,0,0,1]), 0)

print("##Jacobi##")
Jacobi = Matrix([ [trans1, trans2, trans3, trans4],
                  [Or1, Or2, Or3, Or4]])
pprint(simplify(Jacobi))
print("####Jacobi in Nullage###")
Jacobi = Jacobi.subs([(l1, 2), (l2, 1), (l3, 2), (q1, 0), (q2, 0), (q3, 0)])
pprint(Jacobi)