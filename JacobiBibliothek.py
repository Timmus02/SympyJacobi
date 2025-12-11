from sympy import symbols, Matrix, cos, sin, pprint, trigsimp, pi
import numpy as np

class Cdh:
    phi = np.empty(2)
    alpha = np.empty(2)
    a = 0
    d = 0
    phi, d, a, alpha = symbols('phi d a alpha')

    def __init__(self, _phi:np.ndarray, _d, _a, _alpha: np.ndarray): #array of phi and alph should be ["value or symbol", Offset[int]]
        self.alpha = _alpha
        self.a = _a
        self.d = _d
        self.phi = _phi
        self.calcTrans()

    def calcTrans(self):
        phi_offset, alpha_offset = symbols("phi_offset alpha_offset")
        phi, phi_offset, alpha, alpha_offset, a, d = self.phi[0], self.phi[1], self.alpha[0], self.alpha[1] , self.a, self.d

        phi_offset = self._convertToRad(phi_offset)
        alpha_offset = self._convertToRad(alpha_offset)

        if not type(phi) == np.str_:
            phi = self._convertToRad(phi)
        if not type(alpha) == np.str_:
            alpha = self._convertToRad(alpha)

        cPhi = self._AdditionstheoremCOS(phi, phi_offset)
        sPhi = self._AdditionstheoremSIN(phi, phi_offset)
        cAlpha = self._AdditionstheoremCOS(alpha, alpha_offset)
        sAlpha = self._AdditionstheoremSIN(alpha, alpha_offset)
        self.tran = Matrix([
            [cPhi, -sPhi*cAlpha,  sPhi*sAlpha,  a*cPhi],
            [sPhi,  cPhi*cAlpha, -cPhi*sAlpha,  a*sPhi],
            [0,         sAlpha,           cAlpha,           d],
            [0,         0,                    0,                    1]
        ])
 
    def getTrans(self):
        return self.tran
    def _AdditionstheoremCOS(self, x1, x2):
        return cos(x1) * cos(x2) + -1 * sin(x1) * sin(x2)
    def _AdditionstheoremSIN(self, x1, x2):
        return sin(x1) * cos(x2) + cos(x1) * sin(x2)
    def _convertToRad(self, deg):
        print(type(deg))
        if type(deg) == np.float16 or type(deg) == np.int64:
            return (pi/180)*deg
        if type(deg) == np.str_:
            if deg.isdigit():
               return (pi/180)*float(deg) 
        print("Wert nicht convertierbar in RAD")
        exit()

_0T1 = Cdh(np.array(["q1", 0]), 2000, 0, np.array([90, 0]))
pprint(_0T1.getTrans())
_0T1.phi = np.array([0, 0])
_0T1.calcTrans()
pprint(_0T1.getTrans())
