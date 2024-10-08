"""Fit to Deuteron Lineshape NMR Signal

Translated from original C code by C. Dulya into Python by J. Maxwell in 2021.
"A line-shape analysis for spin-1 NMR signals", C. Dulya et. al.,
SMC Collaboration, NIM A 398 (1997) 109-125.

Called as:
result = DFits(freqs, sweep, params)

where freqs is a list of frequecny points, sweep is a list of
signal magnitudes, and params is a dict of initial parameters.

"result" is then a results object in the form of lmfit.
"""

import numpy as np
from lmfit import Model

class DFits():
    '''Fit to Deuteron lineshape, translated to Python by J.Maxwell from C code by C.Dulya.
    "A line-shape analysis for spin-1 NMR signals", C. Dulya et. al., SMC Collaboration, NIM A 398 (1997) 109-125.
    '''

    def __init__(self, freqs, signal, p):
        '''Fits on signal

        Args:
            freqs: list of frequency points (X axis)
            signal: list of signal points   (Y axis)
            p: dict of initial parameters (A, G, r, wQ, wL, eta, xi)

        Returns:
            result object from lmfit
        '''
        mod = Model(self.FitFunc)
        params = mod.make_params(A=p['A'], G=p['G'], r=p['r'], wQ=p['wQ'], wL=p['wL'], eta=p['eta'], xi=p['xi'])
        self.result = mod.fit(signal, params=params, w=freqs)
        return

    def FitFunc(self, w, A, G, r, wQ, wL, eta, xi):
        '''Overall deuteron lineshape function'''
        R = (w - wL) / (3 * wQ)

        Ip, dIpdr = self.Iplus(r, wQ / wL, R)
        Im, dImdr = self.Iminus(r, wQ / wL, R)

        Fm, dFm_dR, dFm_dA, dFm_dEta = self.FandDerivs(R, A, -1, eta)
        Fp, dFp_dR, dFp_dA, dFp_dEta = self.FandDerivs(R, A, 1, eta)

        Fm /= wQ
        dFm_dR /= wQ
        dFm_dA /= wQ
        dFm_dEta /= wQ

        Fp /= wQ
        dFp_dR /= wQ
        dFp_dA /= wQ
        dFp_dEta /= wQ

        F = G * (Im * Fm + Ip * Fp)  # Lineshape
        Fm = G * (Im * Fm)  # Lineshape from minus
        Fp = G * (Ip * Fp)  # Lineshape from plus

        fAsym = 1 + 0.5 * xi * (1 + R)  # False Asymmetry xi = a[7]
        dF_dXi = 0.5 * (1 + R)
        bg = 0  # background

        y = fAsym * F + bg  # total
        ym = fAsym * Fm
        yp = fAsym * Fp
        return y

    def Iplus(self, r, Q, R):
        '''Returns: II, dI_dr '''
        r3QR = np.power(r, -3 * Q * R)
        NN = r * (r + r3QR) + 1
        II = r * (r - r3QR) / NN
        dI_dr = (2 * r * (1 - II) - (1 - 3 * Q * R) * r3QR * (1 + II)) / NN
        return II, dI_dr

    def Iminus(self, r, Q, R):
        '''Returns: II, dI_dr '''
        r3QR = np.power(r, 3 * Q * R)
        NN = r * (r + r3QR) + 1
        II = (r * r3QR - 1) / NN
        dI_dr = ((1 + 3 * Q * R) * r3QR * (1 - II) - 2 * r * II) / NN
        return II, dI_dr

    def Integrals(self, R, A, eps, Y2, etac2p):
        ''' Returns: ans1, ans2, ans3, ans4'''
        Y = np.sqrt(Y2)
        Yx2 = 2 * Y
        z2 = 1 - eps * R - etac2p
        A2 = A * A
        q4 = z2 * z2 + A2
        q2 = np.sqrt(q4)
        qq = np.sqrt(q2)

        cosa = z2 / q2
        cosa_2 = 1 / np.sqrt(2) * np.sqrt(1 + cosa)
        sina_2 = 1 / np.sqrt(2) * np.sqrt(1 - cosa)

        fTmp = Y2 + q2
        fVal = Yx2 * qq * cosa_2

        La = 0.5 * sina_2 * np.log((fTmp + fVal) / (fTmp - fVal))
        Ta = cosa_2 * (np.pi / 2 + np.arctan((Y2 - q2) / (Yx2 * qq * sina_2)))
        Arg = (Y2 * (Y2 - 2 * z2) + q4)

        ans1 = (Ta + La) / (2 * qq * A)
        ans2 = (Ta - La) * qq / (2 * A)
        ans3 = z2 * (ans2) + (2 * A2 + q4) * (ans1) + (Y / Arg) * (Y2 * z2 + 2 * A2 - q4)
        ans4 = ((Y / Arg) * (Y2 - z2) + z2 * (ans1) + (ans2)) / (4 * A2)

        return ans1, ans2, ans3, ans4

    def FandDerivs(self, R, A, eps, eta):
        '''Returns FF, dFdA, dFdR, dFdEta'''
        if eta < 0.001:
            Y2 = 3
            I1, I2, I3, I4 = self.Integrals(R, A, eps, Y2, 0)
            FF = I1 * A
            dFdA = (I1 - 2.0 * A * A * I3)
            dFdR = ((1 - eps * R) * I3 - I4) * 2 * A * eps
            dFdEta = 0
        else:
            Y2 = 3
            I1, I2, I3, I4 = self.Integrals(R, A, eps, Y2, 0)
            FF, dFdA, dFdR, dFdEta = (0, 0, 0, 0)
            eRm1 = 1 - eps * R
            dphi = 1

            for i in (0, 1):
                c2p = np.cos(np.pi * dphi * i)
                ec2p = eta * c2p
                Y2 = 3 - ec2p
                Y = np.sqrt(Y2)
                z2 = eRm1 - ec2p

                I1, I2, I3, I4 = self.Integrals(R, A, eps, Y2, 0)

                fac = 0.5 * np.sqrt(3) / Y
                FF += fac * I1 * A

                dFdA += fac * (I1 - 2 * A * A * I3)
                dFdR += fac * (z2 * I3 - I4) * 2 * A * eps
                gY = Y2 * (Y2 - 2 * z2) + A * A + z2 * z2
                dFdEta += 2 * A * c2p * fac * (z2 * I3 - I4 + I1 / (4 * Y2) - 1 / (4 * Y * gY))

            order = 5
            for N in [np.power(2, n) for n in range(2, order + 1)]:
                dphi = 1 / N

                for i in range(N - 1, 0, -2):
                    c2p = np.cos(np.pi * dphi * i)
                    ec2p = eta * c2p

                    Y2 = 3 - ec2p
                    Y = np.sqrt(Y2)
                    z2 = eRm1 - ec2p

                    I1, I2, I3, I4 = self.Integrals(R, A, eps, Y2, ec2p)

                    fac = np.sqrt(3) / Y

                    FF += fac * I1 * A

                    dFdA += fac * (I1 - 2 * A * A * I3)
                    dFdR += fac * (z2 * I3 - I4) * 2 * A * eps
                    gY = Y2 * (Y2 - 2 * z2) + A * A + z2 * z2
                    dFdEta += 2 * A * c2p * fac * (z2 * I3 - I4 + I1 / (4 * Y2) - 1 / (4 * Y * gY))

            FF = dphi * FF
            dFdA = dphi * dFdA
            dFdR = dphi * dFdR
            dFdEta = dphi * dFdEta

        return FF, dFdA, dFdR, dFdEta
