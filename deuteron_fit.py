"""Fit to Deuteron Lineshape NMR Signal

Translated from original C code by C. Dulya into Python by J. Maxwell in 2021.
"A line-shape analysis for spin-1 NMR signals", C. Dulya et. al.,
SMC Collaboration, NIM A 398 (1997) 109-125.

Called as:
result = fit_deuteron(freqs, sweep, params)

where freqs is a list of frequency points, sweep is a list of
signal magnitudes, and params is a dict of initial parameters
(A, G, r, wQ, wL, eta, xi).

"result" is then a results object in the form of lmfit.
"""

import numpy as np
from lmfit import Model
from scipy.signal import hilbert


def fit_deuteron(freqs, signal, initial_params):
    """Fit to deuteron lineshape.

    Args:
        freqs: list of frequency points (X axis)
        signal: list of signal points   (Y axis)
        initial_params: dict of initial parameters (A, G, r, wQ, wL, eta, xi)

    Returns:
        result object from lmfit
    """
    mod = Model(_fit_func)
    params = mod.make_params(**initial_params)
    return mod.fit(signal, params=params, w=freqs)


def _fit_func(w, A, G, r, wQ, wL, eta, xi):
    """Overall deuteron lineshape function."""
    R = (w - wL) / (3 * wQ)

    Ip, _ = _iplus(r, wQ / wL, R)
    Im, _ = _iminus(r, wQ / wL, R)

    Fm, *_ = _f_and_derivs(R, A, -1, eta)
    Fp, *_ = _f_and_derivs(R, A, 1, eta)

    Fm /= wQ
    Fp /= wQ

    F = G * (Im * Fm + Ip * Fp)
    fAsym = 1 + 0.5 * xi * (1 + R)
    bg = 0

    return fAsym * F + bg


def _iplus(r, Q, R):
    """Returns: II, dI_dr"""
    r3QR = r ** (-3 * Q * R)
    NN = r * (r + r3QR) + 1
    II = r * (r - r3QR) / NN
    dI_dr = (2 * r * (1 - II) - (1 - 3 * Q * R) * r3QR * (1 + II)) / NN
    return II, dI_dr


def _iminus(r, Q, R):
    """Returns: II, dI_dr"""
    r3QR = r ** (3 * Q * R)
    NN = r * (r + r3QR) + 1
    II = (r * r3QR - 1) / NN
    dI_dr = ((1 + 3 * Q * R) * r3QR * (1 - II) - 2 * r * II) / NN
    return II, dI_dr


def _integrals(R, A, eps, Y2, etac2p):
    """Returns: ans1, ans2, ans3, ans4"""
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


def _f_and_derivs(R, A, eps, eta):
    """Returns FF, dFdA, dFdR, dFdEta"""
    if eta < 0.001:
        I1, _, I3, I4 = _integrals(R, A, eps, 3, 0)
        FF = I1 * A
        dFdA = I1 - 2.0 * A * A * I3
        dFdR = ((1 - eps * R) * I3 - I4) * 2 * A * eps
        dFdEta = 0
    else:
        FF, dFdA, dFdR, dFdEta = 0, 0, 0, 0
        eRm1 = 1 - eps * R
        dphi = 1

        for i in (0, 1):
            c2p = np.cos(np.pi * dphi * i)
            ec2p = eta * c2p
            Y2 = 3 - ec2p
            Y = np.sqrt(Y2)
            z2 = eRm1 - ec2p

            I1, _, I3, I4 = _integrals(R, A, eps, Y2, 0)

            fac = 0.5 * np.sqrt(3) / Y
            FF += fac * I1 * A
            dFdA += fac * (I1 - 2 * A * A * I3)
            dFdR += fac * (z2 * I3 - I4) * 2 * A * eps
            gY = Y2 * (Y2 - 2 * z2) + A * A + z2 * z2
            dFdEta += 2 * A * c2p * fac * (z2 * I3 - I4 + I1 / (4 * Y2) - 1 / (4 * Y * gY))

        order = 5
        for N in [2 ** n for n in range(2, order + 1)]:
            dphi = 1 / N

            for i in range(N - 1, 0, -2):
                c2p = np.cos(np.pi * dphi * i)
                ec2p = eta * c2p

                Y2 = 3 - ec2p
                Y = np.sqrt(Y2)
                z2 = eRm1 - ec2p

                I1, _, I3, I4 = _integrals(R, A, eps, Y2, ec2p)

                fac = np.sqrt(3) / Y
                FF += fac * I1 * A
                dFdA += fac * (I1 - 2 * A * A * I3)
                dFdR += fac * (z2 * I3 - I4) * 2 * A * eps
                gY = Y2 * (Y2 - 2 * z2) + A * A + z2 * z2
                dFdEta += 2 * A * c2p * fac * (z2 * I3 - I4 + I1 / (4 * Y2) - 1 / (4 * Y * gY))

        FF *= dphi
        dFdA *= dphi
        dFdR *= dphi
        dFdEta *= dphi

    return FF, dFdA, dFdR, dFdEta


def fit_deuteron_complex(freqs, signal, initial_params):
    """Fit to complex deuteron lineshape (absorption + dispersion mixed by phase).

    Args:
        freqs: list of frequency points (X axis)
        signal: list of signal points   (Y axis)
        initial_params: dict of initial parameters (A, G, r, wQ, wL, eta, xi, phase)

    Returns:
        result object from lmfit
    """
    mod = Model(_fit_func_complex)
    params = mod.make_params(**initial_params)
    return mod.fit(signal, params=params, w=freqs)


def _fit_func_complex(w, A, G, r, wQ, wL, eta, xi, phase):
    """Deuteron lineshape with phase rotation between absorption and dispersion."""
    absorption = _fit_func(w, A, G, r, wQ, wL, eta, xi)
    dispersion = np.imag(hilbert(absorption))
    return np.cos(phase) * absorption + np.sin(phase) * dispersion
