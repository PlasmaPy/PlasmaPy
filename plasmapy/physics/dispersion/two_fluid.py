# FUNCTION TO CALCULATE PHASE SPEEDS OF THE THREE BRANCHES OF TWO FLUID
# DISPERSION RELATION (e.g. STRINGER JPP 1963, ROGERS PRL 2001,
# and BELLAN JGR 2012)

#
#                   Tulasi Nandan Parashar

import matplotlib.pyplot as plt
import numpy as np


def tfps(beta=0.6, speed=1., de2=0.000545, theta=0, wavenumber=1):
    """
    beta: Total plasma beta,
    speed: Alfven speed based on mean field,
    de2: me/mi,
    theta: Angle of propagation in degrees
    wavenumber: Wavenumber of interest in units of kdi

    Output is frequencies of the roots and the phase speeds w/k
    The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron
    """

    pi_angle = theta * (np.pi / 180.0)
    sine = np.sin(pi_angle)
    cosine = np.cos(pi_angle)
    tangent = sine / cosine
    cs = np.sqrt(beta / 2.0) * speed
    di = 1
    caksq = (cosine ** 2) * (speed ** 2)
    cmsq = (speed ** 2) + (cs ** 2)
    D = 1 + wavenumber ** 2 * de2

    term1 = -(speed ** 2 / D + cs ** 2 + caksq * (1 + wavenumber ** 2 * di ** 2 / D) / D) * wavenumber ** 2
    term2 = caksq * wavenumber ** 4 * (speed ** 2 / D + cs ** 2 + cs ** 2 * (1 + wavenumber ** 2 * di ** 2 / D)) / D
    term3 = -(cs ** 2 * wavenumber ** 6 * caksq ** 2) / D ** 2

    dispersion_relation = (1, term1, term2, term3)

    roots = np.roots(dispersion_relation)
    sq_roots = np.sqrt(roots)

    phase_speeds = sq_roots / wavenumber
    return sq_roots, phase_speeds


def single_theta_value(beta=0.6, ca=1., de2=0.000545, theta=0., kmin=1e-2, kmax=10., npoints=200):
    """
    Compute the dispersion relation w(k) vs k for a single theta value.

    beta: Total plasma beta,
    ca: Alfven speed based on mean field,
    de2: me/mi,
    theta: Angle of propagation in degrees
    kkmin: Minimum Wavenumber of interest in units of kdi
    kkmax: Maximum wavenumber of interest in units of kdi

    Output is an array with 4 columns, k, w-fast, w-alf, w-slow
    """

    kmmn = np.log10(kmin)
    kmmx = np.log10(kmax)
    wavenumber = np.logspace(kmmn, kmmx, npoints)
    warray = np.zeros((3, npoints))
    for i in range(0, npoints):
        f, s = tfps(beta, ca, de2, theta, wavenumber[i])
        warray[:, i] = f
    plt.loglog(wavenumber, warray[0, :], label='Fast/Magnetosonic')
    plt.loglog(wavenumber, warray[1, :], label='Alfven/KAW')
    plt.loglog(wavenumber, warray[2, :], label='Slow')
    plt.xlabel('$kd_i$')
    plt.ylabel('$\omega/\omega_{ci}$')
    plt.legend(loc='best', fancybox=True, framealpha=0.2)
    plt.title('Dispersion Relation for beta=' + str(beta) + ' and me/mi=' + str(de2))
    plt.show()


class TFEV:
    def __init__(self, beta=0.6, ca=1., de2=0.000545, theta=0., k=1., aa=0.1):
        """
        beta: Total plasma beta,
        ca: Alfven speed based on mean field,
        de2: me/mi,
        theta: Angle of propagation in degrees
        k: wavenumber of interest in units of kdi

        Output: Prints the two fluid eigenvector to the screen
        """

        self.beta = beta
        self.ca = ca
        self.de2 = de2
        self.theta = theta
        self.k = k
        self.aa = aa
        self.f, _ = tfps(self.beta, ca, de2, theta, k)

    def amp(self, w):
        th = self.theta * np.pi / 180.
        bb = 1 - w ** 2 * (1 + self.de2 * self.k ** 2) / (self.k ** 2 * np.cos(th) ** 2)
        sk = 'sin(' + str(round(self.k, 3)) + 'x)'
        ck = 'cos(' + str(round(self.k, 3)) + 'x)'

        def st(a):
            return str(round(a, 3))

        return {'bx': 0.,
                'by': st(round(2 * self.aa, 3)) + ck,
                'bz': st(-np.cos(th) * bb * 2 * self.aa / w) + sk,
                'ux': st(self.aa * self.k * bb * np.sin(2 * th) / (w ** 2 - self.beta * self.k ** 2)) + sk,
                'uy': st(-2 * self.aa * self.k * np.cos(th) / w) + ck,
                'uz': st(2 * self.aa * self.k * bb * np.cos(th) ** 2/w ** 2) + sk,
                'n': st((self.k * np.cos(th) / w) * self.aa * self.k * bb * np.sin(2 * th) / (w ** 2 - self.beta * self.k ** 2)) + sk}

    # The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron

    @property
    def whistler(self):
        return self.amp(self.f[0])

    @property
    def kaw(self):
        return self.amp(self.f[1])

    @property
    def cyclotron(self):
        return self.amp(self.f[2])
