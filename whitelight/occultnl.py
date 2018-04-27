import numpy as np
import scipy.special as special
from mpmath import mp

def occultnonlin(z,p0, cn):
    """Nonlinear limb-darkening light curve; cf. Section 3 of Mandel & Agol (2002).

    :INPUTS:
        z -- sequence of positional offset values

        p0 -- planet/star radius ratio

        cn -- four-sequence. nonlinear limb darkening coefficients

    :EXAMPLE:
        ::

         # Reproduce Figure 2 of Mandel & Agol (2002):
         from pylab import *
         import transit
         z = linspace(0, 1.2, 50)
         cns = vstack((zeros(4), eye(4)))
         figure()
         for coef in cns:
             f = transit.occultnonlin(z, 0.1, coef)
             plot(z, f)

    :SEE ALSO:
       :func:`t2z`, :func:`occultnonlin_small`, :func:`occultuniform`, :func:`occultquad`

    :NOTES: 
        Scipy is much faster than mpmath for computing the Beta and
        Gauss hypergeometric functions.  However, Scipy does not have
        the Appell hypergeometric function -- the current version is
        not vectorized.

    :CREDIT: Ian Crossfield
    """
    # 2011-04-15 15:58 IJC: Created; forking from occultquad
    #import pdb

    # Initialize:
    cn0 = np.array(cn, copy=True)
    z = np.array(z, copy=True)
    F = np.ones(z.shape, float)

    p = np.abs(p0) # Save the original input


    # Test the simplest case (a zero-sized planet):
    if p==0:
        ret = np.ones(z.shape, float)
        return ret

    # Define useful constants:
    c0 = 1. - np.sum(cn0)
    # Row vectors:
    c = np.concatenate(([c0], cn0))
    n = np.arange(5, dtype=float)
    # Column vectors:
    cc = c.reshape(5, 1)
    nn = n.reshape(5,1)  
    np4 = n + 4.
    nd4 = n / 4.
    twoOmega = 0.5*c[0] + 0.4*c[1] + c[2]/3. + 2.*c[3]/7. + 0.25*c[4]

    a = (z - p)**2
    b = (z + p)**2
    am1 = a - 1.
    bma = b - a
    
    k = 0.5 * np.sqrt(-am1 / (z * p))
    p2 = p**2
    z2 = z**2


    # Define the many necessary indices for the different cases:
    i01 = (p > 0) * (z >= (1. + p))
    i02 = (p > 0) * (z > (.5 + np.abs(p - 0.5))) * (z < (1. + p))
    i03 = (p > 0) * (p < 0.5) * (z > p) * (z <= (1. - p))  # also contains Case 4
    #i04 = (z==(1. - p))
    i05 = (p > 0) * (p < 0.5) * (z == p)
    i06 = (p == 0.5) * (z == 0.5)
    i07 = (p > 0.5) * (z == p)
    i08 = (p > 0.5) * (z >= np.abs(1. - p)) * (z < p)
    i08a = (p == 1) * (z == 0)
    i09 = (p > 0) * (p < 1) * (z > 0) * (z < (0.5 - np.abs(p - 0.5)))
    i10 = (p > 0) * (p < 1) * (z == 0)
    i11 = (p > 1) * (z >= 0.) * (z < (p - 1.))

    iN = i02 + i08
    iM = i03 + i09

    # Compute N and M for the appropriate indices:
    #  (Use the slow, non-vectorized appellf1 function:)
    myappellf1 = np.frompyfunc(mp.appellf1, 6, 1)
    N = np.zeros((5, z.size), float)
    M = np.zeros((3, z.size), float)
    if iN.any():
        termN = myappellf1(0.5, 1., 0.5, 0.25*nn + 2.5, am1[iN]/a[iN], -am1[iN]/bma[iN])
        N[:, iN] = ((-am1[iN])**(0.25*nn + 1.5)) / np.sqrt(bma[iN]) * \
            special.beta(0.25*nn + 2., 0.5) * \
            (((z2[iN] - p2) / a[iN]) * termN - \
                 special.hyp2f1(0.5, 0.5, 0.25*nn + 2.5, -am1[iN]/bma[iN]))

    if iM.any():
        termM = myappellf1(0.5, -0.25*nn[1:4] - 1., 1., 1., -bma[iM]/am1[iM], -bma[iM]/a[iM]) 
        M[:, iM] = ((-am1[iM])**(0.25*nn[1:4] + 1.)) * \
            (((z2[iM] - p2)/a[iM]) * termM - \
                 special.hyp2f1(-0.25*nn[1:4] - 1., 0.5, 1., -bma[iM]/am1[iM]))


    # Begin going through all the cases:

    # Case 1:
    F[i01] = 1.

    # Case 2: (Gauss and Appell hypergeometric functions)
    F[i02] = 1. - (1. / (np.pi*twoOmega)) * \
        (N[:, i02] * cc/(nn + 4.) ).sum(0)

    # Case 3 : (Gauss and Appell hypergeometric functions)
    F[i03] = 1. - (0.5/twoOmega) * \
        (c0*p2 + 2*(M[:, i03] * cc[1:4]/(nn[1:4] + 4.)).sum(0) + \
             c[-1]*p2*(1. - 0.5*p2 - z2[i03]))

    # Case 5: (Gauss hypergeometric function)
    F[i05] = 0.5 + \
        ((c/np4) * special.hyp2f1(0.5, -nd4 - 1., 1., 4*p2)).sum() / twoOmega

    # Case 6:  Gamma function
    F[i06] = 0.5 + (1./(np.sqrt(np.pi) * twoOmega)) * \
        ((c/np4) * special.gamma(1.5 + nd4) / special.gamma(2. + nd4)).sum()

    # Case 7: Gauss hypergeometric function, beta function
    F[i07] = 0.5 + (0.5/(p * np.pi * twoOmega)) * \
        ((c/np4) * special.beta(0.5, nd4 + 2.) * \
             special.hyp2f1(0.5, 0.5, 2.5 + nd4, 0.25/p2)).sum()

    # Case 8: (Gauss and Appell hypergeometric functions)
    F[i08a] = 0.
    F[i08] =  -(1. / (np.pi*twoOmega)) * (N[:, i02] * cc/(nn + 4.) ).sum(0)

    # Case 9: (Gauss and Appell hypergeometric functions)
    F[i09] = (0.5/twoOmega) * \
        (c0 * (1. - p2) + c[-1] * (0.5 - p2*(1. - 0.5*p2 - z2[i09])) - \
             2*(M[:, i09] * cc[1:4] / (nn[1:4] + 4.)).sum(0))

    # Case 10: 
    F[i10] = (2. / twoOmega) * ((c/np4) * (1. - p2)**(nd4 + 1.)).sum()

    # Case 11:
    F[i11] = 0.


    # We're done!

    return F
