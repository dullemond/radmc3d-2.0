def bplanck(freq,temp):
    """
    This function computes the Planck function

                   2 h nu^3 / c^2
       B_nu(T)  = ------------------    [ erg / cm^2 s ster Hz ]
                  exp(h nu / kT) - 1

    Arguments:
         freq  [Hz]            = Frequency in Herz
         temp  [K]             = Temperature in Kelvin
    """
    const1 = 4.7991598e-11
    const2 = 1.4745284e-47
    x      = const1*freq/(temp+1e-99)
    xp     = np.exp(-x)
    xxp    = 1-xp
    mask   = x<1e-6
    xxp[mask] = x[mask]
    bpl    = xp*const2*(freq**3)/xxp
    return bpl

