"""Functions to calculate transport coefficients."""

from astropy import units

from ..constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, ion_mass,
                         charge_state)

import numpy as np

from ..utils import _check_quantity


def Coulomb_logarithm():
    r"""Returns the Coulomb logarithm.

    Parameters
    ----------

    Returns
    -------
    The number 15.

    Raises
    ------

    Notes
    -----

    Examples
    --------

    See also
    --------
    Debye_number

    """

    return 15.0
