"""
Objects representing stacks of film and/or filter layers for charged particle
detectors.
"""

__all__ = [
    "Stack",
    "Layer",
]

import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d



class Layer:
    def __init__(self, thickness: u.m, energy_axis: u.J, stopping_power: u.MeV/u.cm, active: bool = True, name: str =""):
        r"""
        A layer in a detector film stack. The layer could either be an active
        layer (the actual film medium) or an inum_active layer (a filter or
        inum_active part of the film, such as a substrate.)

        Tabulated stopping powers for protons and electrons can be found in the
        `NIST PSTAR database
        <https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html>`_
        and the
        `NIST ESTAR database
        <https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html>`_.

        Parameters
        ----------

        thickness : `~astropy.units.Quantity`
            The thickness of the layer, in units convertible to meters.

        energy_axis : `~astropy.units.Quantity`
            The energies corresponding to the stopping power array.

        stopping_power : `~astropy.units.Quantity`
            The stopping power in the material, multiplied by the material
            mass density. The stopping power is tabulated in units of
            MeV cm\ :sup:`^2` / g, so this variable has units of MeV/cm.

        active : `bool`, optional
            If `True`, this layer is marked as an active layer. The default is `True`.

        name : `str`, optional
            An optional name for the layer.

        """
        self.thickness = thickness
        self.energy_axis = energy_axis
        self.stopping_power = stopping_power
        self.active = active
        self.name = name


class Stack:
    r"""
    An ordered list of
    `~plasmapy.diagnostics.charged_particle_radiography.Layer` objects.


    Parameters
    ----------

    layers : list of `~plasmapy.diagnostics.charged_particle_radiography.Layer` objects
        A list of the `~plasmapy.diagnostics.charged_particle_radiography.Layer` objects that make up the film stack.
    """

    def __init__(self, layers):
        self._layers = layers
        self._energy_bands = None

    @property
    def num_layers(self):
        r"""
        The number of layers in the stack.
        """

        return len(self._layers)

    @property
    def num_active(self):
        r"""
        The number of layers in the stack marked 'active'
        """

        return len([layer for layer in self._layers if layer.active])

    @property
    def thickness(self):
        r"""
        The total thickness of the stack.
        """
        thickness = np.array([layer.thickness.to(u.mm).value for layer in self._layers])
        return np.sum(thickness) * u.mm

    def deposition_curves(self, energies: u.MeV, dx=1 * u.um, return_only_active=True):
        """
        Calculates the deposition of an ensemble of particles over a range of
        energies in a stack of films and filters.

        Parameters
        ----------

        stack : list of `~plasmapy.diagnostics.charged_particle_radiography.Layer` objects
            This list of `~plasmapy.diagnostics.charged_particle_radiography.Layer`
            objects defines the composition of the film stack.

        energies : `~astropy.units.Quantity` array, shape [nenergies,]
            Energies axis over which to calculate the deposition. Units convertible
            to eV.

        dx : `~astropy.units.Quantity`, optional
            The spatial resolution of the numerical integration of the stopping power.
            Defaults to 1 um.

        return_only_active : `bool`, optional
            If `True`, only the deposition in layers in which the active property
            is `True` will be returned. This is usually desirable, since particles
            captured in other layers will not be measured. If `False`, deposition in
            all layers of the stack are returned. The default is `True`.

        Returns
        -------

        deposited : `~numpy.ndarray`, shape [num_layers, nenergies]
            The fraction of an ensemble of each energy that will be deposited in
            each layer of film. The array is normalized such that the sum
            along the first dimension (all of the layers) for each population
            is unity.

        """

        energies = energies.to(u.MeV).value

        # Deposited energy in MeV
        deposited = np.zeros([len(self._layers), energies.size])

        for i, layer in enumerate(self._layers):

            # Interpolate stopping power for each energy
            # stopping power here is in MeV/cm
            sp_fcn = interp1d(
                layer.energy_axis.to(u.MeV).value,
                layer.stopping_power.to(u.MeV / u.cm).value,
                fill_value=(0, np.inf),
                bounds_error=False,
            )

            # Slice the layer into sublayer dx thick
            nsublayers = int(
                np.floor(layer.thickness.to(u.um).value / dx.to(u.um).value)
            )
            sublayers = np.ones(nsublayers) * dx.to(u.um)
            # Include any remainder in the last sublayer
            sublayers[-1] += layer.thickness.to(u.um) % dx.to(u.um)

            # Calculate the energy deposited in each sublayer
            # This is essentially numerically integrating the stopping power
            for ds in sublayers:
                # Interpolate the stopping power at the current energies
                interpolated_stopping_power = sp_fcn(energies)
                # dE is in MeV
                dE = interpolated_stopping_power * ds.to(u.cm).value

                # If dE > E for a given energy, set dE=E (stop the particle)
                dE = np.where(dE > energies, energies, dE)

                energies += -dE
                deposited[i, :] += dE

        # Normalize the deposited energy array so that each number represents
        # the fraction of a population of particles of that energy stopped
        # in that layer.
        deposited /= np.sum(deposited, axis=0)

        # If this flag is set, return only the layers that correspond to active
        # medium, ignoring the filter and substrate layers
        if return_only_active:
            active_ind = [i for i in range(len(self._layers)) if self._layers[i].active]
            deposited = deposited[active_ind, :]

        return deposited

    def energy_bands(self, energy_range, dE, dx=1 * u.um, return_only_active=True):
        """
        Calculate the energy bands in each of the active layers of a film stack.

        Parameters
        ----------
        stack : list of `~plasmapy.diagnostics.charged_particle_radiography.Layer` objects
            This list of `~plasmapy.diagnostics.charged_particle_radiography.Layer`
            objects defines the composition of the film stack.

        energy_range : `~astropy.units.Quantity` list, shape [2,]
            A range of energies to include in the calculation. Units convertible
            to eV.

        dE :  `~astropy.units.Quantity`
            Spacing between energy bins in the calculation. Units convertible
            to eV.

        dx : `~astropy.units.Quantity`, optional
            The spatial resolution of the numerical integration of the stopping power.
            Passed directly to the `~deposition_curves` method. Defaults to 1 um.

        return_only_active : `bool`, optional
            If `True`, only the energy bands of layers in which the active property
            is `True` will be returned. This is usually desirable, since particles
            captured in other layers will not be measured. If `False`, energy bands in
            all layers of the stack are returned. The default is `True`.

        Returns
        -------

        energy_bands : `~astropy.units.Quantity`, shape [num_layers, 2]
            The full-width-half-max energy range of the Bragg peak in each
            active layer of the film stack, in MeV.

        """

        energies = (
            np.arange(
                *energy_range.to(u.MeV).value,
                dE.to(u.MeV).value,
            )
            * u.MeV
        )

        deposited = self.deposition_curves(
            energies, return_only_active=return_only_active
        )

        energy_bands = np.zeros([deposited.shape[0], 2]) * u.MeV

        for i in range(deposited.shape[0]):
            bragg_curve = deposited[i, :]

            # Find the indices corresponding to half the maximum value
            # on either side of the peak
            halfmax = np.max(bragg_curve) / 2

            inds = np.argwhere(bragg_curve > halfmax)
            # Store those energies

            energy_bands[i, 0] = energies[inds[0][0]]
            energy_bands[i, 1] = energies[inds[-1][0]]

        self._energy_bands = energy_bands

        return energy_bands