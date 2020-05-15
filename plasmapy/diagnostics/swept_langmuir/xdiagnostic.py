__all__ = ["XSweptLangmuirDiagnostic"]

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from plasmapy.analysis.swept_langmuir.floating_potential import find_floating_potential
from plasmapy.diagnostics import XAbstractDiagnostic
from plasmapy.diagnostics.swept_langmuir import SweptLangmuirProbe
from warnings import warn


class XSweptLangmuirDiagnostic(XAbstractDiagnostic):
    def __init__(self, dataset: xr.Dataset):
        super().__init__(dataset)
        self._probe_class = SweptLangmuirProbe
        # self._check_dataset()

    def _check_dataset(self):
        # check for XArray Variables
        if 'voltage' not in self._ds.data_vars.keys():
            raise ValueError(f"Dataset does not have the 'voltage' Variable.")
        if 'current' not in self._ds.data_vars.keys():
            raise ValueError(f"Dataset does not have the 'current' Variable.")

    @property
    def floating_potential(self):
        try:
            return self._ds.floating_potential[0].data[()]
        except AttributeError:
            warn('Floating potential has not be calculated.')
            return np.nan

    def analyze(self):
        self.find_floating_potential()

    def find_floating_potential(self, **kwargs):
        self._check_dataset()

        # find floating potential
        vf, vf_err, fit = find_floating_potential(self._ds.voltage.data,
                                                  self._ds.current.data,
                                                  **kwargs)
        # update Dataset
        self._ds['floating_potential'] = xr.DataArray(
            [vf, vf_err],
            dims=('vf',),
            coords={'vf': ['value', 'err']}
        )
        self._ds['floating_potential_fit'] = xr.DataArray(
            [[fit['slope'], fit['slope_err']],
             [fit['intercept'], fit['intercept']],
             [fit['indices'], np.nan]],
            dims=('vf_fit_params', 'vf_fit_values'),
            coords={'vf_fit_params': ['slope', 'intercept', 'indices'],
                    'vf_fit_values': ['value', 'error']}
        )

        return vf, vf_err, fit

    #     def plot(what='all'):
    #         plot_methods = {
    #             'all': None,
    #             'floating_potential': self.plot_floating_potential,
    #         }
    #         try:
    #             plot_methods[what]()
    #         except KeyError:
    #             pass

    def plot_floating_potential_fit(self):
        self._check_dataset()
        fig, ax = plt.subplots()

        # get data subset
        sub = self._ds.floating_potential_fit.sel(vf_fit_params='indices',
                                                  vf_fit_values='value')
        Vsub = self._ds.voltage[sub].data
        Isub = self._ds.current[sub].data

        # generate fit curve
        a = self._ds.floating_potential_fit.sel(vf_fit_params='slope',
                                                vf_fit_values='value').data[()]
        b = self._ds.floating_potential_fit.sel(vf_fit_params='intercept',
                                                vf_fit_values='value').data[()]
        Isub_fit = a * Vsub + b

        #
        vf = self._ds.floating_potential.sel(vf='value').data[()]
        vf_err = self._ds.floating_potential.sel(vf='err').data[()]

        # calc plot limits
        xpad = 0.4 * np.abs(Vsub.max() - Vsub.min())
        ypad = 0.4 * np.abs(Isub.max() - Isub.min())
        xlim = [Vsub.min() - xpad, Vsub.max() + xpad]
        ylim = [Isub.min() - ypad, Isub.max() + ypad]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel('Probe Bias (V)')
        ax.set_ylabel('Current (A)')
        ax.set_title('Floating Potential V$_f$ Fit')

        # add plot elemnts
        ax.plot(ax.get_xlim(), [0.0, 0.0], 'r--', label='I = 0')
        ystop = ax.transData.inverted().transform(ax.transAxes.transform([0.0, 0.85]))[
            1]
        ax.plot([vf, vf],
                [ylim[0], ystop],
                'r',
                label='vf')
        ax.fill_between([vf - vf_err, vf + vf_err],
                        ylim[0],
                        ylim[1],
                        color='r', alpha=0.1)
        ax.plot(Vsub, Isub_fit, label='fit', color='orange')
        ax.scatter(self._ds['voltage'], self._ds['current'], label='data')
        ax.scatter(Vsub, Isub,
                   s=12 ** 2,
                   label='fit data',
                   facecolors='none',
                   edgecolors='orange')

        # labels
        txt = f"V$_f$ = {vf.data[()]:.2f} V"
        txt_loc = [vf, ylim[0]]
        txt_loc = ax.transAxes.inverted().transform(ax.transData.transform(txt_loc))
        txt_loc[1] = 0
        txt_loc += 0.03
        ax.text(txt_loc[0], txt_loc[1], txt, color='r', fontsize='large',
                transform=ax.transAxes)

        vals = [a, b]
        imax = np.abs(np.log10(vals)).argmax()
        power = np.log10(vals[imax])
        power = int(np.sign(power)) * int(np.ceil(np.abs(power)))
        txt = (f"I = ( {a * (10 ** -power):.1f} x 10$^{{{power}}}$ ) * V "
               f"+ ( {b * (10 ** -power):.1f} x 10$^{{{power}}}$ )")
        txt_loc = [0.0 + 0.03, 1.0 - 0.1]
        ax.text(txt_loc[0], txt_loc[1], txt, color='black', fontsize='large',
                transform=ax.transAxes)

        plt.show()

        return fig, ax
