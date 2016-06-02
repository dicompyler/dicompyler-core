#!/usr/bin/env python
# -*- coding: utf-8 -*-
# dvh.py
"""Class that stores dose volume histogram (DVH) data."""
# Copyright (c) 2016 Aditya Panchal
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

import numpy as np
import logging
logger = logging.getLogger('dicompyler.dvh')

# Set default absolute dose and volume  units
abs_dose_units = 'gy'
abs_volume_units = 'cm3'
relative_units = '%'


class DVH:
    """Class that stores dose volume histogram (DVH) data in cGy."""

    def __init__(self, counts, bins,
                 dvh_type='cumulative',
                 dose_units=abs_dose_units,
                 volume_units=abs_volume_units):
        """Initialization for a DVH from existing histogram counts and bins.

        Parameters
        ----------
        counts : iterable or numpy array
            An iterable of volume or percent count data
        bins : iterable or numpy array
            An iterable of dose bins
        dvh_type : str, optional
            Choice of 'cumulative' or 'differential' type of DVH
        dose_units : str, optional
            Absolute dose units, i.e. 'gy' or relative units '%'
        volume_units : str, optional
            Absolute volume units, i.e. 'cm3' or relative units '%'
        """
        self.counts = np.array(counts)
        self.bins = np.array(bins) if bins[0] == 0 else np.append([0], bins)
        self.dvh_type = dvh_type.lower()
        self.dose_units = dose_units.lower()
        self.volume_units = volume_units.lower()

    @classmethod
    def from_dicom_dvh(cls, dataset, sequence_num):
        """Initialization for a DVH from a pydicom RT Dose DVH sequence."""
        dvh = dataset.DVHSequence[sequence_num]
        data = np.array(dvh.DVHData)
        return cls(counts=data[1::2] * dvh.DVHDoseScaling,
                   bins=data[0::2].cumsum(),
                   dvh_type=dvh.DVHType,
                   dose_units=dvh.DoseUnits,
                   volume_units=dvh.DVHVolumeUnits)

    @classmethod
    def from_data(cls, data, binsize=1):
        """Initialization for a DVH from raw data.

        Parameters
        ----------
        data : iterable or numpy array
            An iterable of dose data that is used to create the histogram
        binsize : int, optional
            Bin width size (in cGy used to create the histogram)
        """
        data = np.array(data)
        bins = np.arange(0, data.max() + 1, binsize)
        if bins.size == 1:
            bins = np.array([0, data.max()])
        if data.max() not in bins:
            bins = np.append(bins, data.max())
        counts, bins = np.histogram(data, bins)

        return cls(counts, bins)

    def __repr__(self):
        """String representation of the class."""
        return 'DVH(%s, %r bins: [%r:%r] %s, volume: %r %s)' % \
            (self.dvh_type, self.counts.size, self.bins.min(),
                self.bins.max(), self.dose_units.capitalize(),
                self.volume, self.volume_units.lower())

    def __eq__(self, other):
        """Comparison method between two DVH objects.

        Parameters
        ----------
        other : DVH
            Other DVH object to compare with

        Returns
        -------
        Bool
            True or False if the DVHs have equal attribs and via numpy.allclose
        """
        attribs = ('dvh_type', 'dose_units', 'volume_units')
        attribs_eq = {k: self.__dict__[k] for k in attribs} == \
            {k: other.__dict__[k] for k in attribs}
        return attribs_eq and \
            np.allclose(self.counts, other.counts) and \
            np.allclose(self.bins, other.bins)

# ============================= DVH properties ============================= #

    @property
    def bincenters(self):
        """Return a numpy array containing the bin centers."""
        return 0.5 * (self.bins[1:] + self.bins[:-1])

    @property
    def differential(self):
        """Return a differential DVH from a cumulative DVH."""
        dvh_type = 'differential'
        if self.dvh_type == dvh_type:
            return self
        else:
            return DVH(counts=np.append(abs(np.diff(self.counts) * -1), [0]),
                       bins=self.bins,
                       dvh_type=dvh_type,
                       dose_units=self.dose_units,
                       volume_units=self.volume_units)

    @property
    def cumulative(self):
        """Return a cumulative DVH from a differential DVH."""
        dvh_type = 'cumulative'
        if self.dvh_type == dvh_type:
            return self
        else:
            return DVH(counts=self.counts[::-1].cumsum()[::-1],
                       bins=self.bins,
                       dvh_type=dvh_type,
                       dose_units=self.dose_units,
                       volume_units=self.volume_units)

    def absolute_dose(self, rx_dose, dose_units=abs_dose_units):
        """Return an absolute dose DVH."""
        if self.dose_units == dose_units:
            return self
        else:
            return DVH(counts=self.counts,
                       bins=self.bins * rx_dose / 100,
                       dvh_type=self.dvh_type,
                       dose_units=dose_units,
                       volume_units=self.volume_units)

    def relative_dose(self, rx_dose):
        """Return a relative dose DVH based on a prescription dose.

        Parameters
        ----------
        rx_dose : number
            Prescription dose value used to normalize dose bins
        """
        dose_units = relative_units
        if self.dose_units == dose_units:
            return self
        else:
            return DVH(counts=self.counts,
                       bins=100 * self.bins / rx_dose,
                       dvh_type=self.dvh_type,
                       dose_units=dose_units,
                       volume_units=self.volume_units)

    def absolute_volume(self, volume, volume_units=abs_volume_units):
        """Return an absolute volume DVH.

        Parameters
        ----------
        volume : number
            Absolute volume of the structure
        volume_units : str, optional
            Units for the absolute volume
        """
        if self.volume_units == volume_units:
            return self
        else:
            return DVH(counts=volume * self.counts / 100,
                       bins=self.bins,
                       dvh_type=self.dvh_type,
                       dose_units=self.dose_units,
                       volume_units=volume_units)

    @property
    def relative_volume(self):
        """Return a relative volume DVH."""
        volume_units = relative_units
        if self.volume_units == relative_units:
            return self
        # Convert back to cumulative before returning a relative volume
        elif self.dvh_type == 'differential':
            return self.cumulative.relative_volume.differential
        else:
            return DVH(counts=100 * self.counts / self.counts.max(),
                       bins=self.bins,
                       dvh_type=self.dvh_type,
                       dose_units=self.dose_units,
                       volume_units=volume_units)

    @property
    def max(self):
        """Return the maximum dose."""
        diff = self.differential
        # Find the the maximum non-zero dose bin
        return diff.bins[1:][diff.counts > 0][-1]

    @property
    def min(self):
        """Return the minimum dose."""
        diff = self.differential
        # Find the the minimum non-zero dose bin
        return diff.bins[1:][diff.counts > 0][0]

    @property
    def mean(self):
        """Return the mean dose."""
        diff = self.differential
        # Find the area under the differential histogram
        return (diff.bincenters * diff.counts).sum() / diff.counts.sum()

    @property
    def volume(self):
        """Return the volume of the structure."""
        return self.differential.counts.sum()

    def plot(self):
        """Plot the DVH using Matplotlib if present."""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('Matplotlib could not be loaded. Install and try again.')
        else:
            plt.plot(self.bincenters, self.counts)
            # plt.axis([0, self.bins[-1], 0, self.counts[0]])
            plt.xlabel('Dose [%s]' % self.dose_units.capitalize())
            plt.ylabel('Volume [%s]' % self.volume_units.lower())
        return self
