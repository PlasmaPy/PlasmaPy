import abc


class DistributionFunction(metaclass=abc.ABCMeta):
    """
    Class to represent either a single or multiple distribution functions.
    """
    @property
    def n_spatial_dims(self):
        return self._n_spatial_dims

    @n_spatial_dims.setter
    def n_spatial_dims(self, n_spatial_dims):
        if not 0 <= n_spatial_dims <= 3:
            raise ValueError('Number of spatial dimensions must be bewteen '
                             f'0 and 3 (got {n_spatial_dims}).')
        self._n_spatial_dims = n_spatial_dims

    @property
    def n_velocity_dims(self):
        return self._n_velocity_dims

    @n_velocity_dims.setter
    def n_velocity_dims(self, n_velocity_dims):
        if not 1 <= n_velocity_dims <= 3:
            raise ValueError('Number of velocity dimensions must be bewteen '
                             f'1 and 3 (got {n_velocity_dims}).')
        self._n_velocity_dims = n_velocity_dims

    @property
    def n_time_dims(self):
        return self._n_time_dims

    @n_time_dims.setter
    def n_time_dims(self, n_time_dims):
        if not 0 <= n_time_dims <= 1:
            raise ValueError('Number of time dimensions must be bewteen '
                             f'0 and 1 (got {n_time_dims}).')
        self._n_time_dims = n_time_dims


class DiscreteDistributionFunction(DistributionFunction):
    """
    Class to represent a discretely sampled distribution function.

    Parameters
    ----------
    data : array
        Distribution function samples.
    indexes : list of 1D array
        The index values for each dimension in *data*.
        The length of ``indexes[i]`` must equal ``data.shape[i]``.
    index_types : list of {'time', 'space', 'velocity'}
        The type of each dimension.
    """
    def __init__(self, data, indexes, index_types):
        # Take the shape of data as 'ground truth' to do input checks against
        index_lens = [len(index) for index in indexes]
        if tuple(index_lens) != tuple(data.shape):
            raise ValueError(f'Index lengths ({index_lens}) does not match '
                             f'the data shape ({data.shape})')

        if len(index_types) != len(indexes):
            raise ValueError(f'Number of index units ({len(index_types)}) '
                             'does not match number of indexes '
                             f'({len(indexes)}')

        self.n_spatial_dims = index_types.count('space')
        self.n_velocity_dims = index_types.count('velocity')
        self.n_time_dims = index_types.count('time')

        self._vel_indexes = []
        vel_idxs = index_types.index('velocity')
        for i in range(vel_idxs):
            self._vel_indexes.append(indexes[i])

        if self.n_time_dims:
            self._time_index = indexes[index_types.index('time')]

        if self.n_spatial_dims:
            self._spatial_indexes = []
            spatial_idxs = index_types.index('space')
            for i in range(vel_idxs):
                self._spatial_indexes.append(indexes[i])

        self.data = data

        @property
        def time_index(self):
            return self._time_index

        @property
        def spatial_indexes(self):
            return self._spatial_indexes

        @property
        def velocity_indexes(self):
            return self._vel_indexes


class AnalyticDistributionFunction(DistributionFunction):
    """
    Class to represent an analytic distribution function.

    Parameters
    ----------
    f : callable
        Function representing the analytic distribution function. The signature
        of f must be ...
    """
    def __init__(self, f):
        # TODO: Inspect and check signature of f to extract dimensions
        pass
