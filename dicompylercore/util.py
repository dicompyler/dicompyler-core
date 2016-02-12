#!/usr/bin/env python
# -*- coding: utf-8 -*-
# util.py
"""Several utility functions that don't really belong anywhere."""
# Copyright (c) 2009-2016 Aditya Panchal
# This file is part of dicompyler, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at http://code.google.com/p/dicompyler/

import numpy as np
import collections
import sys


def piecewise(x, condlist, funclist, *args, **kw):
    """
    Evaluate a piecewise-defined function.
    (Borrowed from numpy 1.11 for compatibility for numpy 1.9-1.10)

    Given a set of conditions and corresponding functions, evaluate each
    function on the input data wherever its condition is true.

    Parameters
    ----------
    x : ndarray
        The input domain.
    condlist : list of bool arrays
        Each boolean array corresponds to a function in `funclist`.  Wherever
        `condlist[i]` is True, `funclist[i](x)` is used as the output value.

        Each boolean array in `condlist` selects a piece of `x`,
        and should therefore be of the same shape as `x`.

        The length of `condlist` must correspond to that of `funclist`.
        If one extra function is given, i.e. if
        ``len(funclist) - len(condlist) == 1``, then that extra function
        is the default value, used wherever all conditions are false.
    funclist : list of callables, f(x,*args,**kw), or scalars
        Each function is evaluated over `x` wherever its corresponding
        condition is True.  It should take an array as input and give an array
        or a scalar value as output.  If, instead of a callable,
        a scalar is provided then a constant function (``lambda x: scalar``) is
        assumed.
    args : tuple, optional
        Any further arguments given to `piecewise` are passed to the functions
        upon execution, i.e., if called ``piecewise(..., ..., 1, 'a')``, then
        each function is called as ``f(x, 1, 'a')``.
    kw : dict, optional
        Keyword arguments used in calling `piecewise` are passed to the
        functions upon execution, i.e., if called
        ``piecewise(..., ..., lambda=1)``, then each function is called as
        ``f(x, lambda=1)``.

    Returns
    -------
    out : ndarray
        The output is the same shape and type as x and is found by
        calling the functions in `funclist` on the appropriate portions of `x`,
        as defined by the boolean arrays in `condlist`.  Portions not covered
        by any condition have a default value of 0.


    See Also
    --------
    choose, select, where

    Notes
    -----
    This is similar to choose or select, except that functions are
    evaluated on elements of `x` that satisfy the corresponding condition from
    `condlist`.

    The result is::

            |--
            |funclist[0](x[condlist[0]])
      out = |funclist[1](x[condlist[1]])
            |...
            |funclist[n2](x[condlist[n2]])
            |--

    Examples
    --------
    Define the sigma function, which is -1 for ``x < 0`` and +1 for ``x >= 0``.

    >>> x = np.linspace(-2.5, 2.5, 6)
    >>> np.piecewise(x, [x < 0, x >= 0], [-1, 1])
    array([-1., -1., -1.,  1.,  1.,  1.])

    Define the absolute value, which is ``-x`` for ``x <0`` and ``x`` for
    ``x >= 0``.

    >>> np.piecewise(x, [x < 0, x >= 0], [lambda x: -x, lambda x: x])
    array([ 2.5,  1.5,  0.5,  0.5,  1.5,  2.5])

    """

    # Use the built-in numpy piecewise function if not on numpy 1.9-1.10
    version = np.version.version.split('.')
    if (version[0] == '1') and (version[1] not in ['9', '10']):
        return np.piecewise(x, condlist, funclist, *args, **kw)

    x = np.asanyarray(x)
    n2 = len(funclist)
    if (np.isscalar(condlist) or not (isinstance(condlist[0], list) or
                                      isinstance(condlist[0], np.ndarray))):
        condlist = [condlist]
    condlist = np.array(condlist, dtype=bool)
    n = len(condlist)
    # This is a hack to work around problems with NumPy's
    #  handling of 0-d arrays and boolean indexing with
    #  numpy.bool_ scalars
    zerod = False
    if x.ndim == 0:
        x = x[None]
        zerod = True
        if condlist.shape[-1] != 1:
            condlist = condlist.T
    if n == n2 - 1:  # compute the "otherwise" condition.
        totlist = np.logical_or.reduce(condlist, axis=0)
        # Only able to stack vertically if the array is 1d or less
        if x.ndim <= 1:
            condlist = np.vstack([condlist, ~totlist])
        else:
            condlist = [np.asarray(c, dtype=bool) for c in condlist]
            totlist = condlist[0]
            for k in range(1, n):
                totlist |= condlist[k]
            condlist.append(~totlist)
        n += 1

    y = np.zeros(x.shape, x.dtype)
    for k in range(n):
        item = funclist[k]
        if not isinstance(item, collections.Callable):
            y[condlist[k]] = item
        else:
            vals = x[condlist[k]]
            if vals.size > 0:
                y[condlist[k]] = item(vals, *args, **kw)
    if zerod:
        y = y.squeeze()
    return y


def platform():
    if sys.platform.startswith('win'):
        return 'windows'
    elif sys.platform.startswith('darwin'):
        return 'mac'
    return 'linux'
