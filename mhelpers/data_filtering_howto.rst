====================================
How to filter profile data with Mean
====================================

1. Select the algorithm and create the corresponding object (E.G.
   PLGaussianMean with 5 point interval and sigma of 1.0)::

    from mhelpers.pgmean import PLGaussianMean
    m = PLGaussianMean(5,1.0)

2. Pass the newly created object to the read() method of a Bio_Float
   object. E.G.::

    press, profile = float.read('TEMP', m)

3. You can also set this object as the default for the Bio_Float class::

    Bio_Float.default_mean = m

   If you do so every time you DON'T pass a mean object to read() the
   default will be used.

4. If you want to read the raw data without any filter call the
   read_raw() method::

    press, profile = float.read_raw('TEMP')

