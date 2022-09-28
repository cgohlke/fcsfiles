Read fluorescence correlation spectroscopy (FCS) data files
===========================================================

Fcsfiles is a Python library to read Carl Zeiss(r) ConfoCor(r) RAW and ASCII
measurement data files.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD 3-Clause
:Version: 2022.9.28

Requirements
------------

This release has been tested with the following requirements and dependencies
(other versions may work):

- `CPython 3.8.10, 3.9.13, 3.10.7, 3.11.0rc2 <https://www.python.org>`_
- `Numpy 1.22.4 <https://pypi.org/project/numpy/>`_

Revisions
---------

2022.9.28

- Update metadata.

2022.2.2

- Add type hints.
- Use float64 or int64 for ConfoCor3Fcs arrays.
- Drop support for Python 3.7 and numpy < 1.19 (NEP29).

2021.6.6

- Remove support for Python 3.6 (NEP 29).

2020.9.18

- Relax ConfoCor3Raw header requirement.
- Support os.PathLike file names.

2020.1.1

- Remove support for Python 2.7 and 3.5.

Notes
-----

"Carl Zeiss" and "ConfoCor" are registered trademarks of Carl Zeiss, Inc.

The use of this implementation may be subject to patent or license
restrictions.

The API is not stable yet and is expected to change between revisions.

This module does *not* read flow cytometry standard FCS files.

Examples
--------

Read the CountRateArray from a ConfoCor3 ASCII file as a numpy array:

>>> fcs = ConfoCor3Fcs('ConfoCor3.fcs')
>>> fcs['FcsData']['FcsEntry'][0]['FcsDataSet']['CountRateArray'].shape
(60000, 2)

Read data and metadata from a ConfoCor3 RAW file:

>>> fcs = ConfoCor3Raw('ConfoCor3.raw')
>>> fcs.filename()
'f5ee4f36488fca2f89cb6b8626111006_R1_P1_K1_Ch1.raw'
>>> fcs.frequency
20000000
>>> times = fcs.asarray()
>>> times[10858]
1199925494
>>> times, bincounts = fcs.asarray(bins=1000)
>>> times.shape
(1000,)
>>> bincounts[618]
23
>>> fcs.close()

Read data and metadata from a ConfoCor2 RAW file:

>>> fcs = ConfoCor2Raw('ConfoCor2.raw')
>>> fcs.frequency
20000000
>>> ch0, ch1 = fcs.asarray()
>>> ch1[4812432]
999999833
>>> times, ch0, ch1 = fcs.asarray(bins=1000)
>>> times.shape
(1000,)
>>> ch1[428]
10095
>>> fcs.close()
