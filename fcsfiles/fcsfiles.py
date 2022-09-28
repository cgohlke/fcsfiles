# fcsfiles.py

# Copyright (c) 2012-2022, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Read fluorescence correlation spectroscopy (FCS) data files.

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

"""

from __future__ import annotations

__version__ = '2022.9.28'

__all__ = ['ConfoCor3Fcs', 'ConfoCor3Raw', 'ConfoCor2Raw', 'fcs_bincount']

import os
import struct
from typing import Any, BinaryIO, Sequence

import numpy


class ConfoCor3Fcs(dict):
    """Carl Zeiss ConfoCor3 ASCII data file.

    No specification is available. The encoding is 'Windows-1252'.

    """

    HEADER = 'Carl Zeiss ConfoCor3 - measurement data file - version 3.0 ANSI'

    _filepath: str
    _filename: str

    def __init__(self, filename: os.PathLike | str, /) -> None:
        """Read file content and parse into dictionary."""
        dict.__init__(self)
        filename = os.path.abspath(os.fspath(filename))
        self._filepath, self._filename = os.path.split(filename)
        with open(filename, encoding='Windows-1252') as fh:
            header = fh.read(63)
            if header != ConfoCor3Fcs.HEADER:
                raise ValueError('not a ConfoCor3 measurement data file')
            fh.readline()
            current = self
            stack = []
            array: list[str] | None = []
            key = None
            for line in fh:
                line = line.lstrip()
                if line[0].isdigit():
                    if array is None:
                        array = [line]
                    else:
                        array.append(line)
                    continue
                if array:
                    shape = tuple(int(s) for s in current[key].split())
                    dtype = 'f8' if '.' in array[0] else 'i8'
                    ndarray = numpy.fromstring(
                        ''.join(array), dtype=dtype, sep=' '
                    ).reshape(*shape)
                    current[key] = ndarray
                    del ndarray
                    array = []
                elif array is None:
                    # work around https://github.com/numpy/numpy/issues/1714
                    shape = tuple(int(s) for s in current[key].split())
                    current[key] = numpy.zeros(shape)
                    array = []
                if line.startswith('BEGIN'):
                    stack.append(current)
                    _, key, valstr = line.split()
                    key = key.strip()
                    val = {'_value': int(valstr.strip())}
                    if key[-1].isdigit():
                        key = key[:-1]
                        if key[-1].isdigit():
                            key = key[:-1]
                        if key in current:
                            current[key].append(val)
                        else:
                            current[key] = [val]
                        current = current[key][-1]
                    else:
                        current[key] = val
                        current = current[key]
                elif line.startswith('END'):
                    current = stack.pop()
                elif '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    if key.endswith('ArraySize'):
                        continue
                    if key.endswith('Array'):
                        array = None
                    value = value.strip()
                    if value:
                        for type_t in (int, float, str):
                            try:
                                value = type_t(value)
                                break
                            except (ValueError, TypeError):
                                pass
                    if key[-1].isdigit():
                        key = key[:-1]
                        if key[-1].isdigit():
                            key = key[:-1]
                        if key in current:
                            current[key].append(value)
                        else:
                            current[key] = [value]
                    else:
                        current[key] = value

    def __str__(self) -> str:
        result = [ConfoCor3Fcs.HEADER]

        def append(key: str, value, indent: str = '', index: int = -1):
            # recursively append formatted keys and values to result
            idxstr = str(index + 1) if index >= 0 else ''
            if isinstance(value, dict):
                result.append(f'{indent}BEGIN {key}{idxstr} {value["_value"]}')
                for k, v in value.items():  # sorted(..., key=sortkey) ?
                    append(k, v, indent + '\t')
                result.append(f'{indent}END')
            elif isinstance(value, (list, tuple)):
                for i, val in enumerate(value):
                    append(key, val, indent, i)
            elif isinstance(value, numpy.ndarray):
                size = value.shape[0]
                if size != 1:
                    result.append(f'{indent}{key}Size = {size}')
                result.append(
                    '{}{} = {}'.format(
                        indent, key, ' '.join(str(i) for i in value.shape)
                    )
                )
                fmt = '{}' if value.dtype.kind in 'iu' else '{:.8f}'
                for i in range(size):
                    result.append(
                        '{}{}'.format(
                            indent, '\t '.join(fmt.format(v) for v in value[i])
                        )
                    )
            elif key != '_value':
                result.append(f'{indent}{key}{idxstr} = {value}')

        def sortkey(item: tuple[str, Any]) -> str:
            # sort dictionary items by key string and value type
            key, value = item
            key = key.lower()
            if isinstance(value, (list, tuple)):
                return '~' + key
            if isinstance(value, numpy.ndarray):
                return '~~' + key
            if isinstance(value, dict):
                return '~~~' + key
            return key

        for key, val in self.items():  # sorted(..., key=sortkey) ?
            append(key, val)
        return '\n'.join(result)

    def close(self) -> None:
        """Close open file."""

    def __enter__(self) -> ConfoCor3Fcs:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class ConfoCor3Raw:
    """Carl Zeiss ConfoCor3 raw data file.

    Based on "Confocor3 Raw Data Specification. 2 May 2007. by efg, Stowers
    Institute for Medical Research."

    Attributes
    ----------
    measurement_identifier : str
        Hex form of 4 32-bit little endian integers without leading 0.
        It is the same for all channels of the same repetition.
    measurement_position : int
        Zero based measurement position.
    kinetic_index : int
        Zero based kinetic index.
    repetition_number : int
        Zero based repetition_number.
    frequency : int
        Sampling frequency in Hz.

    """

    HEADER = b'Carl Zeiss ConfoCor3 - raw data file - version 3.000 - '

    file_identifier: bytes
    channel: int
    measurement_position: int
    kinetic_index: int
    repetition_number: int
    frequency: int
    measurement_identifier: str
    _fh: BinaryIO
    _filepath: str
    _filename: str

    def __init__(self, filename: os.PathLike | str) -> None:
        """Read file header."""
        filename = os.path.abspath(filename)
        self._filepath, self._filename = os.path.split(filename)
        self._fh = open(filename, 'rb')
        header = self._fh.read(64)
        if not header.startswith(ConfoCor3Raw.HEADER):
            self._fh.close()
            raise ValueError('not a ConfoCor3 raw data file')
        self.file_identifier = header
        self.channel = int(header.strip().rsplit(b' ', 1)[-1]) - 1
        (
            measureid,
            self.measurement_position,
            self.kinetic_index,
            self.repetition_number,
            self.frequency,
            _,
        ) = struct.unpack('<16sIIII32s', self._fh.read(64))
        measureid = struct.unpack('<IIII', measureid)
        measureid = ''.join(hex(int(i))[2:] for i in measureid)
        self.measurement_identifier = measureid.replace('L', '')

    def filename(self) -> str:
        """Return normalized file name from file content."""
        return '{}_R{}_P{}_K{}_Ch{}.raw'.format(
            self.measurement_identifier,
            self.repetition_number + 1,
            self.measurement_position + 1,
            self.kinetic_index + 1,
            self.channel + 1,
        )

    def asarray(
        self, count: int = -1, skip: int = 0, **kwargs
    ) -> numpy.ndarray | tuple[numpy.ndarray, ...]:
        """Read data from file, perform optional binning, and return as array.

        Parameters
        ----------
        count : int (optional)
            Number of data words to process. Default: -1 (all words).
        skip : int (optional)
            Number of data words to skip at beginning of stream. Default: 0.
        **kwargs
            Optional argument to the 'fcs_bincount' function, specifying
            size or number of bins.

        Returns
        -------
        times : ndarray
            The lower range of the bins in units of seconds, or, if 'binsize'
            is 0, the times of the events in detector clock units.
        bincounts : ndarray
            The number of events in each bin.

        """
        self._fh.seek(128 + skip * 4)
        times = numpy.fromfile(self._fh, dtype='<u4', count=count)
        times = times.astype('u8')
        times = numpy.cumsum(times, out=times)
        if kwargs:
            result = fcs_bincount((times,), self.frequency, **kwargs)
            return result[0], result[1][0]
        return times

    def __str__(self) -> str:
        return '\n '.join(
            (
                self.__class__.__name__,
                os.path.normpath(os.path.normcase(self.filename())),
                f'measurement identifier: {self.measurement_identifier}',
                f'sampling frequency: {self.frequency} Hz',
                f'repetition number {self.repetition_number}',
                f'measurement position: {self.measurement_position}',
                f'kinetic index: {self.kinetic_index}',
                f'channel number: {self.channel}',
            )
        )

    def close(self) -> None:
        """Close open file."""
        self._fh.close()

    def __enter__(self) -> ConfoCor3Raw:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class ConfoCor2Raw:
    """Carl Zeiss ConfoCor2 raw data file.

    The LSM 510 META - ConfoCor 2 system: an integrated imaging and
    spectroscopic platform for single-molecule detection.
    Weisshart K, Jungel V, and Briddon SJ.
    Curr Pharm Biotechnol. 2004 Apr;5(2):135-54.

    ConfoCor 2 - Description of the Raw Data Format. 03.05.01.
    www.zeiss.at/C12567BE00472A5C/EmbedTitelIntern/RawData/$File/RawData.pdf

    Carl Zeiss Advanced Imaging Microscopy claims a patent on the file format
    and data compression scheme.

    Attributes
    ----------
    frequency : int
        Sampling frequency in Hz.

    """

    HEADER = b'ConfoCor 2 - Raw data file 1.0'

    channels: int
    frequency: int
    file_identifier: bytes
    _fh: BinaryIO
    _filepath: str
    _filename: str

    def __init__(self, filename: os.PathLike | str, /) -> None:
        """Read file content and parse into a dictionary."""
        filename = os.path.abspath(filename)
        self._filepath, self._filename = os.path.split(filename)
        self._fh = open(filename, 'rb')
        header = self._fh.read(30)
        if not header.startswith(ConfoCor2Raw.HEADER):
            self._fh.close()
            raise ValueError('not a ConfoCor 2 raw data file')
        self.file_identifier = header
        self.channels = 2
        self.frequency = 20000000

    def asarray(
        self, count: int = -1, skip: int = 0, **kwargs
    ) -> tuple[numpy.ndarray, ...]:
        """Read data from file, perform optional binning, and return as arrays.

        Parameters
        ----------
        count : int (optional)
            Number of data words to process. Default: -1 (all words).
        skip : int (optional)
            Number of data words to skip at beginning of stream. Default: 0.
        **kwargs
            Optional argument to the 'fcs_bincount' function, specifying
            size or number of bins.

        Returns
        -------
        times : ndarray
            If binning, the lower range of the bins in units of seconds.
            Otherwise (default), the times of the events in detector clock
            units for each channel.
        Bincounts : ndarray
            If binning, the number of events in each bin for each channel.

        """
        self._fh.seek(30 + skip * 2)
        data = numpy.fromfile(self._fh, dtype='u1', count=count * 2)
        # accumulate clock time
        times = numpy.empty((len(data) // 2, 4), dtype='u8')
        times[:, 0] = data[::2]
        times[:, 1:] = 1
        times = times.flatten()
        numpy.cumsum(times, out=times)
        # extract events for each channel
        data = numpy.repeat(data[1::2], 8).reshape(-1, 8)
        data &= numpy.array([1, 2, 4, 8, 16, 32, 64, 128], dtype='u1')
        ch0 = numpy.take(times, numpy.where(data[:, 0::2].flatten() != 0)[0])
        ch1 = numpy.take(times, numpy.where(data[:, 1::2].flatten() != 0)[0])
        del data
        del times
        if kwargs:
            result = fcs_bincount((ch0, ch1), self.frequency, **kwargs)
            return result[0], result[1][0], result[1][1]
        return ch0, ch1

    def __str__(self) -> str:
        return '\n '.join(
            (
                self.__class__.__name__,
                os.path.normpath(os.path.normcase(self._filename)),
                f'sampling frequency: {self.frequency} Hz',
            )
        )

    def close(self) -> None:
        """Close open file."""
        self._fh.close()

    def __enter__(self) -> ConfoCor2Raw:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def fcs_bincount(
    data: Sequence[numpy.ndarray],
    frequency: float,
    binsize: int | None = None,
    bins: int | None = None,
    binspm: int | None = None,
) -> tuple[numpy.ndarray, tuple[numpy.ndarray, ...]]:
    """Count number of events in bins.

    Parameters
    ----------
    data : tuple of ndarray
        For each channel, the times of events in detector clock units.
    frequency : int
        Sampling frequency in Hz.
    binsize : int
        Time interval, as a multiple of the sampling frequency, in which
        to bin events.
    bins : int
        If 'binsize' is None, the number of bins.
    binspm : int
        If 'binsize' and 'bins' are None, the number of bins per minute (60s).

    Returns
    -------
    times : ndarray
        The lower range of the bins in units of seconds.
    bincounts : tuple of ndarray
        For each channel, the number of events in each bin.

    """
    assert type(data) in (tuple, list)
    if binsize is None:
        if bins is None:
            if binspm is None:
                raise ValueError('missing parameter binsize, bins, or binspm')
            binsize = int(60 * frequency) // binspm
        else:
            size = int(max((ch[-1] if ch.size else 0) for ch in data))
            binsize = size // bins + 1
    for ch in data:
        ch //= int(binsize)  # ch is ndarray
    size = int(max((ch[-1] if ch.size else 0) for ch in data) + 1)
    times = numpy.linspace(
        0, size * binsize / float(frequency), size, endpoint=False
    )
    # FIXME: work around https://github.com/numpy/numpy/issues/823
    if size < 2**31:
        # use 32 bit signed int if possible
        dataiter = (ch.astype('i4') for ch in data)
    else:
        # else use a signed view; fails on 32 bit
        dataiter = (ch.view('i8') for ch in data)
    bincounts = tuple(numpy.bincount(ch, minlength=size) for ch in dataiter)
    return times, bincounts


if __name__ == '__main__':
    import doctest

    numpy.set_printoptions(suppress=True, precision=5)
    doctest.testmod()
