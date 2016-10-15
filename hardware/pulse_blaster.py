"""
    This file is part of pi3diamond, a toolkit for
    confocal scanning, anti-bunching, FLIM, pulsed ODMR / NMR,
    and more sophisticated quantum physics experiments,
    typically performed with NV centers in diamond,
    written in python using the enthought traits packages.

    pi3diamond is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pi3diamond is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diamond. If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2009-2016 Helmut Fedder <helmut@fedder.net>
"""

import ctypes
import numpy
import traceback

map = {'laser':2, 'mw':1, 'microwave':1, 'trigger':3, 'SequenceTrigger':0, 'awgTrigger':0}

CLOCK = 400.0
DT = 2.5

dll = ctypes.cdll.LoadLibrary('spinapi.dll')

PULSE_PROGRAM  = 0
CONTINUE       = 0
STOP           = 1
LOOP           = 2
END_LOOP       = 3
LONG_DELAY     = 7
BRANCH         = 6
#ON             = 6<<21 # this doesn't work even though it is according to documentation
ON             = 0xE00000

def chk(err):
    """a simple error checking routine"""
    if err < 0:
        dll.pb_get_error.restype = ctypes.c_char_p
        err_str = dll.pb_get_error()
        #raise RuntimeError('PulseBlaster error: %s' % err_str)
        print err_str
        pi3d.get_logger().error('PulseBlaster error: ' + err_str + ''.join(traceback.format_stack()))    

def High(channels):
    """Set specified channels to high, all others to low."""
    pi3d.get_logger().debug(str(channels))
    dll.pb_close()
    chk(dll.pb_init())
    chk(dll.pb_set_clock(ctypes.c_double(CLOCK)))
    chk(dll.pb_start_programming(PULSE_PROGRAM))
    chk(dll.pb_inst_pbonly(ON|flags(channels), STOP, None, ctypes.c_double(100)))
    chk(dll.stop_programming())
    chk(dll.start_pb())
    chk(dll.pb_close())

def Sequence(sequence, loop=True):
    """Run sequence of instructions"""
    pi3d.get_logger().debug(str(sequence))
    #we pb_close() without chk(): this might create an error if board was already closed, but resets it if it was still open
    dll.pb_close()
    chk(dll.pb_init())
    chk(dll.pb_set_clock(ctypes.c_double(CLOCK)))
    chk(dll.pb_start_programming(PULSE_PROGRAM))
    start = write(*sequence[0])
    for step in sequence[1:]:
        label = write(*step)
        if label > 2**12 - 2:
            print 'WARNING in PulseBlaster: command %i exceeds maximum number of commands.' % label
    channels, dt = sequence[-1]
#    N = int(numpy.round(dt/DT))
#    if N > 256:
#        raise RuntimeError, 'Length in STOP / BRANCH exceeds maximum value.'
    if loop == False:
        label = chk(dll.pb_inst_pbonly(ON|flags(channels), STOP, None, ctypes.c_double( 12.5 ) ))
    else:
        label = chk(dll.pb_inst_pbonly(ON|flags([]), BRANCH, start, ctypes.c_double( 12.5 ) ))
    if label > 2**12 - 2 :
        print 'WARNING in PulseBlaster: command %i exceeds maximum number of commands.' % label
    chk(dll.stop_programming())
    chk(dll.start_pb())
    chk(dll.pb_close())


def write(channels, dt):
    channel_bits = flags(channels)
    N = int(numpy.round( dt / DT ))
#    if N == 0 or N == 1:
#        label = chk(dll.pb_inst_direct( ctypes.byref(ctypes.c_int(1<<21|channel_bits)), CONTINUE, None, 2 ))
#    elif N == 2:
#        label = chk(dll.pb_inst_pbonly( 2<<21|channel_bits, CONTINUE, None, ctypes.c_double( 2*DT ) ))
#    elif N == 3:
#        label = chk(dll.pb_inst_pbonly( 3<<21|channel_bits, CONTINUE, None, ctypes.c_double( 2*DT ) ))
#    elif N == 4:
#        label = chk(dll.pb_inst_pbonly( 4<<21|channel_bits, CONTINUE, None, ctypes.c_double( 2*DT ) ))
#    elif N == 5:
#        label = chk(dll.pb_inst_direct( ctypes.byref(ctypes.c_int(5<<21|channel_bits)), CONTINUE, None, 3 ))
    if N <= 256:
        label = chk(dll.pb_inst_pbonly( ON|channel_bits, CONTINUE, None, ctypes.c_double( N*DT ) ))
    else:
        # successively try factorization, reducing N, and putting the subtracted amount into an additional short command if necessary
        B = N
        i = 4
        while True:
            M, K = factor(N)
            #print M, K, i
            if M > 4:
                if K == 1:
                    label = chk(dll.pb_inst_pbonly( ON|channel_bits, CONTINUE, None, ctypes.c_double( M*DT ) ))
                elif K < 1048577:
                    label = chk(dll.pb_inst_pbonly( ON|channel_bits, LONG_DELAY, K, ctypes.c_double( M*DT ) ))
                else:
                    raise RuntimeError, 'Loop count in LONG_DELAY exceedes maximum value.'
                if i > 4:
                    chk(dll.pb_inst_pbonly( ON|channel_bits, CONTINUE, None, ctypes.c_double( i*DT )  ))
                break
            i += 1
            N = B - i
    return label

def flags(channels):
    bits = 0
    for channel in channels:
        bits = bits | 1<<map[channel]
    return bits

#def flags(channels):
#    bits = numpy.zeros((24,), dtype=numpy.bool)
#    for channel in channels:
#        bits[map[channel]] = 1
#    s = ''
#    for bit in bits:
#        s = '%i'%bit + s
#    return int(s,2)

def factor(x):
    i = 256
    while i > 4:
        if x % i == 0:
            return i, x/i
        i -= 1
    return 1, x

def Light():
    High(['laser'])

def Night():
    High([])

def Open():
    High(['laser','mw'])

def test():
    Sequence([ (['laser'], 642.5),
               ([       ], 12.5),
               #(['laser'], 642.5),
               #([       ], 12.5),
               #(['laser'], 1282.5),
               #([       ], 12.5),
               #(['laser'], 3000),
               #([       ], 12.5),
            ], loop=True)


if __name__ == '__main__':
    pass