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
import numpy, numpy.fft
import time

dll = ctypes.windll.LoadLibrary('dp7889.dll')

class ACQSETTING(ctypes.Structure):
    _fields_ = [('range',       ctypes.c_ulong),
                ('prena',       ctypes.c_long),
                ('cftfak',      ctypes.c_long),
                ('roimin',      ctypes.c_ulong),
                ('roimax',      ctypes.c_ulong),
                ('eventpreset', ctypes.c_double),
                ('timepreset',  ctypes.c_double),
                ('savedata',    ctypes.c_long),
                ('fmt',         ctypes.c_long),
                ('autoinc',     ctypes.c_long),
                ('cycles',      ctypes.c_long),
                ('sweepmode',   ctypes.c_long),
                ('syncout',     ctypes.c_long),
                ('bitshift',    ctypes.c_long),
                ('digval',      ctypes.c_long),
                ('digio',       ctypes.c_long),
                ('dac0',        ctypes.c_long),
                ('dac1',        ctypes.c_long),
                ('swpreset',    ctypes.c_double),
                ('nregions',    ctypes.c_long),
                ('caluse',      ctypes.c_long),
                ('fstchan',     ctypes.c_double),
                ('active',      ctypes.c_long),
                ('calpoints',   ctypes.c_long), ]

class ACQDATA(ctypes.Structure):
    _fields_ = [('s0', ctypes.POINTER(ctypes.c_ulong)),
                ('region', ctypes.POINTER(ctypes.c_ulong)),
                ('comment', ctypes.c_char_p),
                ('cnt', ctypes.POINTER(ctypes.c_double)),
                ('hs0', ctypes.c_int),
                ('hrg', ctypes.c_int),
                ('hcm', ctypes.c_int),
                ('hct', ctypes.c_int), ]

class ACQSTATUS(ctypes.Structure):
    _fields_ = [('started', ctypes.c_int),
                ('runtime', ctypes.c_double),
                ('totalsum', ctypes.c_double),
                ('roisum', ctypes.c_double),
                ('roirate', ctypes.c_double),
                ('ofls', ctypes.c_double),
                ('sweeps', ctypes.c_double),
                ('stevents', ctypes.c_double),
                ('maxval', ctypes.c_ulong), ]


class FastComTec(object):
    
    def __init__(self):
        self.BINWIDTH = 0.1

    def Configure(self, t, dt, triggers=1, SoftwareStart=False):
        self.SetSlots(triggers)
        self.SetSoftwareStart(SoftwareStart)
        N, DT = self.SetRange(t, dt)
        dll.RunCmd(0, 'ROIMAX=%i' % (N/triggers))
        self.Erase()
        return N, DT

    def SetRange(self, T, dT):
        self.SetBitshift( int(numpy.round(numpy.log2(dT / self.BINWIDTH))) )
        DT = self.BINWIDTH * 2**self.GetBitshift()
        self.SetLength( int(T / DT) )
        return self.GetLength(), DT

    def GetRange(self):
        return self.GetLength(), self.BINWIDTH * 2**self.GetBitshift()

    def SetLength(self, N):
        dll.RunCmd(0, 'RANGE=%i' % N)
        return self.GetLength()

    def GetLength(self):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        return int(setting.range)

    def SetBitshift(self, BITSHIFT):
        #~ setting = ACQSETTING()
        #~ dll.GetSettingData(ctypes.byref(setting), 0)
        #~ setting.bitshift = BITSHIFT
        #~ dll.StoreSettingData(ctypes.byref(setting), 0)
        #~ dll.NewSetting(0)
        dll.RunCmd(0,'BITSHIFT=%i'%BITSHIFT)
        return self.GetBitshift()

    def GetBitshift(self):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        return int(setting.bitshift)

    def SetSlots(self, M):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        if M == 1:
            setting.sweepmode = int('10100000', 2)
            #~ dll.RunCmd(0, 'SWEEPMODE=' + hex(int('10100000', 2)))
        elif M > 1:
            setting.sweepmode = int('10100100',2)
            #~ dll.RunCmd(0, 'SWEEPMODE=' + hex(int('10100100',2)))
            #~ dll.RunCmd(0, 'SWPRESET=1')
            setting.prena = setting.prena |  16
            setting.prena = setting.prena &~ 4
            setting.swpreset = 1
        setting.cycles = M
        #~ dll.RunCmd(0, 'CYCLES=%i' % M)
        dll.StoreSettingData(ctypes.byref(setting), 0)
        dll.NewSetting(0)

    def SetCycles(self, cycles):
        cycles = numpy.min((int('ffffffff', 16)/2, cycles))
        dll.RunCmd(0, 'SEQUENCES=%i' % cycles)

    def SetSweeps(self, sweeps):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        if sweeps == numpy.Inf:
            setting.prena = setting.prena &~ 16
            setting.prena = setting.prena &~ 4
            #~ dll.RunCmd(0, 'PRENA=' + hex((setting.prena &~ 16)))
            #~ dll.RunCmd(0, 'PRENA=' + hex((setting.prena &~ 4)))
        else:
            # if start event generation is set, set 'start presets'
            # else set 'sweep presets'
            if ( setting.sweepmode & int('10000000',2) ) >> 7:
                setting.prena = setting.prena | 16
                setting.prena = setting.prena &~ 4
                #~ dll.RunCmd(0, 'PRENA=' + hex((setting.prena | 16)))
            else:
                setting.prena = setting.prena &~ 16
                setting.prena = setting.prena | 4
                #~ dll.RunCmd(0, 'PRENA=' + hex((setting.prena | 4)))
            setting.swpreset = float(sweeps)
            #~ dll.RunCmd(0, 'SWPRESET=%i' % sweeps)
        dll.StoreSettingData(ctypes.byref(setting), 0)
        dll.NewSetting(0)

    def SetTime(self, time):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        if time == numpy.Inf:
            #~ setting.prena = setting.prena &~1
            dll.RunCmd(0, 'PRENA=' + hex((setting.prena &~ 1)))
        else:
            #~ setting.prena = setting.prena|1
            #~ setting.timepreset = time
            dll.RunCmd(0, 'PRENA=' + hex((setting.prena | 1)))
            dll.RunCmd(0, 'RTPRESET=%f' % time)
        #~ dll.StoreSettingData(ctypes.byref(setting), 0)
        #~ dll.NewSetting(0)

    def SetSoftwareStart(self,b):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        if b:
            setting.sweepmode = setting.sweepmode |  int('10000',2)
            setting.sweepmode = setting.sweepmode &~ int('10000000',2)
        else:
            setting.sweepmode = setting.sweepmode &~ int('10000',2)
            setting.sweepmode = setting.sweepmode |  int('10000000',2)
        dll.StoreSettingData(ctypes.byref(setting), 0)
        dll.NewSetting(0)

    def SetDelay(self, t):
        #~ setting = ACQSETTING()
        #~ dll.GetSettingData(ctypes.byref(setting), 0)
        #~ setting.fstchan = t/6.4
        #~ dll.StoreSettingData(ctypes.byref(setting), 0)
        #~ dll.NewSetting(0)
        dll.RunCmd(0, 'DELAY=%f' % t)
        return self.GetDelay()

    def GetDelay(self):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        return setting.fstchan * 6.4

    def Start(self):
        dll.Start(0)
        status = ACQSTATUS()
        status.started = 0
        while not status.started:
            time.sleep(0.1)
            dll.GetStatusData(ctypes.byref(status), 0)

    def Halt(self):
        dll.Halt(0)

    def Erase(self):
        dll.Erase(0)

    def Continue(self):
        return dll.Continue(0)

    def GetData(self):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        N = setting.range
        M = setting.cycles
        data = numpy.empty((N,), dtype=numpy.uint32 )       
        dll.LVGetDat(data.ctypes.data, 0)
        if M == 1:
            return data
        else:
            return data.reshape((M,N/M))

    def GetFFT(self):
        StartTime = time.time()
        self.Start()
        while self.Running():
            time.sleep(0.1)
        data = numpy.zeros((2*self.GetLength(),), dtype=numpy.uint32 )
        dll.LVGetDat(data.ctypes.data, 0)
        r = data.astype(numpy.float)
        F = numpy.fft.rfft(r)
        return F

    def AccumulateFFT(self, N):
        Y = numpy.zeros((self.GetLength()+1,))
        i = 0
        while i < N:
            self.Start()
            while self.Running():
                time.sleep(0.1)
            Y += abs(self.GetFFT())
            i += 1
        return Y

    def GetState(self):
        status = ACQSTATUS()
        dll.GetStatusData(ctypes.byref(status), 0)
        return status.runtime, status.sweeps

    def Running(self):
        s = self.GetStatus()
        return s.started

    def SetLevel(self, start, stop):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        def FloatToWord(r):
            return int((r+2.048)/4.096*int('ffff',16))
        setting.dac0 = ( setting.dac0 & int('ffff0000',16) ) | FloatToWord(start)
        setting.dac1 = ( setting.dac1 & int('ffff0000',16) ) | FloatToWord(stop)
        dll.StoreSettingData(ctypes.byref(setting), 0)
        dll.NewSetting(0)
        return self.GetLevel()

    def GetLevel(self):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        def WordToFloat(word):
            return (word & int('ffff',16)) * 4.096 / int('ffff',16) - 2.048
        return WordToFloat(setting.dac0), WordToFloat(setting.dac1)

    def ReadSetting(self):
        setting = ACQSETTING()
        dll.GetSettingData(ctypes.byref(setting), 0)
        return setting

    def WriteSetting(self, setting):
        dll.StoreSettingData(ctypes.byref(setting), 0)
        dll.NewSetting(0)

    def GetStatus(self):
        status = ACQSTATUS()
        dll.GetStatusData(ctypes.byref(status), 0)
        return status
