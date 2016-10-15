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

from nidaq import *

import TimeTagger

import numpy

import time

class PulsedAO(object):
    """Provides pulsed analog output. Provides two output forms: 1. set ao channels permanently to a given value.
    2. Output N values along with a square wave with given seconds per point."""

    def __init__(self, channels, square_wave_device, seconds_per_point, duty_cycle=0.5, range=(-10,10), transfer_timeout=1.0):

        # we start with unknown length
        self._length = None
        self._square_wave_device = square_wave_device
        self._seconds_per_point = seconds_per_point 
        self._duty_cycle = duty_cycle

        self._channels = channels
        self._transfer_timeout = transfer_timeout
        
        self._n_written = ctypes.c_long()

        # create nidaq tasks
        self._ao_task = ctypes.c_ulong()
        self._co_task = ctypes.c_ulong()
        CHK( dll.DAQmxCreateTask('', ctypes.byref(self._ao_task)) )
        CHK( dll.DAQmxCreateTask('', ctypes.byref(self._co_task))  )

        CHK( dll.DAQmxCreateAOVoltageChan(self._ao_task, self._channels, '', ctypes.c_double(range[0]), ctypes.c_double(range[1]), DAQmx_Val_Volts,'')  )

        CHK(  dll.DAQmxCreateCOPulseChanFreq( self._co_task,
                                              self._square_wave_device, '',
                                              DAQmx_Val_Hz, DAQmx_Val_Low, ctypes.c_double(0),
                                              ctypes.c_double(1./seconds_per_point),
                                              ctypes.c_double(duty_cycle) )  )

    def output(self, data):
        """If len(data)==1, output one point. Else output sequence of points.
        The shape of data is (n_channels,n_points)."""
        if data.ndim == 1:
            N = 1
        else:
            N = data.shape[1]
        self.setLength(N)
        if N == 1:
            CHK( dll.DAQmxWriteAnalogF64(self._ao_task, ctypes.c_long(1), 1, ctypes.c_double(self._transfer_timeout), DAQmx_Val_GroupByChannel, data.ctypes.data, ctypes.byref(self._n_written), None) )
        else:
            CHK( dll.DAQmxWriteAnalogF64(self._ao_task, ctypes.c_long(N), 0, ctypes.c_double(self._transfer_timeout), DAQmx_Val_GroupByChannel, data.ctypes.data, ctypes.byref(self._n_written), None) )
            CHK( dll.DAQmxStartTask(self._ao_task)  )
            CHK( dll.DAQmxStartTask(self._co_task)  )
            CHK( dll.DAQmxWaitUntilTaskDone(self._ao_task, ctypes.c_double(4*(N+1)*self._seconds_per_point)) )
            time.sleep(4.*self._seconds_per_point)
            CHK( dll.DAQmxStopTask(self._co_task)  )
            CHK( dll.DAQmxStopTask(self._ao_task)  )
        return self._n_written            

    def setLength(self, N):
        if self._length == N:
            return
        if N == 1:
            CHK( dll.DAQmxSetSampTimingType( self._ao_task, DAQmx_Val_OnDemand)  )
        else:
            CHK( dll.DAQmxSetSampTimingType( self._ao_task, DAQmx_Val_SampClk)  )
            CHK( dll.DAQmxCfgSampClkTiming( self._ao_task,
                                            self._square_wave_device+'InternalOutput',
                                            ctypes.c_double(1./self._seconds_per_point),
                                            DAQmx_Val_Rising,
                                            DAQmx_Val_FiniteSamps,
                                            ctypes.c_ulonglong(N)) )
            CHK(  dll.DAQmxCfgImplicitTiming( self._co_task, DAQmx_Val_FiniteSamps, ctypes.c_ulonglong(N+1))  )
        self._length = N

    def getLength(self):
        return self._length

    def setTiming(self, seconds_per_point, duty_cycle=0.5):
        if seconds_per_point == self._seconds_per_point and duty_cycle == self._duty_cycle:
            return
        CHK( dll.DAQmxSetCOPulseFreq( self._co_task, self._square_wave_device, ctypes.c_double(1./seconds_per_point)  )  )
        CHK( dll.DAQmxSetCOPulseDutyCyc( self._co_task, self._square_wave_device, ctypes.c_double(duty_cycle)  )   )
        self._seconds_per_point = seconds_per_point
        self._duty_cycle = duty_cycle

    def getTiming(self):
        return self._seconds_per_point, self._duty_cycle

    def __del__(self):
        CHK( dll.DAQmxClearTask(self._ao_task)  )
        CHK( dll.DAQmxClearTask(self._co_task)  )


class HybridScannerTimeTaggerNI():

    def __init__(self, analog_channels, square_wave_device, time_tagger_serial, time_tagger_count_channel, time_tagger_marker_channel, seconds_per_point=1e-3, duty_cycle=0.5, voltage_range=(0.,10.), x_range=(0.,100.), y_range=(0.,100), z_range=(0.,20.)):
        self._pulsed_ao = PulsedAO(channels=analog_channels, square_wave_device=square_wave_device, seconds_per_point=seconds_per_point, duty_cycle=duty_cycle, range=voltage_range)
        TimeTagger._Tagger.setSerial(time_tagger_serial)
        self._N = None # we do not know the length of scan line yet
        self._x0 = x_range[0]
        self._x1 = x_range[1]
        self._y0 = y_range[0]
        self._y1 = y_range[1]
        self._z0 = z_range[0]
        self._z1 = z_range[1]
        self._position = numpy.zeros(3)
        self.setPosition(self._x0, self._y0, self._z0)
        self._time_tagger_count_channel=time_tagger_count_channel
        self._time_tagger_marker_channel=time_tagger_marker_channel
        
        self.overflow = TimeTagger.OverFlow()
                
    
    def getXRange(self):
        return self._x0, self._x1
    
    def getYRange(self):
        return self._y0, self._y1
    
    def getZRange(self):
        return self._z0, self._z1
    
    def setPosition(self, x, y, z):
        """Move stage to x, y, z"""
        self._pulsed_ao.output(self._convPosToVolt((x, y, z)))
        self._position[0] = x
        self._position[1] = y
        self._position[2] = z
    
    def getPosition(self):
        return self._position

    def _prepareCounter(self, N):
        if self._N == N:
            self._count_between_markers.clean()
        else:
            self._count_between_markers = TimeTagger.CountBetweenMarkers(self._time_tagger_count_channel, self._time_tagger_marker_channel, N)
            time.sleep(0.5)

    def initImageScan(self, rows, lines, seconds_per_point, return_speed=None):
        self.return_speed = return_speed
        self.seconds_per_point = seconds_per_point
        
        self.first_line = 1
        
        if return_speed:
            self._prepareCounter((rows+1)*lines*2-1)
            self.pixel_per_line = (rows+1)*2
        else:
            self._prepareCounter((rows+1)*lines-1)
            self.pixel_per_line = (rows+1)
    
    def doImageLine(self, line):
        if self.first_line:
            self._pulsed_ao.output( self._convPosToVolt(line[:,0]) )
            time.sleep(0.1)
            self.first_line = 0
        
        if self.return_speed:
            if (line.shape[1]+1)*2 != self.pixel_per_line:
                raise ValueError("analog out line length wrong") 
            self._pulsed_ao.setTiming(self.seconds_per_point, 0.5)
            self._pulsed_ao.output( self._convPosToVolt(line) )
            self._pulsed_ao.setTiming(self.seconds_per_point/float(self.return_speed), 0.5)
            self._pulsed_ao.output( self._convPosToVolt(line[:,::-1]) )
        else:
            if line.shape[1]+1 != self.pixel_per_line:
                raise ValueError("analog out line length wrong") 
            self._pulsed_ao.setTiming(self.seconds_per_point, 0.5)
            self._pulsed_ao.output( self._convPosToVolt(line) )
    
    def getImage(self, blocking=0):
        if blocking:
            timeout = 3.
            start_time = time.time()
            while not self._count_between_markers.ready():
                time.sleep(0.1)
                if time.time() - start_time > timeout:
                    print "count between markers timeout in scanner.py"
                    break
                
        image = numpy.append(self._count_between_markers.getData(),0).reshape((-1, self.pixel_per_line))  
        
        if self.return_speed is not None:
            imageout = image[:,:image.shape[1]/2-1].astype(float)            
        else:
            imageout = numpy.zeros((image.shape[0],image.shape[1]-1))
            for i in range(1,image.shape[0],2):
                imageout[i,:] = image[i,-2::-1]
            for i in range(0,image.shape[0],2):
                imageout[i,:] = image[i,0:-1:1]
    
        return imageout / float(self.seconds_per_point)


    def scanLine(self, line, seconds_per_point, return_speed=None):
        """Perform a line scan. If return_speed is not None, return to beginning of line
        with a speed 'return_speed' times faster than the speed currently set.
        """
        
        self._pulsed_ao.output( self._convPosToVolt(line[:,0]) )
        time.sleep(0.1)
        
        N = line.shape[1]
        self._prepareCounter(N)
        self._pulsed_ao.setTiming(seconds_per_point, 0.5)
        start_time = time.time()
        self._pulsed_ao.output( self._convPosToVolt(line) )
        
        timeout = 1.
        start_time = time.time()
        while not self._count_between_markers.ready():
            time.sleep(0.1)
            if time.time() - start_time > timeout:
                break

        if return_speed is not None:
            self._pulsed_ao.setTiming(seconds_per_point/float(return_speed), 0.5)
            self._pulsed_ao.output( self._convPosToVolt(line[:,::-1]) )

        return self._count_between_markers.getData() / float(seconds_per_point)

    def _convPosToVolt(self, Pos):
        x0 = self._x0
        x1 = self._x1
        y0 = self._y0
        y1 = self._y1
        z0 = self._z0
        z1 = self._z1
        return numpy.vstack( ( 10.0 / (x1-x0) * (Pos[0]-x0),
                               10.0 / (y1-y0) * (Pos[1]-y0),
                               10.0 / (z1-z0) * (Pos[2]-z0) ) )

    def Count(self):
        """Return a single count."""
        if not 'counter' in dir(self):
            self.counter = TimeTagger.Counter(self._time_tagger_count_channel,80000000)
        return self.counter.getData()[0]*10


if __name__ == '__main__':
    #PO = PulsedAO(channels='/Dev1/ao0:2', square_wave_device='/Dev1/Ctr0', seconds_per_point=1e-3, duty_cycle=0.5)
    scanner = HybridScannerTimeTaggerNI(analog_channels='/Dev1/ao0:2',
                                        square_wave_device='/Dev1/Ctr0',
                                        time_tagger_serial='VJaAMgsRxh',
                                        time_tagger_count_channel=0,
                                        time_tagger_marker_channel=1,
                                        seconds_per_point=1e-3,
                                        duty_cycle=0.5,
                                        voltage_range=(0.,10.),
                                        x_range=(0.,75.),
                                        y_range=(0.,75),
                                        z_range=(-25.,25.))
    


