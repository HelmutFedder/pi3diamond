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

import numpy as np
import time
import ctypes
from nidaq import dll, CHK, DAQmx_Val_Volts, DAQmx_Val_Hz, DAQmx_Val_Low, DAQmx_Val_GroupByChannel, DAQmx_Val_OnDemand, DAQmx_Val_SampClk, DAQmx_Val_Rising, DAQmx_Val_FiniteSamps, DOTask
from TimeTagger import CountBetweenMarkers, TimeDifferences, CHANNEL_INVALID

class AnalogOut():
    """Provide burst and steady state analog output with a national instrument card."""

    def __init__(self, analog_channels='/Dev1/ao0:2', output_counter='/Dev1/Ctr0', voltage_range=(10,10)):

        self._analog_channels = analog_channels
        self._output_counter = output_counter
        self._voltage_range=voltage_range

        # we start with default duty cycle
        self._duty_cycle = 0.5
        self._seconds_per_point = 1000.
        self._transfer_timeout = 1.0
        
        # we start with unknown length
        self._length = None

        # create nidaq tasks
        self._n_written = ctypes.c_long()        
        self._ao_task = ctypes.c_ulong()
        self._co_task = ctypes.c_ulong()
        CHK( dll.DAQmxCreateTask('', ctypes.byref(self._ao_task)) )
        CHK( dll.DAQmxCreateTask('', ctypes.byref(self._co_task))  )

        CHK( dll.DAQmxCreateAOVoltageChan(self._ao_task,
                                          self._analog_channels,
                                          '',
                                          ctypes.c_double(self._voltage_range[0]),
                                          ctypes.c_double(self._voltage_range[1]),
                                          DAQmx_Val_Volts,
                                          '')  )

        CHK(  dll.DAQmxCreateCOPulseChanFreq( self._co_task,
                                              self._output_counter, '',
                                              DAQmx_Val_Hz, DAQmx_Val_Low, ctypes.c_double(0),
                                              ctypes.c_double(self._seconds_per_point),
                                              ctypes.c_double(self._duty_cycle) )  )

    def set_output(self, data):
        """If len(data)==1, output one point. Else output sequence of points.
        The shape of data is (n_channels,n_points)."""
        if data.ndim == 1:
            N = 1
        else:
            N = data.shape[1]
        self.set_length(N)
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

    def set_length(self, N):
        if self._length == N:
            return
        if N == 1:
            CHK( dll.DAQmxSetSampTimingType( self._ao_task, DAQmx_Val_OnDemand)  )
        else:
            CHK( dll.DAQmxSetSampTimingType( self._ao_task, DAQmx_Val_SampClk)  )
            CHK( dll.DAQmxCfgSampClkTiming( self._ao_task,
                                            self._output_counter+'InternalOutput',
                                            ctypes.c_double(1./self._seconds_per_point),
                                            DAQmx_Val_Rising,
                                            DAQmx_Val_FiniteSamps,
                                            ctypes.c_ulonglong(N)) )
            CHK(  dll.DAQmxCfgImplicitTiming( self._co_task, DAQmx_Val_FiniteSamps, ctypes.c_ulonglong(N+1))  )
        self._length = N

    def get_length(self):
        return self._length

    def set_timing(self, seconds_per_point, duty_cycle=0.5):
        if seconds_per_point == self._seconds_per_point and duty_cycle == self._duty_cycle:
            return
        CHK( dll.DAQmxSetCOPulseFreq( self._co_task, self._output_counter, ctypes.c_double(1./seconds_per_point)  )  )
        CHK( dll.DAQmxSetCOPulseDutyCyc( self._co_task, self._output_counter, ctypes.c_double(duty_cycle)  )   )
        self._seconds_per_point = seconds_per_point
        self._duty_cycle = duty_cycle

    def get_timing(self):
        return self._seconds_per_point, self._duty_cycle

    def __del__(self):
        CHK( dll.DAQmxClearTask(self._ao_task)  )
        CHK( dll.DAQmxClearTask(self._co_task)  )

        
class Scanner( AnalogOut ):
    
    def __init__(self,
                 analog_channels='/Dev1/ao0:2',
                 output_counter='/Dev1/Ctr0',
                 x_range=(0.,100.),
                 y_range=(0.,100.),
                 z_range=(-10.,10.),
                 voltage_range=(0.,10.),
                 invert_x=False, invert_y=False, invert_z=False, swap_xy=False,
                 TriggerChannels=None):
        
        AnalogOut.__init__(self,
                           analog_channels=analog_channels,
                           output_counter=output_counter,
                           voltage_range=voltage_range)

        if TriggerChannels is not None:
            self._trigger_task = DOTask(TriggerChannels)
        self._x_range = x_range
        self._y_range = y_range
        self._z_range = z_range
        self._v_range = voltage_range 
        self._invert_x = invert_x
        self._invert_y = invert_y
        self._invert_z = invert_z
        self._swap_xy = swap_xy

        self._position = np.zeros(3)

        self.set_position(x_range[0], y_range[0], z_range[0])
        
    def get_x_range(self):
        return self._x_range
        
    def get_y_range(self):
        return self._y_range

    def get_z_range(self):
        return self._z_range

    def set_position(self, x, y, z):
        """Move to x, y, z"""
        self.set_output(self.pos_to_volt((x, y, z)))
        self._position[:] = (x, y, z)

    def get_position(self):
        """Return the current position""" 
        return self._position

    def scan_line(self, line):
        """
        Move along a line.
        """                
        if hasattr(self, '_trigger_task'):
            self._trigger_task.Write(np.array((1,0), dtype=np.uint8) )
            time.sleep(0.001)
            self._trigger_task.Write(np.array((0,0), dtype=np.uint8) )
        self.set_output(self.pos_to_volt(line))
        self._position[:] = line[:,-1]

    def pos_to_volt(self, r):
        x = self._x_range
        y = self._y_range
        z = self._z_range
        v = self._v_range
        v0 = v[0]
        dv = v[1]-v[0]
        if self._invert_x:
            vx = v0+(x[1]-r[0])/(x[1]-x[0])*dv            
        else:
            vx = v0+(r[0]-x[0])/(x[1]-x[0])*dv            
        if self._invert_y:
            vy = v0+(y[1]-r[1])/(y[1]-y[0])*dv            
        else:
            vy = v0+(r[1]-y[0])/(y[1]-y[0])*dv            
        if self._invert_z:
            vz = v0+(z[1]-r[2])/(z[1]-z[0])*dv            
        else:
            vz = v0+(r[2]-z[0])/(z[1]-z[0])*dv
        if self._swap_xy:
            vt = vx
            vx = vy
            vy = vt            
        return np.vstack( (vx,vy,vz) )


"""    
    def pos_to_volt(self, pos):
        v0, v1 = self._voltage_range
        x0, x1 = self._x_range
        y0, y1 = self._y_range
        z0, z1 = self._z_range
        return np.vstack( ( v0 + (v1-v0) / (x1-x0) * (pos[0]-x0),
                               v0 + (v1-v0) / (y1-y0) * (pos[1]-y0),
                               v0 + (v1-v0) / (z1-z0) * (pos[2]-z0) ) )
"""

class TimeTaggerImager():
    """
    An imager based on a scanner and a time tagger
    """
    
    def __init__(self, scanner, tagger, click_channel, begin_channel):
        self.scanner = scanner
        self.tagger = tagger
        self.click_channel = click_channel
        self.begin_channel = begin_channel
        self._n = None

    def scan_line(self, line, seconds_per_point):
        n = line.shape[1]
        if n == self._n:
            self._count.clear()
        else:
            self._n = n
            self._count = CountBetweenMarkers(self.time_tagger, self.click_channel[0], self.begin_channel, CHANNEL_INVALID, n)
        self.scanner.set_timing(seconds_per_point, 0.5)
        ext_line = np.append(line, line[:,-1].reshape((3,1)),axis=1)
        self.scanner.scan_line(ext_line)
        timeout = 2*seconds_per_point*n
        start_time = time.time()
        while not self._count.ready() and start_time < timeout:
            time.sleep(timeout/20.)
        return self._count.getData() / float(seconds_per_point)

    def set_position(self, *args, **kwargs):
        """Move to x, y, z"""
        self.scanner.set_position(*args, **kwargs)

    def get_position(self):
        """Return the current position""" 
        return self.scanner.get_position()

    def set_timing(self, *args, **kwargs):
        self.scanner.set_timing(*args, **kwargs)

    def get_timing(self):
        """Return the current position""" 
        return self.scanner.get_timing()
    def get_x_range(self):
        return self.scanner.get_x_range()

    def get_y_range(self):
        return self.scanner.get_y_range()

    def get_z_range(self):
        return self.scanner.get_z_range()
    
    def new_image(self):
        pass


class HistogramImager( Imager ):
    """Acquire a multidimensional histogram in every pixel."""
    
    def __init__(self, scanner, time_tagger, start_channel=2, next_channel=4, binwidth=5000000, n_bins=2):
        self.scanner = scanner
        self.time_tagger = time_tagger
        self.start_channel = start_channel
        self.next_channel = next_channel
        self.binwidth = binwidth
        self.n_bins = n_bins
        self._n = None
        
    def scan_line(self, line, seconds_per_point):
        n = line.shape[1]
        
        args_0 = (self.time_tagger,
                  0,
                  self.start_channel,
                  self.next_channel,
                  4294967295,
                  self.binwidth,
                  self.n_bins,
                  n)
        
        args_1 = (self.time_tagger,
                  1,
                  self.start_channel,
                  self.next_channel,
                  4294967295,
                  self.binwidth,
                  self.n_bins,
                  n)
        
        if n == self._n:
            self._count_0.clear()
            self._count_1.clear()
        else:
            self._n = n
            self._count_0 = TimeDifferences(*args_0)
            self._count_1 = TimeDifferences(*args_1)
        self.scanner.set_timing(seconds_per_point, 0.5)
        ext_line = np.append(line, line[:,-1].reshape((3,1)),axis=1)
        self.scanner.scan_line(ext_line)
        timeout = 2*seconds_per_point*n
        start_time = time.time()
        while self._count_0.getCounts() < 1 and self._count_1.getCounts() < 1 and time.time() - start_time < timeout:
            time.sleep(timeout/20.)
        return self._count_0.getData() +self._count_1.getData()

class DifferenceImager( HistogramImager ):
    """Acquire a multidimensional histogram in every pixel."""
            
    def scan_line(self, *args, **kwargs):
        return self.process_data(HistogramImager.scan_line(self, *args, **kwargs))

    def process_data(self, data):
        n = data.shape[1]/2
        return (data[:,:n].sum(axis=1)-data[:,n:].sum(axis=1))/self.scanner._seconds_per_point

class FourPointImager( HistogramImager ):
    """Acquire a multidimensional histogram in every pixel."""
            
    def scan_line(self, *args, **kwargs):
        return self.process_data(HistogramImager.scan_line(self, *args, **kwargs))

    def process_data(self, data):
        n = data.shape[1]/4
        newline = data.reshape(len(data),4,n).sum(axis=2)
        self.image.append(newline)
        return np.sign(newline[:,0]-newline[:,2])*np.sqrt(np.abs((newline[:,0]-newline[:,2])*(newline[:,0]-newline[:,1])))/(self.scanner._seconds_per_point)
    
    def new_image(self):
        self.image = []

if __name__ == '__main__':

    """    
    
    binwidth = hold*1000
    n_bins = 2
    
    args = (tagger,
            0,
            2,
            4,
            4294967295,
            binwidth,
            n_bins,
            2)
    """

    hold = 100000
    voltage = 3
    microwave.setOutput(-14,2.87e9)
    
    diff_imager = DifferenceImager(scanner, tagger, binwidth=hold*1000, n_bins=2)
      
    confocal.imager = diff_imager
  
#     fourpoint_imager = FourPointImager(scanner, tagger, binwidth=hold*1000, n_bins=4)
#      
#     confocal.imager = fourpoint_imager
    
    """
    # ODMR
    sequence = 100*[(['aom','detect'],100,0),
                    (['aom'],hold,0),
                    (['aom','microwave'],hold,0), ]
    """
    
    sequence = 100*[(['aom','detect'],100,voltage),
                    (['aom'],hold,voltage),
                    (['aom'],hold,-voltage), ]
   
#     sequence = 100*[(['aom','detect'],100,voltage),
#                     (['aom'],hold,voltage),
#                     (['aom','microwave'],hold,voltage), ]
    
#     sequence = 100*[(['aom','detect'],100,voltage),
#                     (['aom'],hold,voltage),
#                     (['aom','microwave'],hold,voltage),
#                     (['aom'],hold,-voltage),
#                     (['aom','microwave'],hold,-voltage),  ]
    
    analog_pulser.setSequence(sequence)
   
    #td = TimeDifferences(*args)
    
    """
    scanner = Scanner(analog_channels='/Dev1/ao0:2',
                      output_counter='/Dev1/Ctr1',
                      x_range=(0.,20.),
                      y_range=(0.,20.),
                      z_range=(0.,20.),
                      voltage_range=(0.,10.),
                      invert_y=True
                      )    
    
    xline = np.arange(1.,10.,1.)
    yline = np.arange(1.,10.,1.)
    zline = np.arange(1.,10.,1.)
    
    line = np.vstack((xline,yline,zline))
    
    scanner.set_timing(0.01, 0.5)
    scanner.scan_line(line)
    
    from TimeTagger import createTimeTagger
    tagger = createTimeTagger(serial='12220003BR')

    imager = Imager(scanner, tagger, 0, 4)
    imager.scan_line(line, 0.1)
    """