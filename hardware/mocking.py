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
import logging
import time

class Imager():
    """
    A confocal imager.

    To acquire a confocal image, we use a piezo scanner
    that moves the laser focus in x-y-z,
    and some sort of photon counter that detects the
    intensity at every pixel.
    
    We call the combination of the two an 'Imager'.
    """
    def get_x_range(self):
        return (0.,100.)
    def get_y_range(self):
        return (0.,100.)
    def get_z_range(self):
        return (-20.,20.)
    def set_x(self, x):
        pass
    def set_y(self, y):
        pass
    def set_z(self, z):
        pass
    def set_position(self, x, y, z):
        """Move stage to x, y, z"""
        pass
    def scan_line(self, Line, SecondsPerPoint, return_speed=None):
        """
        Move the scanner along a Line at a velocity defined by 'SecondsPerPoint'
        and return the intensity at every pixel (in counts/s).
        """
        time.sleep(0.1)
        return (1000*np.sin(Line[0,:])*np.sin(Line[1,:])*np.exp(-Line[2,:]**2)).astype(int)
    def new_image(self):
        """
        Helper method that is used in special imaging modes where image data is stored
        at the imager. In this case, raw data that is generated with successive line sweeps
        is accumulated at the imager and can be accessed later. The accumuated data is cleared
        with this method.
        """
        pass
        
class Sweeper(  ):
    def configure(self, n, SecondsPerPoint, DutyCycle=0.8):
        x = np.arange(n)
        a = 100.
        c = 50.
        x0 = n/2.
        g = n/10.
        y = np.int32( c - a / np.pi * (  g**2 / ( (x-x0)**2 + g**2 )  ) )
        self._sweeps = 0
        self._y = y
    def run(self):
        time.sleep(1)
        self._sweeps+=1
        return np.random.poisson(self._sweeps*self._y)
    def clear(self):
        pass

class Microwave(  ):
    def setPower(self, power):
        logging.getLogger().debug('Setting microwave power to '+str(power)+'.')
    def setOutput(self, power, frequency):
        logging.getLogger().debug('Setting microwave to p='+str(power)+' f='+str(frequency)+'.')
    def initSweep(self, f, p):
        logging.getLogger().debug('Setting microwave to sweep between frequencies %e .. %e with power %f.'%(f[0],f[-1],p[0]))
    def resetListPos(self):
        pass

class PulseGenerator():
    def Sequence(self, sequence, loop=True):
        pass
    def Light(self):
        pass
    def Night(self):
        pass
    def Open(self):
        pass
    def High(self):
        pass
    def checkUnderflow(self):
        return False
        #return np.random.random()<0.1

class Laser():
    """Provides control of the laser power."""
    voltage = 0.

class PowerMeter():
    """Provides an optical power meter."""
    power = 0.
    def getPower(self):
        """Return the optical power in Watt."""
        PowerMeter.power += 1
        return PowerMeter.power*1e-3

class Coil():
    pass

class RotationStage():
    def set_angle(self, angle):
        pass
        
def createMockingTagger( ):
    """
    Mocking factory function class for the TimeTagger
    """
    return None
    
class TimeDifferences( ):
    """
    TimeDifferences mocking class
    
    Uncomment the desired code to simulate data for Rabi oscillations, Hahn-Echo measurements, etc.    
    """
    def __init__(self, time_tagger,
                 click_channel, start_channel, next_channel, sync_channel,
                 bin_width, n_bins, n_histograms):
        self.n_bins = n_bins
        self.n_histograms = n_histograms
        self.clear()
                    
    def clear(self):
        n_histograms = self.n_histograms
        n_bins = self.n_bins       
        data = np.zeros((n_histograms,n_bins))
        m0 = int(n_bins/5)
        m = float(n_bins-m0)
        M = np.arange(m, dtype=float)
        n = float(n_histograms)
        k = n_histograms/2
        for i in range(n_histograms):
            """Rabi Data"""
            data[i,m0:] = 30*np.cos(3*2*np.pi*i/n)*np.exp(-5*M/m)+100
            """Hahn Data        
            data[i,m0:] = 30*np.exp(-9*i**2/n**2)*np.exp(-5*M/m)+100
            """
            """Hahn 3pi2 Data
            if i < k:
                data[i,m0:] = 30*np.exp(-9*i**2/float(k**2))*np.exp(-5*M/m)+100
            else:
                data[i,m0:] = -30*np.exp(-9*(i-k)**2/float(k**2))*np.exp(-5*M/m)+100
            """
            """T1 Data
            data[i,m0:] = 30*np.exp(-3*i/n)*np.exp(-5*M/m)+100
            """
        self.data = data
        self.counter = 1
    def setMaxCounts(self,arg):
        pass
    def ready(self):
        time.sleep(0.1)
        return True
    def getData(self):
        self.counter += 1
        return np.random.poisson(self.counter*self.data)
    def getCounts(self):
        return self.counter
    def start(self):
        pass
    def stop(self):
        pass

class Countrate( ):
    def __init__(self, channel):
        self.rate = 0.
    def getData(self):
        self.rate += 1.
        return 1e5/(1+20./self.rate)
    def clear(self):
        pass

#        Counter(self.tagger, self.channels, int(self.seconds_per_point*1e12), self.number_of_points)
        
class Counter( ):
    """
    Counter mocking class.
    """
    def __init__(self, tagger, channels, binwidth, n_values):
        self.n_channels = len(channels)
        self.seconds_per_point = float(binwidth)/800000000
        self.n_values = n_values
        self.binwidth = binwidth
    def getData(self):
        return np.random.random_integers(100000,120000, (self.n_channels,self.n_values))*self.seconds_per_point
    def getIndex(self):
        return np.arange(self.n_values)*self.binwidth
    def clear(self):
        pass

