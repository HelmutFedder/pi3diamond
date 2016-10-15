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

'''
Created on 20.04.2012

author: Helmut Fedder
'''

import time
import numpy as np

# helper class to represent a visa instrument via a socket
class SocketInstrument():
    def __init__(self, device):
        import socket
        host,port = device.split(':')
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((host, port))
        self.sock = sock
        
    def write(self, cmd):
        """Sends a command over the socket"""
        cmd_string = cmd + '\n'
        sent = self.sock.sendall(cmd_string)
        if sent != None:
            raise RuntimeError('Transmission failed')
        time.sleep(.1) #add a timeout for the transfer to take place. Should be replaced by a proper error handling at some point
        
    def ask(self,question):
        """sends the question and receives the answer"""
        self.write(question)
        answer = self.sock.recv(2048)#2000
        return answer[:-2]
        
    def close(self):
        self.sock.close()


class HMP2030():

    def __init__(self,device, voltage_max=20.0, current_max=2.0, fuse_voltage_max=20.0, fuse_current_max=2.5):
        """
        Provides communication with a HMP2030 power supply
        via USB (virtual COM) or LAN.
        
        Usage Examples:
        
            hmp = HMP2030('ASRL11::INSTR')
            hmp = HMP2030('192.168.0.11:50000')
        
        Parameters:
        
            device:        string that describes the device (see examples above)
            
        Optional Parameters:
            voltage_max:   maximum allowed voltage
            current_max:   maximum allowed current
            fuse_voltage_max:    maximum allowed fuse voltage
            fuse_current_max:    maximum allowed fuse current
        """
        if '::' in device:
            self._connect_serial(device)
        else:
            self._connect_lan(device)
            
        self.voltage_max=voltage_max
        self.current_max=current_max
        self.fuse_voltage_max=fuse_voltage_max
        self.fuse_current_max=fuse_current_max

    def _connect_serial(self, device):
        import visa
        instr=visa.instrument('ASRL11::INSTR')
        instr.term_chars='\n'
        instr.chunk_size=4096
        instr.timeout=1.0
        self.instr = instr

    def _connect_lan(self, device):
        """connects to the hameg powersupply"""
        self.instr=SocketInstrument(device)    
    
    # convenience method
    def set_output(self,channel,current):
        """Set the current on the given channel. Turn output on or off depending on the specified current."""
        self.set_ch(channel)
        if current<=0 or current is None:
            self.stop()
        else:
            self.set_current(current)
            self.run()
    
    #functions to perform different SCPI-commands            
    def set_ch(self,ch):
        """sets the channel 1, 2 or 3"""
        if ch in [1,2,3]:
            self.instr.write('INST OUPT' + str(ch))
        else:
            raise ValueError('Wrong channel number. Chose 1, 2 or 3.')
    
    def get_ch(self):
        """asks for the selected channel"""
        channel = int(self.instr.ask('INST:NSEL?'))
        return channel
    
    def status(self,ch):
        """gets the current status of the selected channel (CC or CV)"""
        state = self.instr.ask('STAT:QUES:INST:ISUM' + str(ch) + ':COND?')
        if state == 1:
            return 'CC'
        elif state == 2:
            return 'CV'
        else:
            print "Couldn't read the status of the selected channel."
    
    def set_voltage(self,volt):
        """sets the voltage to the desired value"""
        if volt < 0:
            print 'The selected voltage cannot be set.'
        elif volt > self.voltage_max :     #the voltage_max will be set on the power supply if volt exceed voltage_max 
            self.instr.write('VOLT %1.3f' %self.voltage_max)
            print 'The set voltage exceed the maximum voltage: %1.3f' %self.voltage_max 
        else:
            self.instr.write('VOLT %1.3f' %volt)
     
    def set_voltage_step(self,vstep):
        """increases the voltage by a desired step"""
        vset = get_voltage()
        set_voltage(vset + vstep)
    
    def get_voltage(self):
        """measures the voltage"""
        voltage = float(self.instr.ask('MEAS:VOLT?'))
        return voltage
        
    def set_current(self,curr):
        """sets the current to the desired value"""
        if curr < 0:
            print 'The selected current cannot be set.'
        elif curr > self.current_max :    #the voltage_max will be set on the power supply if volt exceed voltage_max 
            self.instr.write('CURR %1.3f' %self.current_max)
            print 'The set current exceed the maximum current: %1.3f' %self.current_max 
        else:
            self.instr.write('CURR %1.3f' %curr)
    
    def set_current_step(self,cstep):
        """increases the current by a desired step"""
        cset = get_current()
        set_current(cset + cstep)
    
    def get_current(self):
        """measures the current"""
        current = float(self.instr.ask('MEAS:CURR?'))
        return current
        
    def set_arbitrary(self,ch, seq, N):
        """performs a sequence of voltage and current values for a given time on one channel with a number of repetitions.
           ch: channel for output
           seq: sequence to be set in form of a nested list = [(voltage,current,time),(..),(..),...]
           N: number of repetitions [1..255]. 0 means infinite repetitions."""
        seq_ary = np.array(seq)
        if max(seq_ary[:,0]) > self.voltage_max:
            print 'The set voltage exceed the maximum voltage: %1.3f' %self.voltage_max
        elif max(seq_ary[:,1]) > self.current_max:
            print 'The set current exceed the maximum current: %1.3f' %self.current_max
        elif min(seq_ary[:,2]) < .5:
            print 'The set time is shorter than 0.5s.'
        elif seq >= 0:
            print 'Negative value of voltage, current or time.'
        elif ch != [1,2,3]:
            print 'Wrong channel number. Chose 1, 2 or 3.'
        elif N != range(0,256):
            print 'The set repetitions are outside the range [0,255].'
        else:
            self.instr.write('ARB:DATA' + ' ' + str(seq).translate(None, '[()] '))
            self.instr.write('ARB:REP' + ' ' + str(N))
            self.instr.write('ARB:TRANS' + ' ' + str(ch))
            self.instr.write('ARB:STAR' + ' ' + str(ch))
            set_ch(ch)
            run()
    
    def stop_arbitrary(self,ch):
        """stops the arbitrary sequence of a specified channel ch, but leafs the output on."""
        self.instr.write('ARB:STOP' + ' ' + str(ch))
    
    def get_arbitrary(self,ch):
        """gets the number of performed repetitions of the arbitrary sequence"""
        set_ch(ch)
        num = int(self.instr.ask('ARB:REP?'))
        return num
        
    def get_all(self):
        """gets the measured values for all channels in the form [(ch,V,A),]"""
        l = []
        for i in [1,2,3]:
            set_ch(i)
            vset = get_voltage()
            cset = get_current()
            l.append((i,vset,cset))
        return l
    
    def run(self):
        """turns the output from the chosen channel on"""
        self.instr.write('OUTP ON')
        
    def run_all(self):
        """turns the output from all channels on"""
        set_ch(1)
        self.instr.write('OUTP:SEL ON')
        set_ch(2)
        self.instr.write('OUTP:SEL ON')
        set_ch(3)
        self.instr.write('OUTP:SEL ON')
        self.instr.write('OUTP:GEN ON')
    
    def stop(self):
        """turns the output from the chosen channel off"""
        self.instr.write('OUTP OFF')
            
    def stop_all(self):
        """stops the output of all channels"""
        set_ch(1)
        self.instr.write('OUTP:SEL OFF')
        set_ch(2)
        self.instr.write('OUTP:SEL OFF')
        set_ch(3)
        self.instr.write('OUTP:SEL OFF')
        self.instr.write('OUTP:GEN OFF')
    
    def start(self):
        """starts up the whole system"""
        self.instr.write('*RST') #resets the device
        self.instr.write('SYST:REM') #sets the instrument to remote control
        
    def close(self):
        """stops and disconnects the device"""
        stop_all()
        self.instr.close()
        
    def beep(self):
        """gives an acoustical signal from the device"""
        self.instr.write('SYST:BEEP')
        
    def error_list(self):
        """prints all errors from the error register."""
        error = str(self.instr.ask('SYST:ERR?'))
        return error
    
    def OVP(self,fuse_voltage_max):
        """sets the Over-Voltage-Protection to the value fuse_voltage_max for a selected channel"""
        if fuse_voltage_max < 0:
            print 'The selected value for voltage protection cannot be set.'
        elif fuse_voltage_max > 32.0: #the maximal voltage which the HMP2030 supplies
            print 'The set voltage exceed the maximum voltage: 32V' 
        else:
            self.instr.write('VOLT:PROT %1.3f' %fuse_voltage_max)
        
    def FUSE(self,fuse_current_max):
        """sets the fuse to the value fuse_current_max and the delay time to 0ms for a selected channel"""
        self.instr.write('FUSE ON')
        if fuse_current_max < 0:
            print 'The selected value for current fuse cannot be set.'
        elif fuse_current_max > 5.0: #the maximal current which the HMP2030 supplies
            print 'The set current exceed the maximum current: 5A' 
        else:
            self.instr.write('CURR %1.3f' %fuse_current_max)
        self.instr.write('FUSE:DEL 0')


class BipolarHMP2030(HMP2030):
    #funktions to reverse the polarity
    def set_polarity(self,ch,p):
        """sets the polarity p of a given channel"""
        pass
    
    def get_polarity(self,ch):
        """gets the polarity p of a given channel"""
        pass

from traits.api import HasTraits, Range
from traitsui.api import View, Item

class HMP2030Traits( HMP2030, HasTraits ):

    def __init__(self, device, voltage_max=20.0, current_max=2.0, fuse_voltage_max=20.0, fuse_current_max=2.5, **kwargs):
        HMP2030.__init__(self, device, voltage_max, current_max, fuse_voltage_max, fuse_current_max)
        self.add_trait('current', Range(low=0.0, high=current_max, value=self.get_current(), label='Current [A]', auto_set=False, enter_set=True))
        HasTraits.__init__(self, **kwargs)

    def _current_changed(self, new):
        self.set_output(1, new)

    traits_view = View(Item('current'),title='HMP2030')


#-------------------------------------------------------------------------------------------------------------------------

#define the sub-function including:set_channel, set_voltage, set_current and run    
#def set(ch):
#    write('SYST REM')
#    write('INST OUPT1')
#    write('OUTP OFF')

#
#write('SYST REM')   #to remote control 
#write('INST OUPT1') #select channel 1
#stop()
