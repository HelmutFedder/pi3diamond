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
import time
import struct

dll=ctypes.WinDLL('ftd2xx.dll')

class RotationStage():
    """Provide control over a thorlabs rotation stage."""

    def __init__(self, serial='83838311'):
        self._open_usb(serial)

    def _open_usb(self, serial):
        """Open a USB connection to the controller with the specified serial number."""
    
        handle=ctypes.c_ulong()

        if dll.FT_OpenEx(serial, 1, ctypes.byref(handle)):
            raise RuntimeError('USB open failed.')
        if dll.FT_SetBaudRate(handle ,ctypes.c_ulong(115200)):
            raise RuntimeError('USB baud rate failed.')
        if dll.FT_SetDataCharacteristics(handle, ctypes.c_uint8(8), ctypes.c_uint8(0), ctypes.c_uint8(0)):
            raise RuntimeError('USB set data format failed.')
        time.sleep(0.1)
        if dll.FT_Purge(handle, ctypes.c_ulong(1|2)):
            raise RuntimeError('USB purge failed.')
        time.sleep(0.1)
        if dll.FT_ResetDevice(handle):
            raise RuntimeError('USB reset failed.')
        if dll.FT_SetFlowControl(handle, ctypes.c_ushort(0x0100), ctypes.c_uint8(0), ctypes.c_uint8(0)):
            raise RuntimeError('USB set flow control failed.')
        if dll.FT_SetRts(handle):
            raise RuntimeError('USB set RTS failed.')
        if dll.FT_SetTimeouts(handle,1000,0):
            raise RuntimeError('USB set timeouts failed.')

        self.usb = handle

    def _close_usb(self):
        """Close the USB connection to the controller."""
        if dll.FT_Close(self.usb):
            raise RuntimeError('USB close connection failed.')

    def get_info(self):
        """Read the device information from the controller USB chip."""
        dev_type=ctypes.c_ulong()
        dev_id=ctypes.c_ulong()
        ser_buf = ctypes.create_string_buffer(16)
        desc_buf = ctypes.create_string_buffer(64)
    
        if dll.FT_GetDeviceInfo(self.usb, ctypes.byref(dev_type), ctypes.byref(dev_id), ser_buf, desc_buf, None):
            raise RuntimeError('USB get device info failed.')
            
        return dev_type.value, dev_id.value, ser_buf.value, desc_buf.value

    def send(self,command):
        """Send a 'set'-type command to the controller."""
        written = ctypes.c_ulong()
        wr_buf = ctypes.create_string_buffer(command,len(command))
        if dll.FT_Write(self.usb, wr_buf, ctypes.c_ulong(len(wr_buf)), ctypes.byref(written)):
            raise RuntimeError('USB write failed.')
        return written.value

    def wait(self,n,timeout=0.1):
        """Wait until n bytes are in the receive buffer or timeout is passed."""
        amount_rx=ctypes.c_ulong(0)
        amount_tx=ctypes.c_ulong(0)
        event_stat=ctypes.c_ulong(0)
        start_time = time.time()
        while amount_rx.value<n and time.time()-start_time<timeout:
            if dll.FT_GetStatus(self.usb, ctypes.byref(amount_rx), ctypes.byref(amount_tx), ctypes.byref(event_stat)):
                raise RuntimeError('USB get status failed.')
            time.sleep(0.1)
        return amount_rx.value
    
    def ack(self):
        self.send('\x92\x04\x00\x00\x05\x01')

    def request(self,command,timeout=0.1):
        """Send a 'get'-type command to the controller and return the reply."""
        written = ctypes.c_ulong()
        read = ctypes.c_ulong()
        wr_buf = ctypes.create_string_buffer(command,len(command))
        time.sleep(0.1)
        if dll.FT_Write(self.usb, wr_buf, ctypes.c_ulong(len(wr_buf)), ctypes.byref(written)):
            raise RuntimeError('USB write failed.')
        if self.wait(6,timeout) < 6: # wait for 6 bytes reply message
            raise RuntimeError('USB request timed out after %fs.'%timeout)
        head_buf = ctypes.create_string_buffer(6)
        if dll.FT_Read(self.usb, head_buf, ctypes.c_ulong(6), ctypes.byref(read)):
            raise RuntimeError('USB read failed.')
        id,length,dest,source=struct.unpack('<2H2B',head_buf.raw) # convert message assuming that it is the header for a data
        if dest>>7: # a data packet will follow
            if self.wait(length,timeout) < length: # wait for expected number of bytes
                raise RuntimeError('USB request timed out after %fs.'%timeout)
            rd_buf = ctypes.create_string_buffer(length)
            if dll.FT_Read(self.usb, rd_buf, ctypes.c_ulong(length), ctypes.byref(read)):
                raise RuntimeError('USB read failed.')
            self.ack()
            return head_buf.raw+rd_buf.raw
        else: # simple 6 byte message
            self.ack()
            return head_buf.raw
        
#        if dll.FT_GetStatus(self.usb, ctypes.byref(amount_rx), ctypes.byref(amount_tx), ctypes.byref(event_stat)):
#            raise RuntimeError('USB get status failed.')
#        time.sleep(0.1)
#        print 'rx', amount_rx.value
#        #print 'tx', amount_tx.value
#        rd_buf = ctypes.create_string_buffer(amount_rx.value)
#        #print 'len(rd_buf)', len(rd_buf)
#        time.sleep(0.1)
#        if dll.FT_Read(self.usb, rd_buf, amount_rx, ctypes.byref(read)):
#            raise RuntimeError('USB read failed.')
#        if read.value != amount_rx.value:
#            raise RuntimeError('USB requested %i bytes but received %i.'%(amount_rx.value,read.value))
        #print read.value
        
#        print rd_buf.value.__repr__()
        
#        id,length,dest,source = struct.unpack('<2H2B',rd_buf.value)
#        print length
#        print dest
#        if dest>>7:
#            print length
#            time.sleep(0.1)
#            rd_buf_dat = ctypes.create_string_buffer(length)
#            if dll.FT_Read(self.usb, rd_buf_dat, ctypes.c_ulong(length), ctypes.byref(read)):
#                raise RuntimeError('USB read failed.')
#            print 'return message and data'
#            return rd_buf.value+rd_buf_dat.value
#        else:
#            print 'return message'
#            return rd_buf.value

    def go_home(self, timeout=60.):
        #start_time=time.time()
        return self.request('\x43\x04\x01\x00\x05\x01', timeout=timeout) # worst case homing takes about 45 seconds with standard speed parameters
        #return time.time()-start_time
        
## not available with TDC001 controller
#    def get_stage_info(self):
#        buf=self.request('\xf1\x04\x01\x00\x05\x01')
#        return struct.unpack('<6s3H16s7i24s',buf)[1:-1]

    def blink(self):
        self.send('\x23\x02\x00\x00\x05\x01')

    def get_status(self):
        buf=self.request('\x90\x04\x01\x00\x05\x01')
        chan, pos, count, status = struct.unpack('<6sH3L',buf)[1:]
        return chan, pos, count, status

    def get_status_bits(self):
        buf=self.request('\x29\x04\x01\x00\x05\x01')
        status = struct.unpack('<6sHL',buf)[2]
        return status

    def get_backlash(self):
        buf=self.request('\x3b\x04\x01\x00\x05\x01')
        return struct.unpack('<6sHL',buf)[1:]

    def get_pid(self):
        buf=self.request('\xa1\x04\x01\x00\x05\x01')
        return struct.unpack('<6sH4LH',buf)[1:]

#    def set_angle(self,angle):
#        """Set the absolute angle [deg]."""
#        counts=1919.641857862339*(angle%360.) # pitch (deg / rev)=17.87, (counts / rev)=34304
#        self.send('\x53\x04\x06\x00\x85\x01'+'\x01\x00'+struct.pack('<l',counts))

    def set_angle(self,angle, timeout=60.):
        """Set the absolute angle [deg]."""
        #start_time=time.time()
        counts=1919.641857862339*(angle%360.) # pitch (deg / rev)=17.87, (counts / rev)=34304
        return self.request('\x53\x04\x06\x00\x85\x01'+'\x01\x00'+struct.pack('<l',counts), timeout=timeout)
        #return time.time()-start_time

from traits.api import HasTraits, Str, Float, Button
from traitsui.api import View, Item

class RotationStageTraits( RotationStage, HasTraits ):
    
    home = Button()
    angle = Float(default_value=0.0, label='angle [deg]', auto_set=False, enter_set=True)
    
    def __init__(self, serial, **kwargs):
        HasTraits.__init__(self, **kwargs)
        RotationStage.__init__(self, serial)

    def _home_changed(self):
        self.go_home()
        
    def _angle_changed(self, new):
        self.set_angle(new)

    traits_view = View(Item('angle'),
                       Item('home', show_label=False),
                       title='Rotation Stage'
                       )

if __name__=='__main__':

    stage=RotationStage()
    #d=stage.move_home()
    
    #print stage.get_pid()
    
    
    #stage.set_angle(0)

