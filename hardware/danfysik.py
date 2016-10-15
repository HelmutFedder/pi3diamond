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

import visa

class DanFysik:
    
    def __init__(self, device='GPIB0::7::INSTR', limit_current=10.0, maximum_current=20.0):
        self.maximum_current=maximum_current
        self.limit_current=limit_current
        instr = visa.instrument(device)
        instr.timeout=1.0
        instr.term_chars='\r'
        self.instr=instr
        
    def set_current(self, current):
        """Sets the current in Ampere."""
        if current > self.limit_current:
            raise ValueError('Limit current exceeded.')
        ppm = int(current / self.maximum_current * 1000000)
        if ( ppm == 0 ):
            self.instr.write('DA 0,%i'%ppm)
            self.instr.write('F')
        else:
            self.instr.write('DA 0,%i'%ppm)
            self.instr.write('N')
            
    def get_current(self):
        s = self.instr.ask('DA 0')
        return float(s.split()[1])*1e-6*self.maximum_current

