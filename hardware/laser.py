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

from traits.api import HasTraits, Range
from traitsui.api import View, Item, RangeEditor

from hardware.nidaq import AOTask

class Laser( HasTraits ):

    def __init__(self, AO_channel='/Dev1/ao3', voltage_range=(0.,5.), **kwargs):
        self.AOTask = AOTask(Channels=AO_channel, range=voltage_range)
        self.AOTask.Write(np.array((float(voltage_range[0]),)))
        self.add_trait('voltage', Range(low=float(voltage_range[0]), high=float(voltage_range[1]), value=float(voltage_range[0]), desc='output voltage', label='Voltage [V]'))
        self.on_trait_change(self.write_voltage, 'voltage')
        HasTraits.__init__(self, **kwargs)
        
    def write_voltage(self, new):
        self.AOTask.Write(np.array((new,)))

    view = View(Item('voltage'),
                title='Laser', width=400, buttons=[], resizable=True)

if __name__=='__main__':
    laser = Laser()
    laser.edit_traits()
    