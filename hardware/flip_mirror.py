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

from traits.api import HasTraits, Button
from traitsui.api import View, Item

from enthought.traits.ui.menu import Action, Menu, MenuBar

from nidaq import DOTask

class FlipMirror( HasTraits ):

    flip = Button(desc="flip the mirror")
    
    def __init__(self, trigger_channels):
        super(HasTraits, self).__init__()
        self._trigger_task = DOTask(trigger_channels)

    def _flip_changed(self):
        self._trigger_task.Write(np.array((1,), dtype=np.uint8) )
        time.sleep(0.001)
        self._trigger_task.Write(np.array((0,), dtype=np.uint8) )

    view = View(Item('flip', show_label=False),
                title='Flip mirror',
                buttons=[],
                resizable=True)
    