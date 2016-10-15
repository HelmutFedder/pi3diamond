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

import logging, logging.handlers
#import os
#import inspect

#path = os.path.dirname(inspect.getfile(inspect.currentframe()))

# start the logger
#file_handler = logging.handlers.TimedRotatingFileHandler(path+'/log/diamond_log.txt', 'W6') # start new file every sunday, keeping all the old ones 
#file_handler.setFormatter(logging.Formatter("%(asctime)s - %(module)s.%(funcName)s - %(levelname)s - %(message)s"))
#file_handler.setLevel(logging.DEBUG)
stream_handler=logging.StreamHandler()
stream_handler.setLevel(logging.INFO) # we don't want the console to be swamped with debug messages
#logging.getLogger().addHandler(file_handler)
logging.getLogger().addHandler(stream_handler) # also log to stderr
logging.getLogger().setLevel(logging.DEBUG)
logging.getLogger().info('Starting logger.')

# start the JobManager (it is currently a singleton)
from tools import emod
emod.JobManager().start()

# start the CronDaemon
#from tools import cron
#cron.CronDaemon().start()

# numerical classes that are used everywhere
import numpy as np 

#########################################
# hardware
#########################################
from hardware.mocking import Imager
imager = Imager()

#from TimeTagger import createMockingTagger
#tagger = createMockingTagger()

#########################################
# example how to define a little panel
# for laser and mirror hardware controls
#########################################

#from traits.api import HasTraits, Instance
#from traitsui.api import View, Item, Group

#class HardwarePanel( HasTraits ):
#    
#    laser   = Instance( hardware.laser.Laser )
#    mirror  = Instance( hardware.flip_mirror.FlipMirror )
#    stage   = Instance( hardware.apt_stage.RotationStageTraits )
#    coil    = Instance( hardware.hameg.HMP2030Traits )
#    
#    traits_view = View(Group(Item('laser', style='custom', show_label=False),
#                             Item('mirror', style='custom', show_label=False),
#                             Item('stage', style='custom', show_label=False),
#                             Item('coil', style='custom', show_label=False),
#                             ),
#                       title='Hardware Panel'
#                       )
#    
#hardware_panel = HardwarePanel(laser=laser, mirror=flip_mirror, stage=rotation_stage, coil=coil)
#hardware_panel.edit_traits()


#########################################
# create measurement widgets
#########################################
from measurements.confocal import Confocal # a confocal scanner imaging tool
from measurements.auto_focus import AutoFocus # auto focus on objects in the image

confocal = Confocal(imager)
confocal.edit_traits()

auto_focus = AutoFocus(imager, confocal)
auto_focus.edit_traits()
