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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diamond. If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2009-2016 Helmut Fedder <helmut@fedder.net>
"""

"""
There are several distinct ways to go through different NVs and perform
certain measurement tasks. 

1. using the queue and 'SwitchTarget' and 'SaveJob' and 'SetJob'

  For each task, create a job and submit it to the queue.
  Provide a 'special' job for switching the NV. I.e., a queue might
  look like this: [ ODMR, Rabi, SwitchTarget, ODMR, Rabi, SwitchTarget, ...]

  pro: - very simple
       - a different set of Jobs can be submitted for individual NVs
       - every part of the 'code' is basically tested separately (uses only
         existing jobs) --> very low chance for errors
       - queue can be modified by user on run time, e.g., if an error in the tasks
         is discovered, it can be corrected
       - the submitted jobs can be run with lower priority than all the usual
         jobs, i.e., the queue can be kept during daily business and will
         automatically resume during any free time
         
  con: - no complicated decision making on how subsequent tasks are executed,
         e.g., no possibility to do first a coarse ESR, then decide in which range
         to do a finer ESR, etc. 
       - it is easy to forget save jobs. If everything goes well this is not a problem,
         because the jobs can be saved later at any time, but if there is a crash,
         unsaved jobs are lost

2. using an independent MissionControl job that is not managed by the JobManager

  Write a new job, that is not managed by the JobManager, i.e., that runs independently
  of the queue. This Job will submit jobs to the queue as needed.
  
  pro: - allows complex ways to submit jobs, e.g., depending on the result of previous
         measurement, with analysis performed in between, etc.

  con: - cannot be changed after started
       - control job will often be 'new code' and thus may have errors. It is
         difficult to test --> error prone
"""

import numpy as np
import threading
import logging
import time

from tools.emod import Job

from tools.utility import GetSetItemsMixin

#ToDo: maybe introduce lock for 'state' variable on each job?

from traits.api import Array, File, Instance, Button
from traitsui.api import View, Item, HGroup, VGroup, InstanceEditor

from chaco.api import ArrayPlotData
from tools.chaco_addons import SavePlot as Plot, SaveTool
from enable.api import ComponentEditor

from measurements.odmr import ODMR

class Zeeman( Job ):#, GetSetItemsMixin ):

    """Zeeman measurement."""

    start_button = Button(label='start', desc='Start the measurement.')
    stop_button  = Button(label='stop', desc='Stop the measurement.')
    
    def _start_button_fired(self):
        """React to submit button. Submit the Job."""
        self.start() 

    def _stop_button_fired(self):
        """React to remove button. Remove the Job."""
        self.stop()

    current = Array(dtype=float)
    basename = File()

    odmr = Instance( ODMR, factory=ODMR )

    frequency = Array()

    line_data   = Instance( ArrayPlotData )
    line_plot   = Instance( Plot, editor=ComponentEditor() )

    traits_view = View(VGroup(HGroup(Item('start_button',   show_label=False),
                                     Item('stop_button',   show_label=False),
                                     Item('state', style='readonly'),
                                     Item('odmr', editor=InstanceEditor(), show_label=False),
                                     ),
                              VGroup(Item('basename'),
                                     #Item('current'), # ToDo: migrate to a custom TabularEditor
                                     ),
                              Item('line_plot', show_label=False, resizable=True),
                              ),
                       title='Zeeman', buttons=[], resizable=True
                       )

    get_set_items = ['current', 'frequency', 'odmr', '__doc__']
    

    def __init__(self, coil, **kwargs):
        self.coil=coil
        super(Zeeman,self).__init__(**kwargs)
        self._create_line_plot()
        self.on_trait_change(self._update_plot, 'frequency', dispatch='ui')

    def _run(self):

        try:
            self.state='run'
            
            if self.basename == '':
                raise ValueError('Filename missing. Please specify a filename and try again.')
            
            odmr = self.odmr
            if odmr.stop_time == np.inf:
                raise ValueError('ODMR stop time set to infinity.')
            #odmr.remove()
            delta_f = (odmr.frequency_end-odmr.frequency_begin)

            self.frequency = np.array(())

            for i,current_i in enumerate(self.current):
                
                self.coil.current = current_i
                # turn off the coil to wait until the drift due to heating is stable
                if current_i > 0.05:
                #    self.coil.current = 0.0
                    threading.currentThread().stop_request.wait(30.0)
                    if threading.currentThread().stop_request.isSet():
                        break
                #odmr.perform_fit=False
                odmr.submit()
                while odmr.state != 'done':
                    threading.currentThread().stop_request.wait(1.0)
                    if threading.currentThread().stop_request.isSet():
                        odmr.remove()
                        break
                if threading.currentThread().stop_request.isSet():
                    break
                time.sleep(2) # make sure that the job is also done from the point of view of the JobManager
                #odmr.perform_fit=True
                basename = self.basename
                try:
                    appendix = basename[-4:]
                    if appendix in ['.pyd','.pys','.asc','.txt']:
                        basename = basename[:-4]
                    else:
                        appendix = '.pys'
                except:
                    appendix = '.pys'
                filename = basename+'_'+str(current_i*1000)+'mA'+appendix
                filename1 = basename+'_'+str(current_i*1000)+'mA.png'
                #filename = basename+'_'+str(current_i)+'A'+appendix
                #filename1 = basename+'_'+str(current_i)+'A.png'
                odmr.save(filename)
                odmr.save_line_plot(filename1)
                # turn off the coil and wait for some time to prevent overheating
                #if current_i > 0.150:
                #    self.coil.current = 0.0
                #    threading.currentThread().stop_request.wait(60.0)
                #    if threading.currentThread().stop_request.isSet():
                #        break
                

                #f = odmr.fit_frequencies[0]
                #self.frequency=np.append(self.frequency,f)
                
                #odmr.frequency_begin = max(10e6,f-0.5*delta_f)
                #odmr.frequency_end = odmr.frequency_begin + delta_f 
                
            self.state='done'
        except:
            logging.getLogger().exception('Error in Zeeman.')
            self.state = 'error'

    def _create_line_plot(self):
        line_data = ArrayPlotData(current=np.array(()), frequency=np.array(()),)
        plot = Plot(line_data, padding=8, padding_left=64, padding_bottom=36)
        plot.plot(('current','frequency'), color='blue', name='zeeman')
        plot.index_axis.title = 'current [mA]'
        plot.value_axis.title = 'frequency [MHz]'
        plot.tools.append(SaveTool(plot))
        self.line_data = line_data
        self.line_plot = plot

    def _update_plot(self,new):
        n = len(new)
        self.line_data.set_data('current',self.current[:n])
        self.line_data.set_data('frequency',new*1e-6)

    def save_line_plot(self, filename):
        self.save_figure(self.line_plot, filename)

if __name__=='__main__':

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    from hardware.mocking import Sweeper, Microwave, Coil
    microwave = Microwave()
    sweeper = Sweeper()
    coil = Coil()

    from tools.emod import JobManager
    JobManager().start()

    from measurements.odmr import ODMR
    
    odmr = ODMR(microwave, sweeper)

    zeeman = Zeeman(coil, odmr=odmr, current=np.arange(0.052, 0.149, 0.002))
    zeeman.edit_traits()
    
    #zeeman.start()

