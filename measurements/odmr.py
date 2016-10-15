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

import numpy as np

from traits.api import Trait, Instance, Property, String, Range, Float, Int, Str, Bool, Array, Enum, Button
from traitsui.api import View, Item, HGroup, VGroup, VSplit, Tabbed, EnumEditor, TextEditor, Group
from enable.api import Component, ComponentEditor
from chaco.api import ArrayPlotData, Plot, Spectral, PlotLabel

from tools.chaco_addons import SavePlot as Plot, SaveTool

import time
import threading
import logging

from tools.emod import ManagedJob

from analysis.fitting import find_peaks, Lorentzian

from tools.utility import GetSetItemsMixin

class ODMR( ManagedJob, GetSetItemsMixin ):
    """
    Implements an Optically Detected Magnetic Resonance (ODMR) measurement.
    
    Here we sweep a microwave source and record
    the photon clicks in every point of the sweep.
    
    This measurement requires a microwwave source
    and a counter. The counter also generates
    a trigger for every counting bin that steps
    the microwave source to the next frequency.
    
    The results from successive sweeps are accumulated.
    
    Optionally one can run a pulsed ODMR sweep with
    microwave pi-pulses.
    
    We provide some basic fitting.
    """

    # starting and stopping
    keep_data = Bool(False) # helper variable to decide whether to keep existing data
    resubmit_button = Button(label='resubmit', desc='Submits the measurement to the job manager. Tries to keep previously acquired data. Behaves like a normal submit if sequence or time bins have changed since previous run.')    
    
    # measurement parameters
    power = Range(low=-100., high=25., value=-20, desc='Power [dBm]', label='Power [dBm]', mode='text', auto_set=False, enter_set=True)
    frequency_begin = Float(default_value=2.85e9,    desc='Start Frequency [Hz]',    label='Begin [Hz]', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%e'))
    frequency_end   = Float(default_value=2.88e9,    desc='Stop Frequency [Hz]',     label='End [Hz]', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%e'))
    frequency_delta = Float(default_value=1e6,       desc='frequency step [Hz]',     label='Delta [Hz]', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%e'))
    t_pi  = Range(low=1., high=100000., value=1000., desc='length of pi pulse [ns]', label='pi [ns]', mode='text', auto_set=False, enter_set=True)
    laser = Range(low=1., high=10000., value=300., desc='laser [ns]', label='laser [ns]', mode='text', auto_set=False, enter_set=True)
    wait  = Range(low=1., high=10000., value=1000., desc='wait [ns]', label='wait [ns]', mode='text', auto_set=False, enter_set=True)
    pulsed = Bool(False, label='pulsed')
    seconds_per_point = Range(low=20e-3, high=1, value=20e-3, desc='Seconds per point', label='Sec per point', mode='text', auto_set=False, enter_set=True)
    stop_time = Range(low=1., value=np.inf, desc='Time after which the experiment stops by itself [s]', label='Stop time [s]', mode='text', auto_set=False, enter_set=True)
    n_lines = Range (low=1, high=10000, value=50, desc='Number of lines in Matrix', label='Matrix lines', mode='text', auto_set=False, enter_set=True)
    
    # control data fitting
    fit_method = Enum('Lorentzian', 'N14 Lorentzian', desc='Fit Method',     label='Fit Method')
    fit = Bool(False, label='fit')
    number_of_resonances = Int(1, desc='Number of peaks to be fitted', label='num peaks', mode='text', auto_set=False, enter_set=True)
    width = Float(5e6, desc='width used for peak detection.', label='width [Hz]', mode='text', auto_set=False, enter_set=True)
    
    # fit result
    fit_parameters = Array(value=np.array(()))
    fit_frequencies = Array(value=np.array(()), label='frequency [Hz]') 
    fit_line_width = Array(value=np.array(()), label='line_width [Hz]') 
    fit_contrast = Array(value=np.array(()), label='contrast [%]')
    fit_peak_ind = Array(value=np.array(()), label='peak pos [Hz]')
    fit_peak_val = Array(value=np.array(()), label='peak val')
    fit_string = Str()
    
    # measurement data    
    frequency = Array()
    counts = Array()
    counts_matrix = Array()
    run_time = Float(value=0.0, desc='Run time [s]', label='Run time [s]')

    # plotting
    line_label  = Instance( PlotLabel )
    line_data   = Instance( ArrayPlotData )
    matrix_data = Instance( ArrayPlotData )
    line_plot   = Instance( Plot, editor=ComponentEditor() )
    matrix_plot = Instance( Plot, editor=ComponentEditor() )

    get_set_items = ['frequency', 'counts', 'counts_matrix',
                     'number_of_resonances', 'width',
                     'fit_parameters', 'fit_contrast', 'fit_line_width', 'fit_frequencies',
                     'fit_peak_ind', 'fit_peak_val', 'fit_string',
                     'fit', 'run_time',
                     'power', 'frequency_begin', 'frequency_end', 'frequency_delta',
                     'laser', 'wait', 'pulsed', 't_pi',
                     'seconds_per_point', 'stop_time', 'n_lines',
                     '__doc__']

    def __init__(self, microwave, counter, pulse_generator=None, **kwargs):
        super(ODMR, self).__init__(**kwargs)
        self.microwave = microwave
        self.counter = counter
        self.pulse_generator = pulse_generator
        self._create_line_plot()
        self._create_matrix_plot()
        self.on_trait_change(self._update_line_data_index,      'frequency',            dispatch='ui')
        self.on_trait_change(self._update_line_data_value,      'counts',               dispatch='ui')
        self.on_trait_change(self._update_matrix_data_value,    'counts_matrix',        dispatch='ui')
        self.on_trait_change(self._update_matrix_data_index,    'n_lines,frequency',    dispatch='ui')
        self.on_trait_change(self._update_fit, 'counts,fit,number_of_resonances,width,fit_method', dispatch='ui')

    def _counts_matrix_default(self):
        return np.zeros( (self.n_lines, len(self.frequency)) )

    def _frequency_default(self):
        return np.arange(self.frequency_begin, self.frequency_end+self.frequency_delta, self.frequency_delta)

    def _counts_default(self):
        return np.zeros(self.frequency.shape)

    # data acquisition
    def apply_parameters(self):
        """Apply the current parameters and decide whether to keep previous data."""
        frequency = np.arange(self.frequency_begin, self.frequency_end+self.frequency_delta, self.frequency_delta)

        if not self.keep_data or np.any(frequency != self.frequency):
            self.frequency = frequency
            self.counts = np.zeros(frequency.shape)
            self.run_time = 0.0

        self.keep_data = True # when job manager stops and starts the job, data should be kept. Only new submission should clear data.

    def _run(self):
                
        try:
            self.state='run'
            self.apply_parameters()

            if self.run_time >= self.stop_time:
                self.state='done'
                return

            # if pulsed, turn on sequence
            if self.pulse_generator:
                if self.pulsed:
                    self.pulse_generator.Sequence( 100 * [ (['detect','aom'],self.laser), ([],self.wait), (['microwave'],self.t_pi) ] )
                else:
                    self.pulse_generator.Open()
            else:
                if self.pulsed:
                    raise ValueError("pulse_generator not defined while running measurement in pulsed mode.")

            n = len(self.frequency)

            """
            self.microwave.setOutput( self.power, np.append(self.frequency,self.frequency[0]), self.seconds_per_point)
            self._prepareCounter(n)
            """
            self.microwave.setPower(self.power)
            self.microwave.initSweep( self.frequency, self.power*np.ones(self.frequency.shape))
            self.counter.configure(n, self.seconds_per_point, DutyCycle=0.8)
            time.sleep(0.5)

            while self.run_time < self.stop_time:
                start_time = time.time()
                if threading.currentThread().stop_request.isSet():
                    break
                self.microwave.resetListPos()
                counts = self.counter.run()                
                self.run_time += time.time() - start_time
                self.counts += counts
                self.counts_matrix = np.vstack( (counts, self.counts_matrix[:-1,:]) )
                self.trait_property_changed('counts', self.counts)
                """
                self.microwave.doSweep()
                
                timeout = 3.
                start_time = time.time()
                while not self._count_between_markers.ready():
                    time.sleep(0.1)
                    if time.time() - start_time > timeout:
                        print "count between markers timeout in ODMR"
                        break
                        
                counts = self._count_between_markers.getData(0)
                self._count_between_markers.clean()
                """
    
            self.microwave.setOutput( None, self.frequency_begin)
            if self.pulse_generator:
                self.pulse_generator.Light()
            self.counter.clear()
        except:
            logging.getLogger().exception('Error in odmr.')
            self.microwave.setOutput( None, self.frequency_begin)
            self.state = 'error'
        else:
            if self.run_time < self.stop_time:
                self.state = 'idle'            
            else:
                try:
                    self.save(self.filename)
                except:
                    logging.getLogger().exception('Failed to save the data to file.')
                self.state='done'

    # fitting
    def _update_fit(self):
        if self.fit:
            try:
                if self.fit_method == 'Lorentzian':
                    fit_res = find_peaks(self.frequency,self.counts,self.width,self.number_of_resonances)
                elif self.fit_method == 'N14 Lorentzian':
                    raise NotImplementedError('you may find it in old code')
                else:
                    raise ValueError('unknown fit method')
            except Exception:
                logging.getLogger().debug('ODMR fit failed.', exc_info=True)
                fit_res = {}
                #p = np.nan*np.empty(4)
        else:
            fit_res = {}
        self.fit_res = fit_res
        if 'p' in fit_res:
            p = fit_res['p']
            x0 = fit_res['x0']
            y0 = fit_res['y0']
            self.fit_peak_ind = x0
            self.fit_peak_val = y0
            self.fit_parameters = p
            self.fit_frequencies = p[1::3]
            self.fit_line_width = p[2::3]

            n = len(p)/3
            contrast = np.empty(n)
            c = p[0]
            pp=p[1:].reshape((n,3))
            for i,pi in enumerate(pp):
                area = pi[2]
                hwhm = pi[1]
                amp = np.abs(area/(np.pi*hwhm))
                if area > 0:
                    contrast[i] = 100*amp/(amp+c)
                else:
                    contrast[i] = 100*amp/c
                    
            self.fit_contrast = contrast

            s = ''
            for i, fi in enumerate(self.fit_frequencies):
                s += 'f %i: %.6e Hz, HWHM %.3e Hz, contrast %.1f%%\n'%(i+1, fi, self.fit_line_width[i], self.fit_contrast[i])
                
            self.fit_string = s

        elif 'x0' in fit_res:
            x0 = fit_res['x0']
            y0 = fit_res['y0']
            self.fit_peak_ind = x0
            self.fit_peak_val = y0
            self.fit_parameters = np.array(())
            self.fit_frequencies = x0
            self.fit_line_width = np.array(())

            s = ''
            for i, fi in enumerate(self.fit_frequencies):
                s += 'f %i: %.6e Hz\n'%(i+1, fi)
                
            self.fit_string = s

        else:
            self.fit_peak_ind = np.array(())
            self.fit_peak_val = np.array(())
            self.fit_parameters = np.array(())
            self.fit_frequencies = np.array(())
            self.fit_line_width = np.array(())
            self.fit_string = ''

    def _fit_string_changed(self, new):
        self.line_label.text = new
    
    # plotting        
    def _create_line_plot(self):
        line_data = ArrayPlotData(frequency=np.array((0.,1.)), counts=np.array((0.,0.)), fit=np.array((0.,0.))) 
        line_plot = Plot(line_data, padding=8, padding_left=64, padding_bottom=32)
        line_plot.plot(('frequency','counts'), style='line', color='blue', name='data')
        line_plot.index_axis.title = 'Frequency [MHz]'
        line_plot.value_axis.title = 'Fluorescence counts'
        line_label = PlotLabel(text='', hjustify='left', vjustify='bottom', position=[64,128])
        line_plot.overlays.append(line_label)
        line_plot.tools.append(SaveTool(line_plot))
        self.line_label = line_label
        self.line_data = line_data
        self.line_plot = line_plot
        
    def _create_matrix_plot(self):
        matrix_data = ArrayPlotData(image=np.zeros((2,2)))
        matrix_plot = Plot(matrix_data, padding=8, padding_left=64, padding_bottom=32)
        matrix_plot.index_axis.title = 'Frequency [MHz]'
        matrix_plot.value_axis.title = 'line #'
        matrix_plot.img_plot('image',
                             xbounds=(self.frequency[0],self.frequency[-1]),
                             ybounds=(0,self.n_lines),
                             colormap=Spectral)
        matrix_plot.tools.append(SaveTool(matrix_plot))
        self.matrix_data = matrix_data
        self.matrix_plot = matrix_plot

    def _fit_parameters_changed(self, new):
        n = len(new)/3
        line_plot = self.line_plot
        line_data = self.line_data
        for name in line_plot.plots.keys():
            if name != 'data':
                line_plot.delplot(name)
        for i in range(n):
            p = np.append(new[0], new[i*3+1:i*3+4])
            name = 'peak_%i'%i
            line_data.set_data(name, Lorentzian(*p)(self.frequency))
            line_plot.plot(('frequency',name), style='line', color='red', name=name)
        line_plot.request_redraw()

    #def _fit_changed(self,new):
        #plot = self.line_plot
        #if new:
            #plot.plot(('frequency','fit'), style='line', color='red', name='fit')
            #self.line_label.visible=True
        #else:
            #plot.delplot('fit')
            #self.line_label.visible=False
        #plot.request_redraw()

    def _update_line_data_index(self):
        self.line_data.set_data('frequency', self.frequency*1e-6)
        self.counts_matrix = self._counts_matrix_default()

    def _update_line_data_value(self):
        self.line_data.set_data('counts', self.counts)

    def _update_matrix_data_value(self):
        self.matrix_data.set_data('image', self.counts_matrix)

    def _update_matrix_data_index(self):
        if self.n_lines > self.counts_matrix.shape[0]:
            self.counts_matrix = np.vstack( (self.counts_matrix,np.zeros((self.n_lines-self.counts_matrix.shape[0], self.counts_matrix.shape[1]))) )
        else:
            self.counts_matrix = self.counts_matrix[:self.n_lines]
        self.matrix_plot.components[0].index.set_data((self.frequency.min()*1e-6, self.frequency.max()*1e-6),(0.0,float(self.n_lines)))

    # react to GUI events
    def submit(self):
        """Submit the job to the JobManager."""
        self.keep_data = False
        ManagedJob.submit(self)

    def resubmit(self):
        """Submit the job to the JobManager."""
        self.keep_data = True
        ManagedJob.submit(self)

    def _resubmit_button_fired(self):
        """React to start button. Submit the Job."""
        self.resubmit() 

    traits_view = View(VGroup(HGroup(Item('submit_button',   show_label=False),
                                     Item('remove_button',   show_label=False),
                                     Item('resubmit_button', show_label=False),
                                     Item('priority', enabled_when='state != "run"'),
                                     Item('state', style='readonly'),
                                     Item('run_time', style='readonly',format_str='%.f'),
                                     Item('stop_time'),
                                     ),
                              HGroup(Item('filename',springy=True),
                                     Item('save_button', show_label=False),
                                     Item('load_button', show_label=False)
                                     ),
                              VGroup(HGroup(Item('power', enabled_when='state != "run"'),
                                            Item('frequency_begin', enabled_when='state != "run"'),
                                            Item('frequency_end', enabled_when='state != "run"'),
                                            Item('frequency_delta', enabled_when='state != "run"'),
                                            Item('seconds_per_point', enabled_when='state != "run"'),
                                            ),
                                     HGroup(Item('pulsed', enabled_when='state != "run"'),
                                            Item('laser', enabled_when='state != "run"'),
                                            Item('wait', enabled_when='state != "run"'),
                                            Item('t_pi', enabled_when='state != "run"'),
                                            ),
                                     HGroup(Item('fit'),
                                            Item('number_of_resonances'),
                                            Item('width'),
                                            Item('fit_method'),
                                            Item('n_lines'),
                                            ),
                                     #Item('fit_string', style='readonly'),
                                     #HGroup(Item('fit_contrast', style='readonly'),
                                     #       Item('fit_line_width', style='readonly'),
                                     #       Item('fit_frequencies', style='readonly'),
                                     #       ),
                                     ),
                              VSplit(Item('matrix_plot', show_label=False, resizable=True),
                                     Item('line_plot', show_label=False, resizable=True),
                                     ),
                              ),
                       title='ODMR', width=900, height=800, buttons=[], resizable=True
                       )


if __name__ == '__main__':

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    from hardware.mocking import Sweeper, Microwave
    microwave = Microwave()
    sweeper = Sweeper()

    from tools.emod import JobManager
    JobManager().start()

    o = ODMR(microwave, sweeper)
    o.edit_traits()
    
    
