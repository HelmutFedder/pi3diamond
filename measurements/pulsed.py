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

from traits.api import Range, Int, Float, Bool, Array, Instance, Enum, on_trait_change, Button
from traitsui.api import View, Item, Tabbed, HGroup, VGroup, VSplit, EnumEditor, TextEditor

import logging
import time

from tools.emod import ManagedJob

from tools.utility import GetSetItemsMixin

#from TimeTagger import TimeDifferences
from hardware.mocking import TimeDifferences

"""
Several options to decide when to start  and when to restart a job, i.e. when to clear data, etc.

1. set a 'new' flag on every submit button

pro: simple, need not to think about anything in subclass

con: continue of measurement only possible by hack (manual submit to JobManager without submit button)
     submit button does not do what it says
     
2. check at start time whether this is a new measurement.

pro: 

con: complicated checking needed
     checking has to be reimplemented on sub classes
     no explicit way to restart the same measurement

3. provide user settable clear / keep flag

pro: explicit

con: user can forget

4. provide two different submit buttons: submit, resubmit

pro: explicit

con: two buttons that user may not understand
     user may use wrong button
     wrong button can result in errors

"""

from analysis.fitting import find_edge

# utility functions
def find_detection_pulses(sequence):
    """
    Find the number of detection triggers in a pulse sequence.
    """
    n = 0
    prev = []
    for channels, t in sequence:
        if 'detect' in channels and not 'detect' in prev:
            n+=1
        prev = channels
        if 'sequence' in channels:
            break
    return n

def sequence_length(sequence):
    """
    Return the total length of a pulse sequence.
    """
    t = 0
    for c,ti in sequence:
        t += ti
    return t

def sequence_union(s1, s2):
    """
    Return the union of two pulse sequences s1 and s2.
    """
    # make sure that s1 is the longer sequence and s2 is merged into it
    if sequence_length(s1) < sequence_length(s2):
        sp = s2
        s2 = s1
        s1 = sp
    s = []
    c1, dt1 = s1.pop(0)
    c2, dt2 = s2.pop(0)
    while True:
        if dt1 < dt2:
            s.append( ( set(c1) | set(c2), dt1) )
            dt2 -= dt1
            try:
                c1, dt1 = s1.pop(0)
            except:
                break
        elif dt2 < dt1:
            s.append( ( set(c1) | set(c2), dt2) )
            dt1 -= dt2
            try:
                c2, dt2 = s2.pop(0)
            except:
                c2 = []
                dt2 = np.inf
        else:
            s.append( ( set(c1) | set(c2), dt1) )
            try:
                c1, dt1 = s1.pop(0)
            except:
                break
            try:
                c2, dt2 = s2.pop(0)
            except:
                c2 = []
                dt2 = np.inf            
    return s

def sequence_simplify(sequence):
    """Merge adjacent pulses that have equal channels"""
    i = 0
    end = len(sequence)
    current = sequence[i]
    while i+1 < end:
        next = sequence[i+1]
        if current[0] == next[0]: # merge next into the current one, pop next, and decrease length
            sequence[i] = (current[0], current[1]+next[1])
            sequence.pop(i+1)
            end -= 1
        else: # move one to the right
            i += 1
            current = next
    return sequence

def sequence_remove_zeros(sequence):
    """remove all pulses with zero length from a sequence"""
    return filter(lambda x: x[1]!=0.0, sequence)

def spin_state(c, dt, T, t0=0.0, t1=-1.):
    
    """
    Compute the spin state from a 2D array of count data.
    
    Parameters:
    
        c    = count data
        dt   = time step
        t0   = beginning of integration window relative to the edge
        t1   = None or beginning of integration window for normalization relative to edge
        T    = width of integration window
        
    Returns:
    
        y       = 1D array that contains the spin state
        profile = 1D array that contains the pulse profile
        edge    = position of the edge that was found from the pulse profile
        
    If t1<0, no normalization is performed. If t1>=0, each data point is divided by
    the value from the second integration window and multiplied with the mean of
    all normalization windows.
    """

    profile = c.sum(0)
    edge = find_edge(profile)
    
    I = int(round(T/float(dt)))
    i0 = edge + int(round(t0/float(dt)))
    y = np.empty((c.shape[0],))
    for i, slot in enumerate(c):
        y[i] = slot[i0:i0+I].sum()
    if t1 >= 0:
        i1 = edge + int(round(t1/float(dt)))    
        y1 = np.empty((c.shape[0],))
        for i, slot in enumerate(c):
            y1[i] = slot[i1:i1+I].sum()
        y = y/y1*y1.mean()
    return y, profile, edge

class Pulsed( ManagedJob, GetSetItemsMixin ):
    
    """
    Defines a pulsed spin resonance measurement.
    
    We output a pulse sequence (consisting of laser
    and microwave control pulses) with some sort
    of pulse pattern generator.
    
    We acquire the pulsed fluorescence repsonsed in
    a two dimensional histogram. Each row in the histogram
    represents the fluorescence response corresponding
    to a laser pulse. The sum over a window of the fluorescence
    response encodes the spin state.
    """
    
    run_time = Float(value=0.0, label='run time [s]',format_str='%.f')
    stop_time = Float(default_value=np.inf, desc='Time after which the experiment stops by itself [s]', label='Stop time [s]', mode='text', auto_set=False, enter_set=True)
    stop_counts = Float(default_value=np.inf, desc='Stop the measurement when all data points of the extracted spin state have at least this many counts.', label='Stop counts', mode='text', auto_set=False, enter_set=True)

    keep_data = Bool(False) # helper variable to decide whether to keep existing data

    resubmit_button = Button(label='resubmit', desc='Submits the measurement to the job manager. Tries to keep previously acquired data. Behaves like a normal submit if sequence or time bins have changed since previous run.')
    
    # acquisition parameters    
    record_length = Float(default_value=3000, desc='length of acquisition record [ns]', label='record length [ns]', mode='text', auto_set=False, enter_set=True)
    bin_width = Range(low=0.1, high=1000., value=1.0, desc='bin width [ns]', label='bin width [ns]', mode='text', auto_set=False, enter_set=True)

    n_laser = Int(2)
    n_bins = Int(2)
    time_bins = Array(value=np.array((0,1)))
    
    sequence = Instance( list, factory=list )

    # measured data
    count_data = Array( value=np.zeros((2,2)) )

    # parameters for calculating spin state
    integration_width   = Float(default_value=200.,   desc='width of integration window [ns]',                     label='integr. width [ns]', mode='text', auto_set=False, enter_set=True)
    position_signal     = Float(default_value=0.,     desc='position of signal window relative to edge [ns]',         label='pos. signal [ns]', mode='text',   auto_set=False, enter_set=True)
    position_normalize  = Float(default_value=-1.,    desc='position of normalization window relative to edge [ns]. If negative, no normalization is performed',  label='pos. norm. [ns]', mode='text',    auto_set=False, enter_set=True)
    
    # analyzed data
    pulse               = Array(value=np.array((0.,0.)))
    edge                = Float(value=0.0)
    spin_state          = Array(value=np.array((0.,0.)))
    
    channel_apd = Int(0)
    channel_detect = Int(2)
    channel_sequence = Int(3)

    def __init__(self, pulse_generator, time_tagger, **kwargs):
        super(Pulsed, self).__init__(**kwargs)
        self.pulse_generator = pulse_generator
        self.time_tagger = time_tagger
    
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

    def generate_sequence(self):
        return []

    def apply_parameters(self):
        """Apply the current parameters and decide whether to keep previous data."""
        n_bins = int(self.record_length / self.bin_width)
        time_bins = self.bin_width*np.arange(n_bins)
        sequence = self.generate_sequence()
        n_laser = find_detection_pulses(sequence)

        if self.keep_data and sequence == self.sequence and np.all(time_bins == self.time_bins): # if the sequence and time_bins are the same as previous, keep existing data
            self.old_count_data = self.count_data.copy()
        else:
            self.old_count_data = np.zeros((n_laser,n_bins))
            self.count_data = np.zeros((n_laser,n_bins))
            self.run_time = 0.0
        
        self.sequence = sequence 
        self.time_bins = time_bins
        self.n_bins = n_bins
        self.n_laser = n_laser
        self.keep_data = True # when job manager stops and starts the job, data should be kept. Only new submission should clear data.

    def start_up(self):
        """Put here additional stuff to be executed at startup."""
        pass

    def shut_down(self):
        """Put here additional stuff to be executed at shut_down."""
        pass

    def _run(self):
        """Acquire data."""

        try: # try to run the acquisition from start_up to shut_down
            self.state='run'
            self.apply_parameters()

            if not (self.run_time < self.stop_time and any(self.spin_state<self.stop_counts)):
                logging.getLogger().debug('Runtime larger than stop_time. Returning')
                self.state='done'
                return

            self.start_up()
            self.pulse_generator.Night()
            pulsed = TimeDifferences(self.time_tagger,
                                     self.channel_apd,
                                     self.channel_detect, 
                                     self.channel_detect,
                                     self.channel_sequence,
                                     int(np.round(self.bin_width*1000)),
                                     self.n_bins,
                                     self.n_laser,
                                     )
            self.pulse_generator.Sequence(self.sequence)
            self.pulse_generator.checkUnderflow()
            
            while self.run_time < self.stop_time and any(self.spin_state<self.stop_counts):
                start_time = time.time()
                self.thread.stop_request.wait(1.0)
                if self.thread.stop_request.isSet():
                    logging.getLogger().debug('Caught stop signal. Exiting.')
                    break
                if self.pulse_generator.checkUnderflow():
                    raise RuntimeError('Underflow in pulse generator.')
                self.count_data = self.old_count_data + pulsed.getData()
                self.run_time += time.time() - start_time

            if self.run_time < self.stop_time and any(self.spin_state<self.stop_counts):
                self.state = 'idle'
            else:
                try:
                    self.save(self.filename)
                except:
                    logging.getLogger().exception('Failed to save the data to file.')
                self.state='done'
            del pulsed
            self.shut_down()
            self.pulse_generator.Light()

        except: # if anything fails, log the exception and set the state
            logging.getLogger().exception('Something went wrong in pulsed loop.')
            self.state='error'

    @on_trait_change('count_data,integration_width,position_signal,position_normalize')
    def _compute_spin_state(self):
        y, profile, edge = spin_state(c=self.count_data,
                                      dt=self.bin_width,
                                      T=self.integration_width,
                                      t0=self.position_signal,
                                      t1=self.position_normalize,
                                      )
        self.spin_state = y
        self.pulse = profile
        self.edge = self.time_bins[edge]

    traits_view = View(VGroup(HGroup(Item('submit_button',   show_label=False),
                                     Item('remove_button',   show_label=False),
                                     Item('resubmit_button', show_label=False),
                                     Item('priority'),
                                     Item('state', style='readonly'),
                                     Item('run_time', style='readonly', format_str='%.f'),
                                     Item('stop_time', format_str='%.f'),
                                     Item('stop_counts'),
                                     ),
                              HGroup(Item('filename',springy=True),
                                     Item('save_button', show_label=False),
                                     Item('load_button', show_label=False)
                                     ),
                              HGroup(Item('bin_width', enabled_when='state != "run"'),
                                     Item('record_length', enabled_when='state != "run"'),
                                     ),
                              ),
                       title='Pulsed Measurement',
                       )

    get_set_items = ['__doc__','record_length','bin_width','n_bins','time_bins','n_laser','sequence','count_data','run_time',
                     'integration_width','position_signal','position_normalize',
                     'pulse','edge','spin_state']


class PulsedTau( Pulsed ):

    """Defines a Pulsed measurement with tau mesh."""

    tau_begin   = Float(default_value=0.,     desc='tau begin [ns]',  label='tau begin [ns]',   mode='text', auto_set=False, enter_set=True)
    tau_end     = Float(default_value=300.,   desc='tau end [ns]',    label='tau end [ns]',     mode='text', auto_set=False, enter_set=True)
    tau_delta   = Float(default_value=3.,      desc='delta tau [ns]',  label='delta tau [ns]',   mode='text', auto_set=False, enter_set=True)

    tau = Array( value=np.array((0.,1.)) )

    def apply_parameters(self):
        """Overwrites apply_parameters() from pulsed. Prior to generating sequence, etc., generate the tau mesh."""
        self.tau = np.arange(self.tau_begin, self.tau_end, self.tau_delta)
        Pulsed.apply_parameters(self)

    get_set_items = Pulsed.get_set_items + ['tau_begin','tau_end','tau_delta','tau']

    traits_view = View(VGroup(HGroup(Item('submit_button',   show_label=False),
                                     Item('remove_button',   show_label=False),
                                     Item('resubmit_button', show_label=False),
                                     Item('priority'),
                                     Item('state', style='readonly'),
                                     Item('run_time', style='readonly',format_str='%.f'),
                                     Item('stop_time'),
                                     ),
                              HGroup(Item('filename',springy=True),
                                     Item('save_button', show_label=False),
                                     Item('load_button', show_label=False)
                                     ),
                              HGroup(Item('bin_width', enabled_when='state != "run"'),
                                     Item('record_length', enabled_when='state != "run"'),
                                     ),
                              ),
                       title='PulsedTau Measurement',
                       )

"""
This module provides analysis of pulsed measurements.

The first part provides simple numeric functions.

The second part provides Trait GUIs
"""

import numpy as np

#########################################
# Trait GUIs for pulsed fits
#########################################

from traits.api import Instance, Any, Property, Range, Float, Int, Bool, Array, List, Str, Tuple, Enum,\
                                 on_trait_change, cached_property, DelegatesTo
from traitsui.api import View, Item, Tabbed, Group, HGroup, VGroup, VSplit, EnumEditor, TextEditor, InstanceEditor
from enable.api import ComponentEditor
from chaco.api import ArrayDataSource, LinePlot, LinearMapper, ArrayPlotData, Spectral, PlotLabel

from tools.chaco_addons import SavePlot as Plot, SaveTool

import threading
import time
import logging

from tools.utility import GetSetItemsMixin

class PulsedTool( GetSetItemsMixin ):

    """
    Widget to plot a pulsed measurement.
    
    We plot
      * the raw data as an image
      * the average fluorescence response as a line plot
        (the average fluorescence response
        is used to determine a trigger point).
      * the spin state as a line plot
    """

    # the measurement to analyze
    measurement = Any(editor=InstanceEditor)
    
    # plotting
    matrix_data = Instance( ArrayPlotData)
    line_data   = Instance( ArrayPlotData )
    pulse_data  = Instance( ArrayPlotData )
    matrix_plot = Instance( Plot, editor=ComponentEditor() )
    pulse_plot  = Instance( Plot, editor=ComponentEditor() )
    line_plot   = Instance( Plot, editor=ComponentEditor() )

    get_set_items = ['__doc__','measurement']

    traits_view = View(VGroup(Item(name='measurement', style='custom', show_label=False),
                              VSplit(Item('matrix_plot', show_label=False, width=500, height=-300, resizable=True),
                                     Item('line_plot', show_label=False, width=500, height=-300, resizable=True),
                                     Item('pulse_plot', show_label=False, width=500, height=-300, resizable=True),
                                     ),
                              ),
                       title='Pulsed Tool',
                       buttons=[],
                       resizable=True,
                       height=-640
                       )

    def __init__(self, **kwargs):
        super(PulsedTool, self).__init__(**kwargs)
        self._create_matrix_plot()
        self._create_pulse_plot()
        self._create_line_plot()
        self.on_trait_change(self._update_matrix_index, 'measurement.time_bins,measurement.n_laser', dispatch='ui')
        self.on_trait_change(self._update_matrix_value, 'measurement.count_data', dispatch='ui')
        self.on_trait_change(self._update_pulse_index, 'measurement.time_bins', dispatch='ui')
        self.on_trait_change(self._update_pulse_value, 'measurement.pulse', dispatch='ui')
        self.on_trait_change(self._update_line_plot_value, 'measurement.spin_state', dispatch='ui')
        self.on_trait_change(self._on_edge_change, 'measurement.edge', dispatch='ui')

    # plotting
    def _create_matrix_plot(self):
        matrix_data = ArrayPlotData(image=np.zeros((2,2)))
        plot = Plot(matrix_data, width=500, height=500, resizable='hv', padding=8, padding_left=48, padding_bottom=36)
        plot.index_axis.title = 'time [ns]'
        plot.value_axis.title = 'laser pulse #'
        plot.img_plot('image',
                      xbounds=(0,1),
                      ybounds=(0,1),
                      colormap=Spectral)[0]
        plot.tools.append(SaveTool(plot))
        self.matrix_data = matrix_data
        self.matrix_plot = plot
    
    def _create_pulse_plot(self):
        pulse_data  = ArrayPlotData(x=np.array((0.,0.1,0.2)),y=np.array((0,1,2)))
        plot = Plot(pulse_data, padding=8, padding_left=64, padding_bottom=36)    
        line = plot.plot(('x','y'), style='line', color='blue', name='data')[0]
        plot.index_axis.title = 'time [ns]'
        plot.value_axis.title = 'intensity'
        edge_marker = LinePlot(index = ArrayDataSource(np.array((0,0))),
                               value = ArrayDataSource(np.array((0,1e9))),
                               color = 'red',
                               index_mapper = LinearMapper(range=plot.index_range),
                               value_mapper = LinearMapper(range=plot.value_range),
                               name='marker')
        plot.add(edge_marker)
        plot.tools.append(SaveTool(plot))
        self.pulse_data = pulse_data
        self.pulse_plot = plot
        
    def _create_line_plot(self):
        line_data   = ArrayPlotData(index=np.array((0,1)), spin_state=np.array((0,0)),)
        plot = Plot(line_data, padding=8, padding_left=64, padding_bottom=36)
        plot.plot(('index','spin_state'), color='blue', name='spin_state')
        plot.index_axis.title = 'pulse #'
        plot.value_axis.title = 'spin state'
        plot.tools.append(SaveTool(plot))
        self.line_data = line_data
        self.line_plot = plot

    def _update_matrix_index(self):
        if self.measurement is not None:
            self.matrix_plot.components[0].index.set_data((self.measurement.time_bins[0], self.measurement.time_bins[-1]),(0.0,float(self.measurement.n_laser)))
        
    def _update_matrix_value(self):
        if self.measurement is not None:
            s = self.measurement.count_data.shape
            if not s[0]*s[1] > 1000000:
                self.matrix_data.set_data('image', self.measurement.count_data)
    def _update_pulse_index(self):
        if self.measurement is not None:
            self.pulse_data.set_data('x', self.measurement.time_bins)        
    def _update_pulse_value(self):
        if self.measurement is not None:
            self.pulse_data.set_data('y', self.measurement.pulse)
    def _on_edge_change(self):
        if self.measurement is not None:
            y=self.measurement.edge
            self.pulse_plot.components[1].index.set_data(np.array((y,y)))
    def _update_line_plot_value(self):
        if self.measurement is not None:
            y=self.measurement.spin_state
            n = len(y)
            old_index = self.line_data.get_data('index')
            if old_index is not None and len(old_index) != n:
                self.line_data.set_data('index',np.arange(n))
            self.line_data.set_data('spin_state',y)

    def save_matrix_plot(self, filename):
        save_figure(self.matrix_plot, filename)
    
    def save_line_plot(self, filename):
        save_figure(self.line_plot, filename)


class PulsedToolTau( PulsedTool ):

    """
    Analysis of a pulsed measurement with a 'tau' as index-data.
    """

    # overwrite __init__ such that change of 'tau' causes plot update 
    def __init__(self, **kwargs):
        super(PulsedToolTau, self).__init__(**kwargs)
        self.on_trait_change(self._on_tau_change, 'measurement.tau', dispatch='ui')
        #self.on_trait_change(self._on_tau_change, 'measurement.tau') # ToDo: fix this

    # overwrite the line_plot such that the x-axis label is time 
    def _create_line_plot(self):
        line_data   = ArrayPlotData(index=np.array((0,1)), spin_state=np.array((0,0)),)
        plot = Plot(line_data, padding=8, padding_left=64, padding_bottom=36)
        plot.plot(('index','spin_state'), color='blue', name='spin_state')
        plot.index_axis.title = 'time [micro s]'
        plot.value_axis.title = 'spin state'
        plot.tools.append(SaveTool(plot))
        self.line_data = line_data
        self.line_plot = plot

    # overwrite this one to throw out setting of index data according to length of spin_state
    def _update_line_plot_value(self):
        y = self.measurement.spin_state
        self.line_data.set_data('spin_state',y)

    # provide method for update of tau
    def _on_tau_change(self):
        self.line_data.set_data('index',self.measurement.tau*1e-3)

    # overwrite this to change the window title
    traits_view = View(VGroup(Item(name='measurement', style='custom', show_label=False),
                              VSplit(Item('matrix_plot', show_label=False, width=500, height=300, resizable=True),
                                     Item('line_plot', show_label=False, width=500, height=300, resizable=True),
                                     Item('pulse_plot', show_label=False, width=500, height=300, resizable=True),
                                     ),
                              ),
                       title='Pulsed Analysis Tau',
                       buttons=[], resizable=True
                       )


class FitToolTau( PulsedToolTau ):

    """
    Base class for PulsedTool with a tau and fit.
    """

    # fit results
    fit_result = Tuple()
    label_text = Str('')

    # add fit results to the get_set_items
    get_set_items = PulsedToolTau.get_set_items + ['fit_result','label_text']

    # overwrite __init__ to trigger update events
    def __init__(self, **kwargs):
        super(FitToolTau, self).__init__(**kwargs)
        self.on_trait_change(self._update_fit, 'measurement.spin_state', dispatch='ui')
        self.on_trait_change(self._on_fit_result_change, 'fit_result', dispatch='ui')
        self.on_trait_change(self._on_label_text_change, 'label_text', dispatch='ui')
    
    def _update_fit(self):
        pass
        
    # overwrite the line_plot to include fit and text label 
    def _create_line_plot(self):
        line_data   = ArrayPlotData(index=np.array((0,1)),
                                    spin_state=np.array((0,0)),
                                    fit=np.array((0,0)))
        plot = Plot(line_data, padding=8, padding_left=64, padding_bottom=36)
        plot.plot(('index','spin_state'), color='blue', name='spin_state')
        plot.plot(('index','fit'), color='red', name='fit')
        plot.index_axis.title = 'time [micro s]'
        plot.value_axis.title = 'spin state'
        plot.overlays.insert(0, PlotLabel(text=self.label_text, hjustify='left', vjustify='bottom', position=[64,32]) )
        plot.tools.append(SaveTool(plot))
        self.line_data = line_data
        self.line_plot = plot

    def _on_fit_result_change(self, new):
        pass
    
    def _on_label_text_change(self, new):
        self.line_plot.overlays[0].text = new    
    
                       

if __name__ == '__main__':
    
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')
    
    from tools.emod import JobManager
    
    JobManager().start()

    from hardware.mocking import PulseGenerator
    pulse_generator = PulseGenerator()

    #from TimeTagger import createMockingTagger
    #tagger = createMockingTagger()

    from mocking import createMockingTagger
    tagger = createMockingTagger()

    pulsed = PulsedTau(pulse_generator, tagger)
    pulsed.edit_traits()
    
