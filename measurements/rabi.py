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

import logging

from pulsed import PulsedTau, sequence_remove_zeros, sequence_union, PulsedToolTau, FitToolTau

from analysis import fitting

from traits.api import Range, Tuple, Int, Float, Bool, Array, Instance, Enum, on_trait_change, Button
from traitsui.api import View, Item, Group, Tabbed, HGroup, VGroup, VSplit, EnumEditor, TextEditor

class Rabi( PulsedTau ):
    
    """Defines a Rabi measurement."""
    
    frequency  = Range(low=1,      high=20e9,  value=2.8705e9, desc='microwave frequency', label='frequency [Hz]', mode='text', auto_set=False, enter_set=True, editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%e'))
    power      = Range(low=-100.,  high=25.,   value=-20,      desc='microwave power',     label='power [dBm]',    mode='text', auto_set=False, enter_set=True)
    
    laser      = Float(default_value=3000., desc='laser [ns]',      label='laser [ns]',       mode='text', auto_set=False, enter_set=True)
    decay_init = Float(default_value=1000., desc='time to let the system decay after laser pulse [ns]',       label='decay init [ns]',        mode='text', auto_set=False, enter_set=True)
    decay_read = Float(default_value=0.,    desc='time to let the system decay before laser pulse [ns]',       label='decay read [ns]',        mode='text', auto_set=False, enter_set=True)
    
    aom_delay  = Float(default_value=0.,    desc='If set to a value other than 0.0, the aom triggers are applied\nearlier by the specified value. Use with care!', label='aom delay [ns]', mode='text', auto_set=False, enter_set=True)
    
    def __init__(self, pulse_generator, time_tagger, microwave, **kwargs):
        super(Rabi, self).__init__(pulse_generator, time_tagger, **kwargs)
        self.microwave = microwave
    
    def start_up(self):
        self.pulse_generator.Night()
        self.microwave.setOutput(self.power, self.frequency)

    def shut_down(self):
        self.pulse_generator.Light()
        self.microwave.setOutput(None, self.frequency)

    def generate_sequence(self):
        tau = self.tau
        laser = self.laser
        decay_init = self.decay_init
        decay_read = self.decay_read
        aom_delay = self.aom_delay
        if aom_delay == 0.0:
            sequence = [ (['aom'],laser) ]
            for t in tau:
                sequence += [ ([],decay_init), (['microwave'],t), ([],decay_read), (['detect','aom'],laser) ]
            sequence += [ (['sequence'], 100) ]
            sequence = sequence_remove_zeros(sequence)
        else:
            s1 = [ (['aom'],laser) ]
            s2 = [ ([],aom_delay+laser) ]
            for t in tau:
                s1 += [ ([], decay_init+t+decay_read), (['aom'], laser) ]
                s2 += [ ([], decay_init), (['microwave'],t), ([],decay_read), (['detect'],laser) ]
            s2 += [ (['sequence'],100) ]
            s1 = sequence_remove_zeros(s1)            
            s2 = sequence_remove_zeros(s2)            
            sequence = sequence_union(s1,s2)
        return sequence

    get_set_items = PulsedTau.get_set_items + ['frequency','power','laser','decay_init','decay_read','aom_delay']

    traits_view = View(VGroup(HGroup(Item('submit_button', show_label=False),
                                     Item('remove_button', show_label=False),
                                     Item('resubmit_button', show_label=False),
                                     Item('priority'),
                                     Item('state', style='readonly'),
                                     Item('run_time', style='readonly',format_str='%.f'),
                                     Item('stop_time',format_str='%.f'),
                                     Item('stop_counts'),
                                     ),
                              HGroup(Item('filename',springy=True),
                                     Item('save_button', show_label=False),
                                     Item('load_button', show_label=False)
                                     ),
                              Tabbed(VGroup(HGroup(Item('frequency', enabled_when='state != "run"'),
                                                   Item('power', enabled_when='state != "run"'),
                                                   Item('aom_delay', enabled_when='state != "run"'),
                                                   ),
                                            HGroup(Item('laser', enabled_when='state != "run"'),
                                                   Item('decay_init', enabled_when='state != "run"'),
                                                   Item('decay_read', enabled_when='state != "run"'),
                                                   ),
                                            HGroup(Item('tau_begin', enabled_when='state != "run"'),
                                                   Item('tau_end', enabled_when='state != "run"'),
                                                   Item('tau_delta', enabled_when='state != "run"'),
                                                   ),
                                            label='stimulation'),
                                     VGroup(HGroup(Item('record_length', enabled_when='state != "run"'),
                                                   Item('bin_width', enabled_when='state != "run"'),
                                                   ),
                                            label='acquisition'),
                                     VGroup(HGroup(Item('integration_width'),
                                                   Item('position_signal'),
                                                   Item('position_normalize'),
                                                   ),
                                            label='analysis'),
                                     ),
                              ),
                       title='Rabi',
                       )


class RabiTool( FitToolTau ):
    
    """
    Analysis of a Rabi measurement.
    """
    
    measurement = Instance( Rabi )
    
    # fit results
    contrast    = Tuple((np.nan,np.nan), editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str=' %.1f+-%.1f %%'))
    period      = Tuple((np.nan,np.nan), editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str=' %.2f+-%.2f'))
    q           = Float(np.nan, editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_func=lambda x:(' %.3f' if x >= 0.001 else ' %.2e')%x))
    t_pi2       = Tuple((np.nan,np.nan), editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str=' %.2f+-%.2f'))
    t_pi        = Tuple((np.nan,np.nan), editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str=' %.2f+-%.2f'))
    t_3pi2      = Tuple((np.nan,np.nan), editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str=' %.2f+-%.2f'))
    
    # add fit results to the get_set_items
    get_set_items = FitToolTau.get_set_items + ['contrast','period','t_pi2','t_pi','t_3pi2']
        
    def _update_fit(self):
        y = self.measurement.spin_state
        try:
            fit_result = fitting.fit_rabi( self.measurement.tau, y, y**0.5 )
        except:
            fit_result = (np.NaN*np.zeros(3), np.NaN*np.zeros((3,3)), np.NaN, np.NaN)

        p, v, q, chisqr = fit_result
        a, T, c = p
        a_var, T_var, c_var = v.diagonal()

        # compute some relevant parameters from fit result
        contrast =  200*a/(c+a)
        contrast_delta = 200./c**2*(a_var*c**2+c_var*a**2)**0.5
        T_delta = abs(T_var)**0.5
        pi2         = 0.25*T
        pi          = 0.5*T
        threepi2    = 0.75*T
        pi2_delta       = 0.25*T_delta
        pi_delta        = 0.5*T_delta
        threepi2_delta  = 0.75*T_delta
        
        # set respective attributes
        self.q = q
        self.period = T, T_delta
        self.contrast = contrast, contrast_delta
        self.t_pi2 = pi2, pi2_delta
        self.t_pi = pi, pi_delta
        self.t_3pi2 = threepi2, threepi2_delta

        # create a summary of fit result as a text string
        s = 'q: %.2e\n'%q
        s += 'contrast: %.1f+-%.1f%%\n'%(contrast, contrast_delta)
        s += 'period: %.2f+-%.2f ns\n'%(T, T_delta)
        #s += 'pi/2: %.2f+-%.2f ns\n'%(pi2, pi2_delta)
        #s += 'pi: %.2f+-%.2f ns\n'%(pi, pi_delta)
        #s += '3pi/2: %.2f+-%.2f ns\n'%(threepi2, threepi2_delta)
        
        self.fit_result = fit_result
        self.text = s
                
    def _on_fit_result_change(self, new):
        if len(new) > 0 and new[0][0] is not np.NaN:
            self.line_data.set_data('fit',fitting.Cosinus(*new[0])(self.measurement.tau))             

    traits_view = View(VGroup(VGroup(Item(name='measurement', style='custom', show_label=False),
                                     VGroup(HGroup(Item('contrast',  style='readonly'),
                                                   Item('period',    style='readonly'),
                                                   Item('q',         style='readonly'),
                                                   ),
                                            HGroup(Item('t_pi2',     style='readonly'),
                                                   Item('t_pi',      style='readonly'),
                                                   Item('t_3pi2',    style='readonly'),
                                                   ),
                                            label='fit_result',
                                            ),
                                     ),
                              Tabbed(Item('line_plot', show_label=False, resizable=True),
                                     Item('matrix_plot', show_label=False, resizable=True),
                                     Item('pulse_plot', show_label=False, resizable=True),
                                     ),
                              ),
                       title='Rabi Tool',
                       buttons=[],
                       resizable=True,
                       )


class PulseExtractTool( PulsedToolTau ):

    """
    Provides extraction of pulses from a Rabi measurement.
    """
        
    period      = Float()
    contrast    = Float()
    t_pi2       = Array()
    t_pi        = Array()
    t_3pi2      = Array()
    t_2pi       = Array()
    
    # add fit results to the get_set_items
    get_set_items = PulsedToolTau.get_set_items + ['contrast','period','t_pi2','t_pi','t_3pi2','t_2pi']

    traits_view = View(VGroup(VGroup(Item(name='measurement', style='custom', show_label=False),
                                     Group(VGroup(HGroup(Item('contrast',  style='readonly', format_str='%.1f'),
                                                         Item('period',    style='readonly'),
                                                         ),
                                                  HGroup(Item('t_pi2',     style='readonly'),
                                                         Item('t_pi',      style='readonly'),
                                                         Item('t_3pi2',    style='readonly'),
                                                         Item('t_2pi',     style='readonly'),
                                                         ),
                                                  label='fit_result',
                                                  ),
                                           VGroup(HGroup(Item('integration_width'),
                                                         Item('position_signal'),
                                                         Item('position_normalize'),
                                                         ),
                                                  label='fit_parameter',
                                                  ),
                                           orientation='horizontal', layout='tabbed', springy=False,
                                           ),
                                     ),
                              VSplit(Item('matrix_plot', show_label=False, width=500, height=300, resizable=True),
                                     Item('line_plot', show_label=False, width=500, height=300, resizable=True),
                                     Item('pulse_plot', show_label=False, width=500, height=300, resizable=True),
                                     ),
                              ),
                       title='Pulse Extraction',
                       buttons=[], resizable=True)

    # overwrite __init__ to trigger update events
    def __init__(self, **kwargs):
        super(PulseExtractTool, self).__init__(**kwargs)
        self.on_trait_change(self._update_fit, 'spin_state', dispatch='ui')
    
    def _update_fit(self, y):
        # try to extract pulses from spin_state and tau.
        tau = self.measurement.tau
        try:
            f,r,p,tp=fitting.extract_pulses(y)
        except:
            logging.getLogger().debug('Failed to compute fit.')
            f=[0]
            r=[0]
            p=[0]
            tp=[0]
            
        pi2         = tau[f]
        pi          = tau[p]
        three_pi2   = tau[r]
        two_pi      = tau[tp]

        # compute some relevant parameters from the result
        mi = y.min()
        ma = y.max()
        contrast =  100*(ma-mi)/ma
        # simple rough estimate of the period to avoid index out of range Error
        T = 4*(pi[0]-pi2[0])
        
        # set respective attributes
        self.period = T
        self.contrast = contrast
        self.t_pi2 = pi2
        self.t_pi = pi
        self.t_3pi2 = three_pi2
        self.t_2pi = two_pi

        # create a summary of the result as a text string
        s = 'contrast: %.1f\n'%contrast
        s += 'period: %.2f ns\n'%T
        
        # markers in the plot that show the result of the pulse extraction
        self.line_data.set_data('pulse_indices',np.hstack((pi2, pi, three_pi2, two_pi))*1e-3)            
        self.line_data.set_data('pulse_values',np.hstack((y[f], y[p], y[r], y[tp])))
        self.line_plot.overlays[0].text = s

    # overwrite the line_plot to include fit and text label 
    def _create_line_plot(self):
        line_data   = ArrayPlotData(index=np.array((0,1)),
                                    spin_state=np.array((0,0)),
                                    pulse_indices=np.array((0,0)),
                                    pulse_values=np.array((0,0)))
        plot = Plot(line_data, padding=8, padding_left=64, padding_bottom=36)
        plot.plot(('index','spin_state'), color='blue', name='spin_state')
        plot.plot(('pulse_indices','pulse_values'),
                  type='scatter',
                  marker='circle',
                  color='none',
                  outline_color='red',
                  line_width=1.0,
                  name='pulses')
        plot.index_axis.title = 'time [micro s]'
        plot.value_axis.title = 'spin state'
        plot.overlays.insert(0, PlotLabel(text=self.label_text, hjustify='left', vjustify='bottom', position=[64,32]) )
        self.line_data = line_data
        self.line_plot = plot

if __name__ == '__main__':
    
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')
    
    from tools.emod import JobManager
    
    JobManager().start()

    from hardware.mocking import PulseGenerator, Microwave
    pulser = PulseGenerator()
    microwave = Microwave()

    #from TimeTagger import createMockingTagger
    #tagger = createMockingTagger()

    from hardware.mocking import createMockingTagger
    tagger = createMockingTagger()

    rabi = Rabi(pulser, tagger, microwave)
    rabi_tool = RabiTool(measurement=rabi)
    rabi_tool.edit_traits()
    
