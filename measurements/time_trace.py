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

import logging

import numpy as np

from traits.api         import Instance, Int, Float, Bool, Array, Button, on_trait_change
from traitsui.api       import View, Item, HGroup, VGroup

from enable.api         import ComponentEditor
from chaco.api          import ArrayPlotData, Plot
from pyface.timer.api   import Timer

#from TimeTagger import Counter
from hardware.mocking import Counter

from tools.utility import GetSetItemsMixin

class TimeTrace( GetSetItemsMixin ):

    number_of_points = Int(200, desc='Length of Count Trace', label='Number of points', mode='text', auto_set=False, enter_set=True)
    seconds_per_point = Float(0.1, desc='Seconds per point [s]', label='Seconds per point [s]', mode='text', auto_set=False, enter_set=True)

    # trace data
    count_rate = Array()
    time = Array()
    
    enable_0 = Bool(True, label='channel 0', desc='enable channel 0')
    enable_1 = Bool(False, label='channel 1', desc='enable channel 1')
    enable_2 = Bool(False, label='channel 2', desc='enable channel 2')
    enable_3 = Bool(False, label='channel 3', desc='enable channel 3')
    enable_4 = Bool(False, label='channel 4', desc='enable channel 4')
    enable_5 = Bool(False, label='channel 5', desc='enable channel 5')
    enable_6 = Bool(False, label='channel 6', desc='enable channel 6')
    enable_7 = Bool(False, label='channel 7', desc='enable channel 7')
    
    channels = Instance( list, factory=list )

    plot = Instance( Plot )
    plot_data = Instance( ArrayPlotData )
    
    start_button = Button(label='start', show_label=False)
    stop_button = Button(label='stop', show_label=False)
    clear_button = Button(label='clear', show_label=False)

    get_set_items= ['count_rate', 'time', 'channels', 'number_of_points', 'seconds_per_point']

    def __init__(self, tagger, **kwargs):
        super(TimeTrace, self).__init__(**kwargs)
        self.tagger = tagger
        self._reset()

    def _channels_default(self):
        return list(np.arange(8)[np.array([self.enable_0, self.enable_1, self.enable_2, self.enable_3, self.enable_4, self.enable_5, self.enable_6, self.enable_7])])

    @on_trait_change('enable_0,enable_1,enable_2,enable_3,enable_4,enable_5,enable_6,enable_7')
    def _update_channels(self):
        self.channels = self._channels_default()

    @on_trait_change('channels,number_of_points,seconds_per_point')
    def _reset(self):
        if hasattr(self,'_timer'):
            self._timer.Stop()
        self._create_plot()
        self._counter = Counter(self.tagger, self.channels, int(self.seconds_per_point*1e12), self.number_of_points)
        self.time = self._counter.getIndex() * 1e-12
        self.count_rate = self._counter.getData() / self.seconds_per_point
        self._timer = Timer(200, self._refresh_data)
        self._timer.Start()
        
    def _refresh_data(self):
        self.count_rate = self._counter.getData() / self.seconds_per_point

    def _count_rate_changed(self, new):
        for i, line_i in enumerate(new):
            self.plot_data.set_data('channel_'+str(self.channels[i]), line_i)

    def _time_changed(self, new):
        self.plot_data.set_data('time', new)

    def _create_plot(self):
        kwargs = {}
        for i, channel_i in enumerate(self.channels):
            kwargs['channel_'+str(channel_i)] = np.array(())
        data = ArrayPlotData(time=np.array(()), **kwargs)
        plot = Plot(data, width=100, height=100, resizable='hv', padding_left=96, padding_bottom=32)
        plot.index_axis.title = 'time [s]'
        plot.value_axis.title = 'count rate [counts/s]'

        color_map = {0:'blue', 1:'red', 2:'green', 3:'black', 4:'brown', 5:'yellow', 6:'magenta', 7:'cyan'}

        for channel in self.channels:
            plot.plot(('time','channel_'+str(channel)), type='line', color=color_map[channel], name='channel '+str(channel))

        plot.legend.align = 'll'
        if len(self.channels) > 1:
            plot.legend.visible = True
        else:
            plot.legend.visible = False
            
        self.plot_data = data
        self.plot = plot
    
    def _clear_button_fired(self):
        self._counter.clear()
        self.count_rate = self._counter.getData()
    
    def _start_button_fired(self):
        self._counter.start()
        self._timer.Start()
    
    def _stop_button_fired(self):
        self._counter.stop()
        self._timer.Stop()

    def __del__(self):
        self._timer.Stop()

    traits_view = View(VGroup(
                       HGroup(Item('clear_button', show_label=False),
                              Item('start_button', show_label=False),
                              Item('stop_button', show_label=False),
                              Item('save_button', show_label=False),
                              ),
                       HGroup(Item('plot', editor=ComponentEditor(size=(100,100)), show_label=False),
                              VGroup(Item('enable_0'),
                                     Item('enable_1'),
                                     Item('enable_2'),
                                     Item('enable_3'),
                                     Item('enable_4'),
                                     Item('enable_5'),
                                     Item('enable_6'),
                                     Item('enable_7'),
                                     ),
                              ),
                       HGroup(Item('number_of_points'),
                              Item ('seconds_per_point'),
                              ),
                       ),
                       title='Counter',
                       width=780,
                       height=520,
                       buttons=[],
                       resizable=True,
                       )


if __name__=='__main__':

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')
    
    #from TimeTagger import createMockingTagger
    from hardware.mocking import createMockingTagger
    tagger = createMockingTagger()
    
    time_trace = TimeTrace(tagger)
    time_trace.edit_traits()
    
