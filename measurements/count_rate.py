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

import numpy

import threading

# enthought library imports
from traits.api import Instance, Int, Float, Range,\
                       Bool, Array, Enum, Button, on_trait_change
from traitsui.api import Handler, View, Item, Group, HGroup, VGroup, Tabbed, EnumEditor

from enable.api import ComponentEditor
from chaco.api import ArrayPlotData
from tools.chaco_addons import SavePlot as Plot, SaveTool

from tools.utility import GetSetItemsMixin

from tools.emod import Job

#from TimeTagger import Counter
from hardware.mocking import Counter

class StartThreadHandler( Handler ):

    def init(self, info):
        info.object.start()
        
class TimeTrace( Job, GetSetItemsMixin ):

    TraceLength = Range(low=10, high=10000, value=100, desc='Length of Count Trace', label='Trace Length')
    SecondsPerPoint = Range(low=0.001, high=1, value=0.1, desc='Seconds per point [s]', label='Seconds per point [s]')
    RefreshRate = Range(low=0.01, high=1, value=0.1, desc='Refresh rate [s]', label='Refresh rate [s]')

    # trace data
    value = Array()
    index = Array()
    
    c_enable0 = Bool(False, label='channel 0', desc='enable channel 0')
    c_enable1 = Bool(False, label='channel 1', desc='enable channel 1')
    c_enable2 = Bool(False, label='channel 2', desc='enable channel 2')
    c_enable3 = Bool(False, label='channel 3', desc='enable channel 3')
    c_enable4 = Bool(False, label='channel 4', desc='enable channel 4')
    c_enable5 = Bool(False, label='channel 5', desc='enable channel 5')
    c_enable6 = Bool(False, label='channel 6', desc='enable channel 6')
    c_enable7 = Bool(False, label='channel 7', desc='enable channel 7')
    
    TracePlot = Instance( Plot )
    plot_data = Instance( ArrayPlotData )

    def __init__(self, time_tagger, **kwargs):
        super(TimeTrace, self).__init__(**kwargs)
        self.time_tagger = time_tagger
        self._reset()
        #self.on_trait_change(self._index_changed, 'index', dispatch='ui')
        #self.on_trait_change(self._value_changed, 'value', dispatch='ui')
        #self._create_counter()

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
        
    def _create_counter(self):
        channels = []
        for i in range(8):
            if getattr(self,'c_enable%i'%i):
                channels.append(i)
        self.counter = Counter(self.time_tagger, channels, int(self.SecondsPerPoint*1e12), self.TraceLength) 
        self.value = self.counter.getData() / self.seconds_per_point
        self.index = self.counter.getIndex() * 1e-12
        self.channels = channels

        
    def _value_default(self):
        return numpy.zeros((self.TraceLength,))            
    
    def _index_default(self):
        return self.SecondsPerPoint*numpy.arange(self.TraceLength)

    def _index_changed(self):
        self.plot_data.set_data('t', self.index)

    def _value_changed(self, new):
        for i, line_i in enumerate(new):
            self.plot_data.set_data('y%i'%self.channels[i], line_i)

    def _TraceLength_changed(self):
        self._create_counter()
        
    def _SecondsPerPoint_changed(self):
        self._create_counter()

    def _plot_data_default(self):
        return ArrayPlotData(t=self.T, y0=self.C0, y1=self.C1, y2=self.C2, y3=self.C3, y4=self.C4, y5=self.C5, y6=self.C6, y7=self.C7, y8=self.C0C1)
    
    def _TracePlot_default(self):
        plot = Plot(self.plot_data, width=500, height=500, resizable='hv')
        plot.plot(('t','y0'), type='line', color='black')
        plot.tools.append(SaveTool(plot))        
        return plot
    
    @on_trait_change('c_enable0,c_enable1,c_enable2,c_enable3,c_enable4,c_enable5,c_enable6,c_enable7)
    def _replot(self):
        
        self.TracePlot = Plot(self.plot_data, width=500, height=500, resizable='hv')
        self.TracePlot.legend.align = 'll'
        
        n=0
        if self.c_enable0:
            self.TracePlot.plot(('t','y0'), type='line', color='blue',  name='channel 0')
            n+=1
        if self.c_enable1:
            self.TracePlot.plot(('t','y1'), type='line', color='red',   name='channel 1')
            n+=1
        if self.c_enable2:
            self.TracePlot.plot(('t','y2'), type='line', color='green', name='channel 2')
            n+=1
        if self.c_enable3:
            self.TracePlot.plot(('t','y3'), type='line', color='black', name='channel 3')
            n+=1
        if self.c_enable4:
            self.TracePlot.plot(('t','y4'), type='line', color='blue',  name='channel 4')
            n+=1
        if self.c_enable5:
            self.TracePlot.plot(('t','y5'), type='line', color='red',   name='channel 5')
            n+=1
        if self.c_enable6:
            self.TracePlot.plot(('t','y6'), type='line', color='green', name='channel 6')
            n+=1
        if self.c_enable7:
            self.TracePlot.plot(('t','y7'), type='line', color='black', name='channel 7')

        if n > 1:
            self.TracePlot.legend.visible = True
        else:
            self.TracePlot.legend.visible = False

    def _run(self):
        """Acquire Count Trace"""
        while True:
            threading.current_thread().stop_request.wait(self.RefreshRate)
            if threading.current_thread().stop_request.isSet():
                break
            self.value = self.counter.getData() / self.SecondsPerPoint

    traits_view = View( HGroup(Item('TracePlot', editor=ComponentEditor(), show_label=False),
                               #VGroup(Item('c_enable0'),Item('c_enable1'),Item('c_enable2'),Item('c_enable3'),Item('c_enable4'),Item('c_enable5'),Item('c_enable6'),Item('c_enable7'),Item('sum_enable'))
                               VGroup(Item('c_enable0'),Item('c_enable1'),Item('c_enable2'),Item('c_enable3'),Item('c_enable4'),Item('c_enable5'),Item('c_enable6'),Item('sum_enable'))
                        ),
                        Item('TraceLength'),
                        Item ('SecondsPerPoint'),
                        Item ('RefreshRate'),
                        title='Counter', width=800, height=600, buttons=[], resizable=True,
                        handler=StartThreadHandler
                  )


if __name__=='__main__':

    import logging
    
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    from tools.emod import JobManager
    JobManager().start()

    p = PhotonTimeTrace()
    p.edit_traits()
    