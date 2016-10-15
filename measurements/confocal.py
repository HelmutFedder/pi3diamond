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
import threading
import time
import logging

from traits.api import HasTraits, Trait, Instance, Property, Float, Range,\
                       Bool, Array, String, Str, Enum, Button, on_trait_change, cached_property, DelegatesTo
from traitsui.api import View, Item, Group, HGroup, VGroup, Tabbed, EnumEditor, TextEditor, Action, Menu, MenuBar

from enable.api import ComponentEditor, Component
from chaco.api import CMapImagePlot, ArrayPlotData, DataRange1D,\
                      RdBu, reverse, ColorBar, LinearMapper, DataLabel, PlotLabel
                      
RdBu_r = reverse(RdBu)
#import chaco.api
from chaco.tools.cursor_tool import CursorTool2D

from tools.emod import ManagedJob

#customized zoom tool to keep aspect ratio #
from tools.chaco_addons import SavePlot as Plot, SaveHPlotContainer as HPlotContainer, SaveTool, AspectZoomTool
from tools.utility import GetSetItemsMixin, GetSettableHistory as History

class Confocal( ManagedJob, GetSetItemsMixin ):
    """
    A confocal imaging widget.
    
    It uses a 3 axis piezo scanner combined with some sort
    of photon counting hardware to acquire and display
    images by point scanning. We call the combination
    of the latter two an "imager".
    
    This widget allows you to acquire, save and load
    microscope images, or to move the laser focus around
    manually, using the cross hair or sliders.
    
    you can do x-y images, as well as x-z and y-z images,
    i.e. horizontal cuts as well as cross sections through
    the sample.
    """
    
    # overwrite default priority from ManagedJob (default 0)
    priority = 9

    resolution = Range(low=1, high=1000, value=100, desc='Number of point in long direction', label='resolution', auto_set=False, enter_set=True)
    seconds_per_point = Range(low=1e-3, high=10, value=0.005, desc='Seconds per point [s]', label='Seconds per point [s]', mode='text', auto_set=False, enter_set=True)
    bidirectional = Bool( True )
    return_speed = Range(low=1.0, high=100., value=10., desc='Multiplier for return speed of Scanner if mode is monodirectional', label='return speed', mode='text', auto_set=False, enter_set=True)
    constant_axis = Enum('z', 'x', 'y',
                         label='constant axis',
                         desc='axis that is not scanned when acquiring an image',
                         editor=EnumEditor(values={'x':'1:x','y':'2:y','z':'3:z',},cols=3),)

    # buttons
    history_back    = Button(label='Back')
    history_forward = Button(label='Forward')
    reset_range     = Button(label='reset range')
    reset_cursor    = Button(label='reset position')

    # plot parameters
    thresh_high = Trait( 'auto', Str('auto'), Float(10000.), desc='High Limit of image plot', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float))
    thresh_low = Trait( 'auto', Str('auto'), Float(0.), desc='Low Limit of image plot', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float))
    #colormap = Enum('Spectral','gray')
    show_labels = Bool(False)

    # scan data
    X = Array()
    Y = Array()
    image = Array()

    # plots
    plot_data           = Instance( ArrayPlotData )
    scan_plot           = Instance( CMapImagePlot )
    cursor              = Instance( CursorTool2D )
    zoom                = Instance( AspectZoomTool )
    figure              = Instance( Plot )
    figure_container    = Instance( HPlotContainer, editor=ComponentEditor() )
    z_label_text    = Str('z:0.0')
    cursor_position = Property(depends_on=['x','y','z','constant_axis'])
    
    get_set_items=['constant_axis', 'X', 'Y', 'thresh_high', 'thresh_low', 'seconds_per_point',
                   'return_speed', 'bidirectional', 'history', 'image', 'z_label_text',
                   'resolution', 'x', 'x1', 'x2', 'y', 'y1', 'y2', 'z', 'z1', 'z2']

    def __init__(self, imager, pulser=None, **kwargs):
        super(Confocal, self).__init__(**kwargs)
        self.imager = imager
        self.pulser = pulser
        
        # imager position
        self.add_trait('x', Range(low=imager.get_x_range()[0], high=imager.get_x_range()[1], value=0.5*(imager.get_x_range()[0]+imager.get_x_range()[1]), desc='x [micron]', label='x [micron]', mode='slider'))
        self.add_trait('y', Range(low=imager.get_y_range()[0], high=imager.get_y_range()[1], value=0.5*(imager.get_y_range()[0]+imager.get_y_range()[1]), desc='y [micron]', label='y [micron]', mode='slider'))
        self.add_trait('z', Range(low=imager.get_z_range()[0], high=imager.get_z_range()[1], value=0.5*(imager.get_z_range()[0]+imager.get_z_range()[1]), desc='z [micron]', label='z [micron]', mode='slider'))
    
        # imagging parameters
        self.add_trait('x1', Range(low=imager.get_x_range()[0], high=imager.get_x_range()[1], value=imager.get_x_range()[0], desc='x1 [micron]', label='x1', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%.2f')))
        self.add_trait('y1', Range(low=imager.get_y_range()[0], high=imager.get_y_range()[1], value=imager.get_y_range()[0], desc='y1 [micron]', label='y1', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%.2f')))
        self.add_trait('z1', Range(low=imager.get_z_range()[0], high=imager.get_z_range()[1], value=imager.get_z_range()[0], desc='z1 [micron]', label='z1', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%.2f')))
        self.add_trait('x2', Range(low=imager.get_x_range()[0], high=imager.get_x_range()[1], value=imager.get_x_range()[1], desc='x2 [micron]', label='x2', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%.2f')))
        self.add_trait('y2', Range(low=imager.get_y_range()[0], high=imager.get_y_range()[1], value=imager.get_y_range()[1], desc='y2 [micron]', label='y2', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%.2f')))
        self.add_trait('z2', Range(low=imager.get_z_range()[0], high=imager.get_z_range()[1], value=imager.get_z_range()[1], desc='z2 [micron]', label='z2', editor=TextEditor(auto_set=False, enter_set=True, evaluate=float, format_str='%.2f')))
        
        self.X = np.linspace(imager.get_x_range()[0], imager.get_x_range()[-1], self.resolution+1)
        self.Y = np.linspace(imager.get_y_range()[0], imager.get_y_range()[-1], self.resolution+1)
        self.image = np.zeros((len(self.X), len(self.Y)))
        self._create_plot()
        self.figure.index_range.on_trait_change(self.update_axis_li, '_low_value', dispatch='ui')
        self.figure.index_range.on_trait_change(self.update_axis_hi, '_high_value', dispatch='ui')
        self.figure.value_range.on_trait_change(self.update_axis_lv, '_low_value', dispatch='ui')
        self.figure.value_range.on_trait_change(self.update_axis_hv, '_high_value', dispatch='ui')
        self.zoom.on_trait_change(self.check_zoom, 'box', dispatch='ui')
        self.on_trait_change(self.set_mesh_and_aspect_ratio, 'X,Y', dispatch='ui')
        self.on_trait_change(self.update_image_plot, 'image', dispatch='ui')
        self.sync_trait('cursor_position', self.cursor, 'current_position')
        self.sync_trait('thresh_high', self.scan_plot.value_range, 'high_setting')
        self.sync_trait('thresh_low', self.scan_plot.value_range, 'low_setting')
        self.on_trait_change(self.scan_plot.request_redraw, 'thresh_high', dispatch='ui')
        self.on_trait_change(self.scan_plot.request_redraw, 'thresh_low', dispatch='ui')
        self.history = History(length = 10)
        self.history.put( self.copy_items(['constant_axis', 'X', 'Y', 'image', 'z_label_text', 'resolution'] ) )
        self.labels = {}
        self.label_list = []

    # imager position

    @on_trait_change('x,y,z')
    def _set_imager_position(self):
        if self.state != 'run':
            self.imager.set_position(self.x, self.y, self.z)
    
    @cached_property
    def _get_cursor_position(self):
        if self.constant_axis == 'x':
            return self.z, self.y
        elif self.constant_axis == 'y':
            return self.x, self.z
        elif self.constant_axis == 'z':
            return self.x, self.y
    
    def _set_cursor_position(self, position):
        if self.constant_axis == 'x':
            self.z, self.y = position
        elif self.constant_axis == 'y':
            self.x, self.z = position
        elif self.constant_axis == 'z':
            self.x, self.y = position
    
    # image acquisition
    
    def _run(self):
        """Acquire a scan"""

        try:
            self.state='run'

            self.imager.new_image()
            self.update_mesh()
            X = self.X
            Y = self.Y
            
            XP = X[::-1]
    
            self.image=np.zeros((len(Y),len(X)))
            
            """
            if not self.bidirectional:
                self.imager.initImageScan(len(X), len(Y), self.seconds_per_point, return_speed=self.return_speed)
            else:
                self.imager.initImageScan(len(X), len(Y), self.seconds_per_point, return_speed=None)
            """
            for i,y in enumerate(Y):
                if threading.current_thread().stop_request.isSet():
                    break
                if i%2 != 0 and self.bidirectional:
                    XL = XP
                else:
                    XL = X
                YL = y * np.ones(X.shape)
                
                if self.constant_axis == 'x':
                      const = self.x * np.ones(X.shape)
                      Line = np.vstack( (const, YL, XL) )
                elif self.constant_axis == 'y':
                      const = self.y * np.ones(X.shape)
                      Line = np.vstack( (XL, const, YL) )
                elif self.constant_axis == 'z':
                      const = self.z * np.ones(X.shape)
                      Line = np.vstack( (XL, YL, const) )
                
                if self.bidirectional:
                    c = self.imager.scan_line(Line, self.seconds_per_point)
                else:
                    #start_time = time.time()
                    c = self.imager.scan_line(Line, self.seconds_per_point)
                    #print 'nominal time: '+str(self.seconds_per_point*Line.shape[1])
                    #print 'actual time: '+str(time.time() - start_time)
                    self.imager.scan_line(Line[:,::-1], self.seconds_per_point/self.return_speed)
                if i%2 != 0 and self.bidirectional:
                    self.image[i,:] = c[-1::-1]
                else:
                    self.image[i,:] = c[:]
                
                """
                self.imager.doImageLine(Line)
                self.image = self.imager.getImage()
                """
                self.trait_property_changed('image', self.image)
          
                if self.constant_axis == 'x':
                    self.z_label_text='x:%.2f'%self.x
                elif self.constant_axis == 'y':
                    self.z_label_text='y:%.2f'%self.y
                elif self.constant_axis == 'z':
                    self.z_label_text='z:%.2f'%self.z

            """
            # block at the end until the image is ready 
            if not threading.current_thread().stop_request.isSet(): 
                self.image = self.imager.getImage(1)
                self._image_changed()
            """
                
            self.imager.set_position(self.x, self.y, self.z)

            #save scan data to history
            self.history.put( self.copy_items(['constant_axis', 'X', 'Y', 'image', 'z_label_text', 'resolution'] ) )

        finally:
            self.state = 'idle'

    # plotting
    
    def _create_plot(self):
        plot_data = ArrayPlotData(image=self.image)
        plot = Plot(plot_data, width=500, height=500, resizable='hv', aspect_ratio=1.0, padding=8, padding_left=32, padding_bottom=32)
        plot.img_plot('image',  colormap=RdBu_r, xbounds=(self.X[0],self.X[-1]), ybounds=(self.Y[0],self.Y[-1]), name='image')
        image = plot.plots['image'][0]
        image.x_mapper.domain_limits = (self.imager.get_x_range()[0],self.imager.get_x_range()[1])
        image.y_mapper.domain_limits = (self.imager.get_y_range()[0],self.imager.get_y_range()[1])
        zoom = AspectZoomTool(image, enable_wheel=False)
        cursor = CursorTool2D(image, drag_button='left', color='yellow', marker_size=1.0, line_width=1.0 )
        image.overlays.append(cursor)
        image.overlays.append(zoom)
        colormap = image.color_mapper
        colorbar = ColorBar(index_mapper=LinearMapper(range=colormap.range),
                            color_mapper=colormap,
                            plot=plot,
                            orientation='v',
                            resizable='v',
                            width=16,
                            height=320,
                            padding=8,
                            padding_left=32)
        container = HPlotContainer()
        container.add(plot)
        container.add(colorbar)
        z_label = PlotLabel(text='z=0.0', color='red', hjustify='left', vjustify='bottom', position=[10,10])
        container.overlays.append(z_label)
        container.tools.append(SaveTool(container))

        self.plot_data = plot_data
        self.scan_plot = image
        self.cursor = cursor
        self.zoom = zoom
        self.figure = plot
        self.figure_container = container
        self.sync_trait('z_label_text', z_label, 'text')

    def set_mesh_and_aspect_ratio(self):
        self.scan_plot.index.set_data(self.X,self.Y)
        x1=self.X[0]
        x2=self.X[-1]
        y1=self.Y[0]
        y2=self.Y[-1]
        
        self.figure.aspect_ratio = (x2-x1) / float((y2-y1))
        self.figure.index_range.low = x1
        self.figure.index_range.high = x2
        self.figure.value_range.low = y1
        self.figure.value_range.high = y2

    def check_zoom(self, box):
        li,lv,hi,hv=box
        if self.constant_axis == 'x':
            if not li<self.z<hi:
                self.z = 0.5*(li+hi)
            if not lv<self.y<hv:
                self.y = 0.5*(lv+hv)
        elif self.constant_axis == 'y':
            if not li<self.x<hi:
                self.x = 0.5*(li+hi)
            if not lv<self.z<hv:
                self.z = 0.5*(lv+hv)
        elif self.constant_axis == 'z':
            if not li<self.x<hi:
                self.x = 0.5*(li+hi)
            if not lv<self.y<hv:
                self.y = 0.5*(lv+hv)
        
    def center_cursor(self):
        i = 0.5 * (self.figure.index_range.low + self.figure.index_range.high) 
        v = 0.5 * (self.figure.value_range.low + self.figure.value_range.high) 
        if self.constant_axis == 'x':
            self.z = i
            self.y = v
        elif self.constant_axis == 'y':
            self.x = i
            self.z = v
        elif self.constant_axis == 'z':
            self.x = i
            self.y = v
    
    def _constant_axis_changed(self):
        self.update_mesh()
        self.image = np.zeros((len(self.X), len(self.Y)))  
        self.update_axis()
        self.set_mesh_and_aspect_ratio()

    def update_image_plot(self):
        self.plot_data.set_data('image', self.image)
    
    """
    def _colormap_changed(self, new):
        data = self.figure.datasources['image']
        func = getattr(chaco.api,new)
        self.figure.color_mapper=func(DataRange1D(data))
        self.figure.request_redraw()
    """

    def _show_labels_changed(self, name, old, new):
        for item in self.scan_plot.overlays:
            if isinstance(item, DataLabel) and item.label_format in self.labels:
                item.visible = new
        self.scan_plot.request_redraw()
        
    def get_label_index(self, key):
        for index, item in enumerate(self.scan_plot.overlays):
            if isinstance(item, DataLabel) and item.label_format == key:
                 return index
        return None

    def set_label(self, key, coordinates, **kwargs):

        plot = self.scan_plot

        if self.constant_axis == 'x':
            point = (coordinates[2],coordinates[1])
        elif self.constant_axis == 'y':
            point = (coordinates[0],coordinates[2])
        elif self.constant_axis == 'z':
            point = (coordinates[0],coordinates[1])

        defaults = {'component':plot,
                    'data_point':point,
                    'label_format':key,
                    'label_position':'top right',
                    'bgcolor':'transparent',
                    'text_color':'black',
                    'border_visible':False,
                    'padding_bottom':8,
                    'marker':'cross',
                    'marker_color':'black',
                    'marker_line_color':'black',
                    'marker_line_width':1.5,
                    'marker_size':6,
                    'arrow_visible':False,
                    'clip_to_plot':False,
                    'visible':self.show_labels}

        defaults.update(kwargs)

        label = DataLabel(**defaults)

        index = self.get_label_index(key)
        if index is None:
            plot.overlays.append(label)
        else:
            plot.overlays[index] = label
        self.labels[key] = coordinates
        plot.request_redraw()

    def remove_label(self, key):
        plot = self.scan_plot
        index = self.get_label_index(key)
        plot.overlays.pop(index)
        plot.request_redraw()
        self.labels.pop(key)
        
    def remove_all_labels(self):
        plot = self.scan_plot
        new_overlays = []
        for item in plot.overlays:
            if not ( isinstance(item, DataLabel) and item.label_format in self.labels ) :
                 new_overlays.append(item)
        plot.overlays = new_overlays
        plot.request_redraw()
        self.labels.clear()
        
    @on_trait_change('constant_axis')
    def relocate_labels(self):
        for item in self.scan_plot.overlays:
            if isinstance(item, DataLabel) and item.label_format in self.labels:
                coordinates = self.labels[item.label_format]
                if self.constant_axis == 'x':
                    point = (coordinates[2],coordinates[1])
                elif self.constant_axis == 'y':
                    point = (coordinates[0],coordinates[2])
                elif self.constant_axis == 'z':
                    point = (coordinates[0],coordinates[1])
                item.data_point = point

    def update_axis_li(self):
        if self.constant_axis == 'x':
            self.z1 = self.figure.index_range.low
        elif self.constant_axis == 'y':
            self.x1 = self.figure.index_range.low
        elif self.constant_axis == 'z':
            self.x1 = self.figure.index_range.low 
    def update_axis_hi(self):
        if self.constant_axis == 'x':
            self.z2 = self.figure.index_range.high
        elif self.constant_axis == 'y':
            self.x2 = self.figure.index_range.high
        elif self.constant_axis == 'z':
            self.x2 = self.figure.index_range.high 
    def update_axis_lv(self):
        if self.constant_axis == 'x':
            self.y1 = self.figure.value_range.low
        elif self.constant_axis == 'y':
            self.z1 = self.figure.value_range.low
        elif self.constant_axis == 'z':
            self.y1 = self.figure.value_range.low 
    def update_axis_hv(self):
        if self.constant_axis == 'x':
            self.y2 = self.figure.value_range.high
        elif self.constant_axis == 'y':
            self.z2 = self.figure.value_range.high
        elif self.constant_axis == 'z':
            self.y2 = self.figure.value_range.high
    
    def update_axis(self):
        self.update_axis_li()
        self.update_axis_hi()
        self.update_axis_lv()
        self.update_axis_hv()

    def update_mesh(self):
        if self.constant_axis == 'x':
            x1=self.z1
            x2=self.z2
            y1=self.y1
            y2=self.y2
        elif self.constant_axis == 'y':
            x1=self.x1
            x2=self.x2
            y1=self.z1
            y2=self.z2
        elif self.constant_axis == 'z':
            x1=self.x1
            x2=self.x2
            y1=self.y1
            y2=self.y2

        if (x2-x1) >= (y2-y1):
            self.X = np.linspace(x1,x2,self.resolution)
            self.Y = np.linspace(y1,y2,int(self.resolution*(y2-y1)/(x2-x1)))
        else:
            self.Y = np.linspace(y1,y2,self.resolution)
            self.X = np.linspace(x1,x2,int(self.resolution*(x2-x1)/(y2-y1)))
        
    # GUI buttons

    def _history_back_fired(self):
        self.stop()
        self.set_items( self.history.back() )

    def _history_forward_fired(self):
        self.stop()
        self.set_items( self.history.forward() )

    def _reset_range_fired(self):
        self.x1 = self.imager.get_x_range()[0]
        self.x2 = self.imager.get_x_range()[1]
        self.y1 = self.imager.get_y_range()[0]
        self.y2 = self.imager.get_y_range()[1]
        self.z1 = self.imager.get_z_range()[0]
        self.z2 = self.imager.get_z_range()[1]

    def _reset_cursor_fired(self):
        self.center_cursor()

    # saving images

    def save_image(self, filename=None):
        self.save_figure(self.figure_container, filename)
        
    traits_view = View(VGroup(HGroup(Item('submit_button', show_label=False),
                                     Item('remove_button', show_label=False),
                                     Item('priority'),
                                     Item('state'),
                                     Item('history_back', show_label=False),
                                     Item('history_forward', show_label=False),
                                     ),
                              HGroup(Item('filename',springy=True),
                                     Item('save_button', show_label=False),
                                     Item('load_button', show_label=False)
                                     ),                        
                              Item('figure_container', show_label=False, resizable=True,
                                   enabled_when='state != "run"'
                                   ),
                              HGroup(Item('thresh_low'),
                                     Item('thresh_high'),
                                     #Item('colormap', width=-100),
                                     Item('show_labels'),
                                    ),
                              HGroup(Item('resolution', enabled_when='state != "run"'),
                                     Item('x1'),
                                     Item('x2'),
                                     Item('y1'),
                                     Item('y2'),
                                     Item('z1'),
                                     Item('z2'),
                                     Item('reset_range', show_label=False),
                                     ),
                              HGroup(Item('constant_axis', style='custom', show_label=False, enabled_when='state != "run"'),
                                     Item('bidirectional', enabled_when='state != "run"'),
                                     Item('seconds_per_point'),
                                     Item('return_speed'),
                                     Item('reset_cursor', show_label=False),
                                     ),
                              Item('x', enabled_when='state != "run" or (state == "run" and constant_axis == "x")'),
                              Item('y', enabled_when='state != "run" or (state == "run" and constant_axis == "y")'),
                              Item('z', enabled_when='state != "run" or (state == "run" and constant_axis == "z")'),
                       ),
                       title='Confocal', buttons=[], resizable=True, x=0, y=0
                       )
    
if __name__ == '__main__':

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    from tools.emod import JobManager
    JobManager().start()

    from hardware.mocking import Imager
    imager = Imager()
    
    confocal = Confocal(imager)
    confocal.edit_traits()
    