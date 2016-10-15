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

from traits.api       import Instance, Float, Range, Bool, Array, Str, Enum, Button, on_trait_change, Trait
from traitsui.api     import View, Item, Group, HGroup, VGroup, VSplit, Tabbed, EnumEditor
from enable.api       import ComponentEditor
from chaco.api        import PlotAxis, CMapImagePlot, ColorBar, LinearMapper, ArrayPlotData, RdBu, reverse, PlotLabel
RdBu_r = reverse(RdBu)

from tools.chaco_addons import SaveHPlotContainer as HPlotContainer, SavePlot as Plot, SaveTool 

# date and time tick marks
from chaco.scales.api import CalendarScaleSystem
from chaco.scales_tick_generator import ScalesTickGenerator

import threading
import time
import logging

from tools.emod import ManagedJob, FreeJob
from tools.cron import CronDaemon, CronEvent

from tools.utility import GetSetItemsMixin, warning

from measurements.confocal import Confocal

#from analysis.gaussfitter import gaussfit, twodgaussian, integer_to_mesh
from analysis.fitting import fit, Gaussian

class AutoFocus( ManagedJob, GetSetItemsMixin ):
    """
    An autofocus tool.
    
    It uses an imager and a confocal widget.
    
    To focus in x-y-z, it aquires a small square
    image and finds the maximum intensity and subsequently
    performs a line scan in z-direction and finds the maximum.
    
    It offers several ways to find the maximum,
    including the absolute maximum, gaussian fit
    and spine fit (the latter only for the z-line scan).
    
    The widget can auto_focus periodically and recall
    the image drift.
    
    The widget offers the possibility to mark target points.
    It stores their positions and marks them in the confocal widget.
    """

    # overwrite default priority from ManagedJob (default 0)
    priority = 9

    confocal = Instance( Confocal )

    x1                  = Float(default_value=-1.,   desc='Size of XY Scan',                 label='x1 [micron]',           mode='text',  auto_set=False, enter_set=True)
    x2                  = Float(default_value=1.,    desc='Size of XY Scan',                 label='x2 [micron]',           mode='text',  auto_set=False, enter_set=True)
    y1                  = Float(default_value=-1.,   desc='Size of XY Scan',                 label='y1 [micron]',           mode='text',  auto_set=False, enter_set=True)
    y2                  = Float(default_value=1.,    desc='Size of XY Scan',                 label='y2 [micron]',           mode='text',  auto_set=False, enter_set=True)
    z1                  = Float(default_value=0.,    desc='Size of Z Scan',                  label='z1 [micron]',            mode='text',  auto_set=False, enter_set=True)
    z2                  = Float(default_value=1.,    desc='Size of Z Scan',                  label='z2 [micron]',            mode='text',  auto_set=False, enter_set=True)
    step_xy                 = Float(default_value=0.1,  label='Step XY [micron]',     desc='Step of XY Scan',           mode='text',  auto_set=False, enter_set=True)
    step_z                  = Float(default_value=0.1,  label='Step Z [micron]',      desc='Step of Z Scan',                  mode='text',  auto_set=False, enter_set=True)
    xy_absolute             = Bool(False, label='absolute xy', desc='whether or not to take x1,x2,y1,y2 absolute or relative to the last focus or cross hair.')
    z_absolute              = Bool(False, label='absolute z', desc='whether or not to take z1 and z2 absolute or relative to the last focus or cross hair.')
    seconds_per_point_xy    = Range(low=1e-3, high=10, value=0.05,  desc='Seconds per point for XY Scan',   label='sec/point XY [s]',   mode='text',    auto_set=False, enter_set=True)
    seconds_per_point_z     = Range(low=1e-3, high=10, value=0.2,   desc='Seconds per point for Z Scan',    label='sec/point Z [s]',    mode='text',    auto_set=False, enter_set=True)

    zfit_p = Array()
    zfit_ind = Float()
    zfit_val = Float()

    enable_xy = Bool(True)
    enable_z = Bool(True)

    fit_method_xy = Enum('Maximum', 'Gaussian', desc='Fit Method for XY Scan',    label='XY Fit Method')
    fit_method_z  = Enum('Maximum', 'Gaussian', 'Spline', desc='Fit Method for Z Scan',     label='Z Fit Method')
    smoothing  = Float(default_value=1e8,    desc='spline smoothing parameter',                  label='smoothing',            mode='text',  auto_set=False, enter_set=True)
    
    fit_text = Str('')
    
    X = Array(value=np.array((0.,1.)) )
    Y = Array(value=np.array((0.,1.)) )
    Z = Array(value=np.array((-1.,1.)) )

    data_xy = Array( )
    data_z = Array( value=np.array((0,0)) )

    data_fit_xy = Array( )
    data_fit_z = Array( value=np.array((0,0)) )

    targets         = Instance( {}.__class__, factory={}.__class__ ) # Dict traits are no good for pickling, therefore we have to do it with an ordinary dictionary and take care about the notification manually 
    target_list     = Instance( list, factory=list, args=([None],) ) # list of targets that are selectable in current_target editor
    current_target  = Enum(values='target_list')
        
    drift               = Array( value=np.array(((0,0,0,),)) )
    drift_time          = Array( value=np.array((0,)) )
    current_drift       = Array( value=np.array((0,0,0)) )

    focus_interval    = Range(low=1, high=6000, value=10, desc='Time interval between automatic focus events', label='Interval [m]', auto_set=False, enter_set=True)
    periodic_focus    = Bool(False, label='Periodic focusing')

    target_name = Str(label='name', desc='name to use when adding or removing targets')
    add_target_button       = Button(label='Add Target', desc='add target with given name')
    remove_current_target_button    = Button(label='Remove Current', desc='remove current target')
    remove_all_targets_button    = Button(label='Remove All', desc='remove all targets')
    forget_drift_button    = Button(label='Forget Drift', desc='forget the accumulated drift and reset drift plot')
    next_target_button      = Button(label='Next Target', desc='switch to next available target')
    undo_button       = Button(label='undo', desc='undo the movement of the stage')
    
    previous_state = Instance( () )
    
    plot_data_image = Instance( ArrayPlotData )
    plot_data_line  = Instance( ArrayPlotData )
    plot_data_drift = Instance( ArrayPlotData )
    figure_image    = Instance( HPlotContainer, editor=ComponentEditor() )
    image_label     = Instance( PlotLabel )
    figure_line     = Instance( Plot, editor=ComponentEditor() )
    figure_drift    = Instance( Plot, editor=ComponentEditor() )
    image_plot      = Instance( CMapImagePlot )

    get_set_items=['confocal','targets','current_target','current_drift','drift','drift_time','periodic_focus',
                   'x1', 'x2', 'y1', 'y2', 'z1', 'z2', 'step_xy', 'step_z', 'seconds_per_point_xy', 'seconds_per_point_z',
                   'enable_xy', 'enable_z', 'xy_absolute', 'z_absolute', 'fit_method_xy', 'fit_method_z', 'smoothing', 'fit_text',
                   'data_xy', 'data_z', 'data_fit_xy', 'data_fit_z', 'X', 'Y', 'Z', 'focus_interval' ]
    get_set_order=['confocal','targets']

    def __init__(self, imager, confocal, **kwargs):
        super(AutoFocus, self).__init__(**kwargs)
        self.imager = imager
        self.confocal = confocal
        self.on_trait_change(self.update_plot_image, 'data_xy', dispatch='ui')
        self.on_trait_change(self.update_plot_fit, 'data_fit_xy', dispatch='ui')
        self.on_trait_change(self.update_plot_line_value, 'data_z', dispatch='ui')
        self.on_trait_change(self.update_plot_line_fit, 'data_fit_z', dispatch='ui')
        self.on_trait_change(self.update_plot_line_index, 'Z', dispatch='ui')
        self.on_trait_change(self.update_plot_drift_value, 'drift', dispatch='ui')
        self.on_trait_change(self.update_plot_drift_index, 'drift_time', dispatch='ui')
        self.sync_trait('fit_text', self.image_label, 'text')
    
    @on_trait_change('next_target_button')
    def next_target(self):
        """Convenience method to switch to the next available target."""
        keys = self.targets.keys()
        key = self.current_target
        if len(keys) == 0:
            logging.getLogger().info('No target available. Add a target and try again!')
        elif not key in keys:
            self.current_target = keys[0]
        else:
            self.current_target = keys[(keys.index(self.current_target)+1)%len(keys)]

    def _targets_changed(self, name, old, new):
        l = new.keys() + [None]      # rebuild target_list for Enum trait
        l.sort()
        self.target_list = l
        self._draw_targets()    # redraw target labels

    def _current_target_changed(self):
        self._draw_targets()    # redraw target labels

    def _draw_targets(self):
        c = self.confocal
        c.remove_all_labels()
        c.show_labels=True
        for key, coordinates in self.targets.iteritems():
            if key == self.current_target:
                c.set_label(key, coordinates, marker_color='red')
            else:
                c.set_label(key, coordinates)

    def _periodic_focus_changed(self, new):
        if not new and hasattr(self, 'cron_event'):
            CronDaemon().remove(self.cron_event)
        if new:
            self.cron_event = CronEvent(self.submit, min=range(0,60,self.focus_interval))
            CronDaemon().register(self.cron_event)

    def fit_xy_max(self):
        index = self.data_xy.argmax()
        a = self.data_xy.max()
        xp = self.X[index%len(self.X)]
        yp = self.Y[index/len(self.X)]
        self.fit_text = 'max: {:.0f}\nx: {:.2f}, y: {:.2f})'.format(a,xp,yp)
        self.data_fit_xy = self.data_xy
        return xp, yp
            
    def fit_xy_gaussian(self):
        pars = gaussfit(self.data_xy, circle=True, rotate=False)
        self.fit_pars = pars
        
        xp, yp, width = integer_to_mesh(self.X, self.Y, pars)

        amp = pars[1]

        self.fit_text = 'max: {:.0f}\nx: {:.2f}, y: {:.2f}\ns: {:.2f}'.format(amp,xp,yp,width)

        func = twodgaussian(pars, circle=True, rotate=False)
        data = np.empty_like(self.data_xy)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                data[i,j]=func(i,j)
        self.data_fit_xy = data
        
        return xp, yp
            
    def fit_xy(self):
        if self.fit_method_xy == 'Maximum':
            xp, yp = self.fit_xy_max()                
        elif self.fit_method_xy == 'Gaussian':
            raise NotImplementedError('you may find it in old code')
            #try:
            #    xp, yp = self.fit_xy_gaussian()
            #except: # fall back to maximum
            #    logging.getLogger().warning('Gaussian fit failed. Falling back to Maximum.', exc_info=True)
            #    xp, yp = self.fit_xy_max()
        else:
            raise ValueError('Unknown fit method')
        self.image_plot.request_redraw()
        self.XYFitParameters = [xp, yp]
        self.xfit = xp
        self.yfit = yp
        return xp, yp

    def fit_z(self):
        if self.fit_method_z == 'Maximum':
            zp = self.Z[self.data_z.argmax()]
            self.zfit = zp
            return zp
        elif self.fit_method_z == 'Gaussian':
            x = self.Z
            y = self.data_z
            def GaussianEstimator(x,y):
                n = len(x)
                dx = np.abs(x[1]-x[0])
                c = y.min()
                center = x[y.argmax()]
                amp = y.max()-c
                integral = (y.sum()-n*c)*dx
                width = integral/amp*(2*np.pi)**-0.5
                return c, amp, center, width
            p = fit(x, y, Gaussian, GaussianEstimator)
            self.zfit_p = p
            offset, amp, center, width = p
            self.zfit = center
            self.zfit_ind = center
            self.zfit_val = Gaussian(*p)(center)
            self.data_fit_z = Gaussian(*p)(x)
        elif self.fit_method_z == 'Spline':
            from scipy.interpolate import UnivariateSpline
            sort_map = self.Z.argsort()
            x = self.Z[sort_map]
            y = self.data_z[sort_map]
            
            x1 = x[::2]
            y1 = y[::2]
            x2 = x[1::2]
            y2 = y[1::2]
            
            s1 = UnivariateSpline (x1, y1, k=3, s=self.smoothing)
            s2 = UnivariateSpline (x2, y2, k=3, s=self.smoothing)
            
            ind1 = s1(x1).argmax()
            ind2 = s2(x2).argmax()
            
            x0 = 0.5*(x1[ind1]+x2[ind2])
            y0 = 0.5*(y1[ind1]+y2[ind2])
            
            center = x0
            self.zfit = x0
            self.zfit_ind = x0
            self.zfit_val = y0
            self.data_fit_z = 0.5*(s1(self.Z)+s2(self.Z))
        else:
            raise ValueError('Unknown fit method')
        return center
            

    def add_target(self, key, coordinates=None):
        if coordinates is None:
            c = self.confocal
            coordinates = np.array((c.x,c.y,c.z))
        if self.targets == {}:
            self.forget_drift()
        if self.targets.has_key(key):
            if warning('A target with this name already exists.\nOverwriting will move all targets.\nDo you want to continue?'):
                self.current_drift = coordinates - self.targets[key]
                self.forget_drift()
            else:
                return
        else:
            coordinates = coordinates - self.current_drift
            self.targets[key] = coordinates
        self.trait_property_changed('targets', self.targets)    # trigger event such that Enum is updated and Labels are redrawn
        self.confocal.show_labels=True

    def remove_target(self, key):
        if not key in self.targets:
            logging.getLogger().info('Target cannot be removed. Target does not exist.')
            return
        self.targets.pop(key)        # remove target from dictionary
        self.trait_property_changed('targets', self.targets)    # trigger event such that Enum is updated and Labels are redrawn
        
    def remove_all_targets(self):
        self.targets = {}

    def forget_drift(self):
        targets = self.targets
        # reset coordinates of all targets according to current drift
        for key in targets:
            targets[key] += self.current_drift
        # trigger event such that target labels are redrawn
        self.trait_property_changed('targets', self.targets)
        # set current_drift to 0 and clear plot
        self.current_drift = np.array((0., 0., 0.))
        self.drift_time = np.array((time.time(),))
        self.drift = np.array(((0,0,0),))
        
    def _add_target_button_fired(self):
        self.add_target( self.target_name )
        
    def _remove_current_target_button_fired(self):
        self.remove_target( self.current_target )

    def _remove_all_targets_button_fired(self):
        if warning('Remove all targets. Are you sure?'):
            self.remove_all_targets()

    def _forget_drift_button_fired(self):
        if warning('Forget accumulated drift. Are you sure?'):
            self.forget_drift()

    def _run(self):
        
        logging.getLogger().debug("trying run.")
        
        try:
            self.state='run'
            if self.current_target is None:
                self.focus()
            else: # focus target
                coordinates = self.targets[self.current_target]
                confocal = self.confocal
                confocal.x, confocal.y, confocal.z = coordinates + self.current_drift
                current_coordinates = self.focus()
                self.current_drift = current_coordinates - coordinates  
                self.drift = np.append(self.drift, (self.current_drift,), axis=0)
                self.drift_time = np.append(self.drift_time, time.time())
                logging.getLogger().debug('Drift: %.2f, %.2f, %.2f'%tuple(self.current_drift))
            try:
                self.save('log/auto_focus.pys')
            except:
                pass
        finally:
            self.state = 'idle'

    def focus_xy(self):
            safety = 0
            imager = self.imager
            
            xp = self.confocal.x
            yp = self.confocal.y
            zp = self.confocal.z
            
            if self.xy_absolute:
                xmin = self.x1
                xmax = self.x2
                ymin = self.y1
                ymax = self.y2
            else:
                xmin = np.clip(xp+self.x1, imager.get_x_range()[0]+safety, imager.get_x_range()[1]-safety)
                xmax = np.clip(xp+self.x2+self.step_xy, imager.get_x_range()[0]+safety, imager.get_x_range()[1]-safety)
                ymin = np.clip(yp+self.y1, imager.get_y_range()[0]+safety, imager.get_y_range()[1]-safety)
                ymax = np.clip(yp+self.y2+self.step_xy, imager.get_y_range()[0]+safety, imager.get_y_range()[1]-safety)
            X = np.arange(xmin, xmax, self.step_xy)
            Y = np.arange(ymin, ymax, self.step_xy)
            ZL = zp * np.ones(X.shape)

            self.X = X
            self.Y = Y

            XP = X[::-1]

            self.data_xy=np.zeros((len(Y),len(X)))
            #self.image_plot.index.set_data(X, Y)  
                        
            for i,y in enumerate(Y):
                if threading.current_thread().stop_request.isSet():
                    break
                if i%2 != 0:
                    XL = XP
                else:
                    XL = X
                YL = y * np.ones(X.shape)
                Line = np.vstack( (XL, YL, ZL) )
                
                c = imager.scan_line(Line, self.seconds_per_point_xy) #changed scanLine to scan_line
                if i%2 == 0:
                    self.data_xy[i,:] = c[:]
                else:
                    self.data_xy[i,:] = c[-1::-1]
                
                self.trait_property_changed('data_xy', self.data_xy)
            else:
                xp, yp = self.fit_xy()
                
            self.confocal.x = xp
            self.confocal.y = yp 

            logging.getLogger().info('Focus x,y: %.2f, %.2f' %(xp,yp))
            
    def focus_z(self):
        
            safety = 0
            imager = self.imager
            
            xp = self.confocal.x
            yp = self.confocal.y
            zp = self.confocal.z
            
            if self.z_absolute:
                Z = np.hstack( ( np.arange(zp, self.z1, -self.step_z),
                                    np.arange(self.z1, self.z2, self.step_z),
                                    np.arange(self.z2, zp, -self.step_z) ) )
            else:
                Z = np.hstack( ( np.arange(zp, zp+self.z1, -self.step_z),
                                    np.arange(zp+self.z1, zp+self.z2, self.step_z),
                                    np.arange(zp+self.z2, zp, -self.step_z) ) )
            Z = np.clip(Z, imager.get_z_range()[0]+safety, imager.get_z_range()[1]-safety)

            X = xp * np.ones(Z.shape)
            Y = yp * np.ones(Z.shape)

            if not threading.current_thread().stop_request.isSet():
                Line = np.vstack( (X, Y, Z) )
                data_z = imager.scan_line(Line, self.seconds_per_point_z) #changed scanLine to scan_line

                self.Z = Z
                self.data_z = data_z

                zp = self.fit_z()

            self.confocal.z = zp

            logging.getLogger().info('Focus z: %.2f' %zp)
            
    def focus(self):
            """
            Focuses around current position in x, y, and z-direction.
            """

            self.previous_state = ((self.confocal.x,
                                    self.confocal.y,
                                    self.confocal.z), self.current_target)

            if self.enable_xy:
                self.focus_xy()

            if self.enable_z:
                self.focus_z()

            return self.confocal.x, self.confocal.y, self.confocal.z

    def undo(self):
        if self.previous_state is not None:
            coordinates, target = self.previous_state
            self.confocal.x, self.confocal.y, self.confocal.z = coordinates
            if target is not None:
                self.drift_time = np.delete(self.drift_time, -1)
                self.current_drift = self.drift[-2]
                self.drift = np.delete(self.drift, -1, axis=0)
            self.previous_state = None
        else:
            logging.getLogger().info('Can undo only once.')

    def _undo_button_fired(self):
        self.remove()
        self.undo()
    
    def _plot_data_image_default(self):
        return ArrayPlotData(image=np.zeros((2,2)), fit=np.zeros((2,2)))
    def _plot_data_line_default(self):
        return ArrayPlotData(x=self.Z, y=self.data_z, fit=self.data_fit_z)
    def _plot_data_drift_default(self):
        return ArrayPlotData(t=self.drift_time, x=self.drift[:,0], y=self.drift[:,1], z=self.drift[:,2])

    def _image_label_default(self):
        return PlotLabel(text='', color='red', hjustify='left', vjustify='bottom', position=[10,10])    
        
    def _figure_image_default(self):
        plot = Plot(self.plot_data_image, width=100, height=100, padding=8, padding_left=48, padding_bottom=32)
        plot.img_plot('image', colormap=RdBu_r, name='image')
        plot.aspect_ratio=1
        plot.index_mapper.domain_limits = (self.imager.get_x_range()[0],self.imager.get_x_range()[1])
        plot.value_mapper.domain_limits = (self.imager.get_y_range()[0],self.imager.get_y_range()[1])
        plot_fit = Plot(self.plot_data_image, width=100, height=100, padding=8, padding_left=48, padding_bottom=32)
        plot_fit.img_plot('fit', colormap=RdBu_r, name='fit')
        plot_fit.aspect_ratio=1
        plot_fit.index_mapper.domain_limits = (self.imager.get_x_range()[0],self.imager.get_x_range()[1])
        plot_fit.value_mapper.domain_limits = (self.imager.get_y_range()[0],self.imager.get_y_range()[1])
        container = HPlotContainer()
        image = plot.plots['image'][0]
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
        container.add(plot_fit)
        container.add(colorbar)
        container.overlays.append(self.image_label)        
        container.tools.append(SaveTool(container))
        return container
        
    def _figure_line_default(self):
        plot = Plot(self.plot_data_line, width=100, height=100, padding=8, padding_left=64, padding_bottom=32)
        plot.plot(('x','y'), color='blue')
        plot.plot(('x','fit'), color='red')
        plot.index_axis.title = 'z [micron]'
        plot.value_axis.title = 'Fluorescence [ counts / s ]'
        plot.tools.append(SaveTool(plot))
        return plot
    def _figure_drift_default(self):
        plot = Plot(self.plot_data_drift, width=100, height=100, padding=8, padding_left=64, padding_bottom=32)
        plot.plot(('t','x'), type='line', color='blue', name='x')
        plot.plot(('t','y'), type='line', color='red', name='y')
        plot.plot(('t','z'), type='line', color='green', name='z')
        bottom_axis = PlotAxis(plot,
                               orientation="bottom",
                               tick_generator=ScalesTickGenerator(scale=CalendarScaleSystem()))
        plot.index_axis=bottom_axis
        plot.index_axis.title = 'time'
        plot.value_axis.title = 'drift [micron]'
        plot.legend.visible=True
        plot.tools.append(SaveTool(plot))
        return plot        

    def _image_plot_default(self):
        return self.figure_image.components[0].plots['image'][0]

    def update_plot_image(self):
        self.plot_data_image.set_data('image', self.data_xy)
    def update_plot_fit(self):
        self.plot_data_image.set_data('fit', self.data_fit_xy)
    def update_plot_line_value(self):
        self.plot_data_line.set_data('y', self.data_z)
    def update_plot_line_fit(self):
        self.plot_data_line.set_data('fit', self.data_fit_z)
    def update_plot_line_index(self):
        self.plot_data_line.set_data('x', self.Z)
    def update_plot_drift_value(self):
        if len(self.drift) == 1:
            self.plot_data_drift.set_data('x', np.array(()))
            self.plot_data_drift.set_data('y', np.array(()))
            self.plot_data_drift.set_data('z', np.array(()))            
        else:
            self.plot_data_drift.set_data('x', self.drift[:,0])
            self.plot_data_drift.set_data('y', self.drift[:,1])
            self.plot_data_drift.set_data('z', self.drift[:,2])
    def update_plot_drift_index(self):
        if len(self.drift_time) == 0:
            self.plot_data_drift.set_data('t', np.array(()))
        else:
            self.plot_data_drift.set_data('t', self.drift_time - self.drift_time[0])

    traits_view = View(VGroup(HGroup(Item('submit_button', show_label=False),
                                     Item('remove_button', show_label=False),
                                     Item('priority'),
                                     Item('state', style='readonly'),
                                     Item('undo_button', show_label=False),
                                     ),
                              HGroup(Item('filename',springy=True),
                                     Item('save_button', show_label=False),
                                     Item('load_button', show_label=False)
                                     ),
                              Group(VGroup(HGroup(Item('target_name'),
                                                  Item('add_target_button', show_label=False),
                                                  ),
                                           HGroup(Item('current_target'),
                                                  Item('next_target_button', show_label=False),
                                                  ),
                                           HGroup(Item('remove_all_targets_button', show_label=False),
                                                  Item('remove_current_target_button', show_label=False),
                                                  Item('forget_drift_button', show_label=False),
                                                  ),
                                           HGroup(Item('periodic_focus'),
                                                  Item('focus_interval', enabled_when='not periodic_focus'),
                                                  ),
                                           label='tracking',
                                           ),
                                    VGroup(HGroup(Item('enable_xy'),
                                                  Item('enable_z'),
                                                  Item('xy_absolute'),
                                                  Item('z_absolute'),
                                                  ),
                                           HGroup(Item('fit_method_xy'),
                                                  Item('fit_method_z'),
                                                  ),
                                           HGroup(Item('x1'), Item('x2'),),
                                           HGroup(Item('y1'), Item('y2'),),
                                           HGroup(Item('z1'), Item('z2'),),
                                           HGroup(Item('step_xy'),Item('step_z'),),
                                           HGroup(Item('seconds_per_point_xy'),Item('seconds_per_point_z'),
                                                  ),
                                           label='Settings',
                                           springy=True,
                                           ),
                                    layout='tabbed'
                                    ),
                              VSplit(Item('figure_image', show_label=False, resizable=True),
                                     Item('figure_line', show_label=False, resizable=True),
                                     Item('figure_drift', show_label=False, resizable=True),
                                     ),
                              ),
                       title='Auto Focus', buttons=[], resizable=True
                       )

if __name__ == '__main__':
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    from tools.emod import JobManager
    JobManager().start()

    from hardware.mocking import Imager
    imager = Imager()
    
    from measurements.confocal import Confocal
    confocal = Confocal(imager)
    confocal.edit_traits()
    
    auto_focus = AutoFocus(imager, confocal)
    auto_focus.edit_traits()
    