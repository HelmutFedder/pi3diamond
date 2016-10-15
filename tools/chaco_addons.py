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

# Enthought library imports
from traits.api import Float, HasTraits, HasPrivateTraits, Str, Tuple, File, Button
from traitsui.api import Handler, View, Item, OKButton, CancelButton

from exceptions import IOError

class ChacoSaveMixin:
    """
    Provides a 'save' method, that saves the enable component as graphics file.
    """

    def save_raster(self, filename):
        """
        Saves an image of a chaco component (e.g. 'Plot' or 'Container')
        to a raster file, such as .jpg or .png. The file type is terermined
        by the extension.
        """
        from chaco.api import PlotGraphicsContext
        gc = PlotGraphicsContext(self.outer_bounds, dpi=72)
        self.draw(gc, mode="normal")
        #gc.render_component(self)
        gc.save(filename)
        return
    
    def save_pdf(self, filename):
        """
        Saves an image of a chaco component (e.g. 'Plot' or 'Container')
        to a .pdf file.
        """
        from chaco.pdf_graphics_context import PdfPlotGraphicsContext
        gc = PdfPlotGraphicsContext(filename=filename)
                #pagesize = self.pagesize,
                #dest_box = self.dest_box,
                #dest_box_units = self.dest_box_units)
        gc.render_component(self)
        gc.save()
    
    def save(self, filename):
        """
        Saves the plot to a graphics file, e.g. .png or .pdf.
        
        Example of usage:
        
            plot = my_instance.line_plot
            filename = 'foo.png'
            save_figure(plot, filename)
        """
        if filename:
            if os.path.splitext(filename)[-1] == ".pdf":
                self.save_pdf(filename)
            else:
                self.save_raster(filename)
        else:
            raise IOError('Empty filename.')

from chaco.api import Plot, HPlotContainer

class SavePlot(Plot, ChacoSaveMixin):
    pass

class SaveHPlotContainer(HPlotContainer, ChacoSaveMixin):
    pass

# Major library imports
import os.path

# Enthought library imports
from enable.api import BaseTool

from utility import save_file

class SaveTool(BaseTool):
    """
    This tool allows the user to press Ctrl+S to save a snapshot image of
    the plot component.
    """
    #-------------------------------------------------------------------------
    # Override default trait values inherited from BaseTool
    #-------------------------------------------------------------------------

    # This tool does not have a visual representation (overrides BaseTool).
    draw_mode = "none"

    # This tool is not visible (overrides BaseTool).
    visible = False

    def normal_key_pressed(self, event):
        """ Handles a key-press when the tool is in the 'normal' state.

        Saves an image of the plot if the keys pressed are Control and S.
        """
        if self.component is None:
            return

        if event.character == "s" and event.control_down:
            filename = save_file(title='Save plot as png or pdf.')
            if filename:
                self.component.save(filename)
            event.handled = True
        return


from chaco.tools.simple_zoom import SimpleZoom 

class AspectZoomTool(SimpleZoom):

    box = Tuple()

    def _do_zoom(self):
        """ Does the zoom operation.
        """
        # Sets the bounds on the component using _cur_stack_index
        low, high = self._current_state()
        orig_low, orig_high = self._history[0]
    
        if self._history_index == 0:
            if self.tool_mode == "range":
                mapper = self._get_mapper()
                mapper.range.low_setting = self._orig_low_setting
                mapper.range.high_setting = self._orig_high_setting
            else:
                x_range = self.component.x_mapper.range
                y_range = self.component.y_mapper.range
                x_range.low_setting, y_range.low_setting = \
                    self._orig_low_setting
                x_range.high_setting, y_range.high_setting = \
                    self._orig_high_setting

                # resetting the ranges will allow 'auto' to pick the values
                x_range.reset()
                y_range.reset()
               
        else:   
            if self.tool_mode == "range":
                mapper = self._get_mapper()
                if self._zoom_limit_reached(orig_low, orig_high, low, high, mapper):
                    self._pop_state()
                    return
                mapper.range.low = low
                mapper.range.high = high
            else:
                for ndx in (0, 1):
                    mapper = (self.component.x_mapper, self.component.y_mapper)[ndx]
                    if self._zoom_limit_reached(orig_low[ndx], orig_high[ndx],
                                                low[ndx], high[ndx], mapper):
                        # pop _current_state off the stack and leave the actual
                        # bounds unmodified.
                        self._pop_state()
                        return
                x_range = self.component.x_mapper.range
                y_range = self.component.y_mapper.range
                x_range.low, y_range.low = low
                x_range.high, y_range.high = high

        plot = self.component.container
        plot.aspect_ratio = (x_range.high - x_range.low) / (y_range.high - y_range.low)
        
        self.box=(low[0],low[1],high[0],high[1])
        
        self.component.request_redraw()
        return

