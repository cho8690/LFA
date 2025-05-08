#!/usr/bin/env python
import tkinter.filedialog
from tkinter import *
import matplotlib

matplotlib.use("TkAgg")  # Need this to avoid iOS crash
from colorCheckerImage import ColorCheckerImage
from colorCheckerSquares import ColorCheckerSquares
from colorCheckerAffine import ColorCheckerAffine
from colorCheckerCalAnalysis import ColorCheckerCalAnalysis
from colorCheckerConstants import *
from colorCalc import colorCalc
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math
import os
import rawpy
from pathlib import Path
import pickle

_author_ = "A.J. Aranyosi"
_project_ = "colorChecker"
_file_ = "colorChecker.py"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>

ColorChecker: the master class for automated analysis of photos containing an X-Rite ColorChecker.
"""


def blackness(color):
    """
    Compute the blackness of a color by multiplying the L* value by the vector length
    of the a* and b* values. Smallest result is blackest.
    """
    return color[0] * math.ceil(math.sqrt(color[1] * color[1] + color[2] * color[2]))


class ColorChecker(object):
    """
    The master class for analyzing ColorChecker images.
    Presents a GUI for loading an image, displays the image, allows selection of the corners of
    the ColorChecker, allows selection of the regions to measure, measures the regions, and writes
    the results to a CSV file. Results include the original measured colors plus the colors after
    applying an affine transform based on the measured and known CIELAB colors of the ColorChecker.
    """

    def __init__(self, master):
        self.master = master
        self.image_frame = None
        self.lineplot_frame = None
        self.image = ColorCheckerImage(self)
        self.squares = ColorCheckerSquares(self)
        self.affine = ColorCheckerAffine(self)
        self.color_calc = colorCalc(self, image=self.image, affine=self.affine)
        self.cal_analysis = ColorCheckerCalAnalysis(self)
        self.corners = None
        self.image_canvas = None
        self.frame = None
        self.control_frame = None
        self.file_button = None
        self.measure_button = None
        self.output_button = None
        self.image_filename = None
        self.click_tolerance = None
        self.current_corner = None
        self.square_colors = None
        self.medians = None
        self.means = None
        self.stdevs = None
        self.mins = None
        self.maxes = None
        self.percentiles = None
        self.A = None
        self.fgnd_polygon = []
        self.bkgd_polygon = []
        self.group_polygons = []
        self.fgnd_new_press = True
        self.bkgd_new_press = True
        self.group_new_press = True
        self.fgnd_set = False
        self.bkgd_set = False
        self.group_set = False
        self.outname = None
        self.outname_addition = ""
        self.persistent_groups = IntVar()
        self.persistent_groups.set(1)
        self.clear_groups_button = None
        self.clear_last_group_button = None
        self.persistent_groups_check = None
        self.zoom_in_scale = CC_ZOOM_IN_SCALE
        self.zoom_out_scale = CC_ZOOM_OUT_SCALE
        self.median_coords = None
        self.key_modifiers = {
            'shift': False,  # foreground selection
            'super': False,  # Fn on Mac - Zoom in/out
            'control': False,  # Background selection, also modifier for zoom
            '1': False,  # Decrease amt to shift regions by with arrow keys
            '2': False,  # Increase amt to shift regions by with arrow keys
            'r': False,  # foreground selection, alternate
            'e': False,  # Background selection, alternate
            'y': False,  # Group selection
            'up': False,  # Move all regions up
            'down': False,  # Move all regions down
            'left': False,  # Move all regions left
            'right': False,  # Move all regions right
            '+': False,  # Scale all region coordinates up
            '-': False,  # Scale all region coordinates down
        }
        self.click_type = None
        self.shift_scale_factor = CC_NORMAL_SCALE  # Scale factor for shifting regions
        self.endpoints = [[0, 0], [0, 0]]
        self.current_coordinate = 0
        self.analysis_output_file = None

        self.load_pickle()
        self.init_gui()

    def init_gui(self):
        """
        Initialize the GUI
        :return: None
        """
        self.master.protocol('WM_DELETE_WINDOW', self.quit)
        self.master.title("ColorChecker Measurement")
        self.frame = Frame(self.master)
        self.frame.pack(anchor=N, fill="both", expand=True)
        self.control_frame = Frame(self.frame, padx=0, pady=0)
        self.control_frame.pack(side=TOP, fill="x", expand=False)
        self.image_frame = Frame(self.frame, padx=0, pady=0)
        self.image_frame.pack(side=BOTTOM, fill="both", expand=True)

        self.file_button = Button(self.control_frame, text='Image', command=self.get_image)
        self.file_button.pack(side=LEFT, padx=3)
        self.measure_button = Button(self.control_frame, text='Output All', command=self.output_all)
        self.measure_button.pack(side=RIGHT, padx=3, pady=0)
        self.clear_groups_button = Button(self.control_frame, text='Clear Groups', command=self.clear_groups)
        self.clear_groups_button.pack(side=RIGHT, padx=6, pady=0)
        self.clear_last_group_button = Button(self.control_frame, text='Clear Last Group', command=self.clear_last_group)
        self.clear_last_group_button.pack(side=RIGHT, padx=3, pady=0)
        self.persistent_groups_check = Checkbutton(self.control_frame, text="Persistent Groups",
                                                 variable=self.persistent_groups)
        self.persistent_groups_check.pack(side=RIGHT, padx=0)

    def quit(self):
        """
        Routine to call when quitting the program
        :return: None
        """
        self.save_pickle()
        self.master.quit()
        self.master.destroy()

    def save_pickle(self):
        with open(CC_COLOR_CHECKER_PICKLE_FILE, 'wb') as f:
            pickle.dump(self.corners, f, -1)
            pickle.dump(self.fgnd_polygon, f, -1)
            pickle.dump(self.bkgd_polygon, f, -1)

    def load_pickle(self):
        try:
            with open(CC_COLOR_CHECKER_PICKLE_FILE, 'rb') as f:
                self.corners = pickle.load(f)
                self.fgnd_polygon = pickle.load(f)
                self.bkgd_polygon = pickle.load(f)
        except FileNotFoundError:
            pass
        except EOFError:
            pass

    def get_image(self):
        """
        Routine for opening and loading an image
        :return: None
        """
        image_filename = tkinter.filedialog.askopenfilename()
        if image_filename:
            self.image_filename = image_filename
            if self.image_filename[-3:].upper() in CC_RAW_FILE_EXTENSIONS:
                no_auto_bright = False
                try:
                    self.image.load_raw_image(self.image_filename, no_auto_bright=no_auto_bright,
                                            use_camera_wb=True, use_auto_wb=False,
                                            use_explicit_wb=False)
                except rawpy._rawpy.LibRawFileUnsupportedError:
                    self.image.load_image(self.image_filename)
            else:
                self.image.load_image(self.image_filename)
            self.master.title(os.path.basename(self.image_filename))
            self.image.draw_image()
            self.reset_key_modifiers()
            self.show_image_canvas()
            self.update_image_canvas()
            self.click_tolerance = max(self.image.image.shape) / CC_TOLERANCE_FRACTION
            self.find_corners()
            self.endpoints = [[0, 0], [0, 0]]
            self.current_coordinate = 0
            self.outname_addition = ""

    def clear_groups(self):
        self.group_polygons = []
        self.update_image_annotations()

    def clear_last_group(self):
        if self.group_polygons is not None and len(self.group_polygons) > 0:
            self.group_polygons = self.group_polygons[:-1]
        self.update_image_annotations()

    def find_corners(self):
        [h, w, _] = self.image.image.shape
        if self.corners is not None:
            for corner in self.corners:
                if not 0 <= corner[0] < h or not 0 <= corner[1] < w:
                    print("Corner {} outside image range ({})! Resetting.".format(corner, self.image.image.shape))
                    self.corners = None
        if self.corners is None:
            self.corners = [[h * b, w * a] for (a, b) in CC_DEFAULT_CORNERS]
        self.update_image_annotations()

    @staticmethod
    def compute_dist(c1, c2):
        """
        Compute the distance between two coordinates
        :param c1: an n-D tuple corresponding to the first coordinate
        :param c2: an n-D tuple corresponding to the second coordinate
        :return: the distance between them
        """
        return math.sqrt(sum([(a - b) * (a - b) for (a, b) in zip(c1, c2)]))

    def is_selected_corner(self, y, x, corner):
        """
        Determines if the x,y coordinate is close enough to be considered the relevant corner.
        Close enough means within 1 / TOLERANCE_FRACTION of the width or height of the corner,
        and not having another corner that's even closer and also within that fraction.
        :param y: y coordinate to be compared
        :param x: x coordinate to be compared
        :param corner: the index of the corner in self.corners to compare to
        :return: True if this coordinate matches, False otherwise
        """
        cor = self.corners[corner]
        cur_xlim = self.image.axis.get_xlim()
        cur_ylim = self.image.axis.get_ylim()
        w = abs(cur_xlim[0] - cur_xlim[1])
        h = abs(cur_ylim[0] - cur_ylim[1])
        if abs(cor[0] - y) > h / CC_TOLERANCE_FRACTION or abs(cor[1] - x) > w / CC_TOLERANCE_FRACTION:
            return False
        ref_dist = self.compute_dist((y, x), cor)
        for c in self.corners:
            dist = self.compute_dist((y, x), c)
            if dist < ref_dist and abs(c[1] - x) < w / CC_TOLERANCE_FRACTION and abs(
                    c[0] - y) < h / CC_TOLERANCE_FRACTION:
                return False
        return True

    def update_selected_corner(self, y, x, corner):
        """
        Updates the coordinates of the specified corner
        :param y: the y coordinate of the new location
        :param x: the x coordinate of the new location
        :param corner: the index into self.corners of the corner to update
        :return: False on error, None otherwise
        """
        self.corners[corner][0] = y
        self.corners[corner][1] = x

    def update_image_annotations(self):
        if self.image.image is None:
            return
        reset_annotated = True
        self.squares.set_corners(self.corners)
        self.image.draw_corners(self.corners, reset_annotated=reset_annotated)
        self.squares.get_analysis_coords()
        reset_annotated = False
        self.image.draw_boxes(self.squares.analysis_coords, reset_annotated=reset_annotated)
        if len(self.bkgd_polygon) > 1:
            close_polygon = False
            if self.bkgd_set is True:
                close_polygon = True
            self.image.draw_polygon(self.bkgd_polygon, color=(CC_BKGD_FG, CC_BKGD_BG),
                                    close_polygon=close_polygon, reset_annotated=reset_annotated)
            reset_annotated = False
        if len(self.fgnd_polygon) > 1:
            close_polygon = False
            if self.fgnd_set is True:
                close_polygon = True
            self.image.draw_polygon(self.fgnd_polygon, color=(CC_FGND_FG, CC_FGND_BG),
                                    close_polygon=close_polygon, reset_annotated=reset_annotated)
            reset_annotated = False
        if self.group_polygons is not None and len(self.group_polygons) > 0:
            for i in range(len(self.group_polygons)):
                close_polygon = True
                if i == len(self.group_polygons) - 1 and self.group_set is False:
                    close_polygon = False
                if len(self.group_polygons[i]) > 1:
                    self.image.draw_polygon(self.group_polygons[i], color=(CC_GROUP_COLOR, CC_GROUP_COLOR),
                                            close_polygon=close_polygon, reset_annotated=reset_annotated)
        self.update_image_canvas(use_transformed=True)

    def measure_colors(self):
        self.squares.measure_colors(use_lab=True)
        # Figure out which corner is our starting point - could be [0][0], [0][-1], [-1][0], [-1][-1],
        # although with some clever logic we could eliminate two of these due to the shape of the
        # colorChecker.
        # Note that our starting corner is the opposite from the blackest corner.
        # Multiply brightness by color distance to determine which one is blackest.
        corner_brights = {(-1, -1): blackness(self.squares.measured_colors[0][0][0]),
                            (1, -1): blackness(self.squares.measured_colors[-1][0][0]),
                            (-1, 1): blackness(self.squares.measured_colors[0][-1][0]),
                            (1, 1): blackness(self.squares.measured_colors[-1][-1][0])}
        xy_dir = sorted(corner_brights, key=lambda x: corner_brights[x])[0]
        if xy_dir[0] == 1:
            xdir = range(CC_SQUARE_NX)
        else:
            xdir = range(CC_SQUARE_NX - 1, -1, -1)
        if xy_dir[1] == 1:
            ydir = range(CC_SQUARE_NY)
        else:
            ydir = range(CC_SQUARE_NY - 1, -1, -1)
        self.medians = []
        self.means = []
        self.stdevs = []
        self.median_coords = []
        self.mins = []
        self.maxes = []
        self.percentiles = []
        for i in xdir:
            for j in ydir:
                self.medians.append(self.squares.measured_colors[i][j][0])
                self.means.append(self.squares.measured_colors[i][j][1])
                self.stdevs.append(self.squares.measured_colors[i][j][2])
                self.median_coords.append(self.squares.analysis_coords[i * CC_SQUARE_NY + j])
                self.mins.append(self.squares.measured_colors[i][j][4])
                self.maxes.append(self.squares.measured_colors[i][j][5])
                self.percentiles.append(self.squares.measured_colors[i][j][6])
        self.squares.measure_foreground(foreground=self.fgnd_polygon, background=self.bkgd_polygon,
                                          use_lab=True)
        if self.group_polygons is not None:
            self.squares.measure_group(group=self.group_polygons, use_lab=True)

    def do_image_output(self, cal_values, reference, extension):
        self.color_calc.set_foregrounds(self.squares.foreground_colors[0])
        self.color_calc.set_backgrounds(self.squares.background_colors[0])
        self.color_calc.set_cal_values(np.asarray(cal_values))
        self.color_calc.compute_affine(reference=reference)
        output_filename = self.image_filename + extension
        if Path(output_filename).is_file():
            output_filename = tkinter.filedialog.asksaveasfilename(
                title="Output file already exists. Choose output file", initialfile=os.path.basename(output_filename))
        if output_filename is not None and output_filename != '':
            imdata = np.reshape(self.affine.apply_affine(data=self.image.lab_image), self.image.lab_image.shape)
            max_scale = [100, -127, 127, -127, 127]
            current_scale = [np.max(imdata[:, :, 0]), np.min(imdata[:, :, 1]), np.max(imdata[:, :, 1]),
                             np.min(imdata[:, :, 2]), np.max(imdata[:, :, 2])]
            # Scale by 0.85 so as not to saturate the output
            best_scale = min([x / y for (x, y) in zip(max_scale, current_scale)]) * 0.85
            print("Writing output image {}".format(output_filename))
            print("Image ranges: L*: {}, {}  a*: {}, {}  b*: {}, {}".format(np.min(imdata[:, :, 0]),
                                                                            np.max(imdata[:, :, 0]),
                                                                            np.min(imdata[:, :, 1]),
                                                                            np.max(imdata[:, :, 1]),
                                                                            np.min(imdata[:, :, 2]),
                                                                            np.max(imdata[:, :, 2])))
            print("Scale: {}".format(best_scale))
            self.image.write_image(output_filename, imdata * best_scale)

    def output_image(self):
        self.measure_colors()
        cal_values = self.medians
        reference = CC_REF_COLORS
        extension = "_corrected.ppm"
        self.do_image_output(cal_values, reference, extension)

    @staticmethod
    def write_color(fh, colors, polygon):
        ys, xs = zip(*polygon)
        fh.write(f"{colors[0][0]},{colors[0][1]},{colors[0][2]},")  # Median color values
        fh.write(f"{min(xs)},{min(ys)},{max(xs)},{max(ys)},")  # Bounding box
        fh.write(f"{colors[1][0]},{colors[1][1]},{colors[1][2]},")  # Mean color values
        fh.write(f"{colors[2][0]},{colors[2][1]},{colors[2][2]},")  # Stdev color values
        fh.write(f"{colors[4][0]},{colors[4][1]},{colors[4][2]},")  # Min color values
        fh.write(f"{colors[5][0]},{colors[5][1]},{colors[5][2]},")  # Max color values
        fh.write(f"{colors[6][0][0]},{colors[6][1][0]},{colors[6][2][0]},")  # 10% color values
        fh.write(f"{colors[6][0][1]},{colors[6][1][1]},{colors[6][2][1]},")  # 25% color values
        fh.write(f"{colors[6][0][2]},{colors[6][1][2]},{colors[6][2][2]},")  # 75% color values
        fh.write(f"{colors[6][0][3]},{colors[6][1][3]},{colors[6][2][3]}")  # 90% color values
        for point in polygon:
            fh.write(",{},{}".format(point[0], point[1]))
        fh.write("\n")

    def output_measured_color_values(self):
        if Path(self.outname).is_file():
            print("Output file {} exists".format(self.outname))
            self.outname = tkinter.filedialog.asksaveasfilename(
                title="Output file already exists. Choose output file", initialfile=os.path.basename(self.outname))
        with open(self.outname, 'w') as fh:
            fh.write("L*,a*,b*,")
            fh.write("x1,y1,x2,y2,")
            fh.write("L* mean,a* mean,b* mean,")
            fh.write("L* stdev,a* stdev,b* stdev,")
            fh.write("L* min,a* min,b* min,")
            fh.write("L* max,a* max,b* max,")
            fh.write("10% L*,10% a*,10% b*,")
            fh.write("25% L*,25% a*,25% b*,")
            fh.write("75% L*,75% a*,75% b*,")
            fh.write("90% L*,90% a*,90% b*,polygon xy values\n")
            self.write_color(fh, self.squares.foreground_colors, self.fgnd_polygon)
            self.write_color(fh, self.squares.background_colors, self.bkgd_polygon)
            if self.squares.group_colors is not None:
                for (color, segment) in zip(self.squares.group_colors, self.group_polygons):
                    self.write_color(fh, color, segment)

    def output_reference_target_color_values(self):
        self.measure_colors()
        self.outname = self.image_filename + "_colors" + self.outname_addition + ".csv"
        self.output_measured_color_values()
        with open(self.outname, 'a') as fh:
            for val in range(len(self.medians)):
                median = self.medians[val]
                mean = self.means[val]
                coords = self.median_coords[val]
                stdevs = self.stdevs[val]
                mins = self.mins[val]
                maxes = self.maxes[val]
                percentiles = self.percentiles[val]
                fh.write(f"{median[0]},{median[1]},{median[2]},")
                fh.write(f"{coords[0][1]},{coords[0][0]},{coords[-1][1]},{coords[-1][0]},")
                fh.write(f"{mean[0]},{mean[1]},{mean[2]},")
                fh.write(f"{stdevs[0]},{stdevs[1]},{stdevs[2]},")
                fh.write(f"{mins[0]},{mins[1]},{mins[2]},")
                fh.write(f"{maxes[0]},{maxes[1]},{maxes[2]},")
                fh.write(f"{percentiles[0][0]},{percentiles[1][0]},{percentiles[2][0]},")
                fh.write(f"{percentiles[0][1]},{percentiles[1][1]},{percentiles[2][1]},")
                fh.write(f"{percentiles[0][2]},{percentiles[1][2]},{percentiles[2][2]},")
                fh.write(f"{percentiles[0][3]},{percentiles[1][3]},{percentiles[2][3]}\n")
        print("Wrote {}".format(self.outname))

    def show_image_canvas(self):
        """
        Shows the matplotlib canvas for the image
        :return: None
        """
        if self.image_canvas is None:
            self.image_canvas = FigureCanvasTkAgg(self.image.f, self.image_frame)
            self.image_canvas.draw()
            self.image_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
            self.image_canvas.mpl_connect('key_press_event', self.key_press_event)
            self.image_canvas.mpl_connect('key_release_event', self.key_release_event)
            self.image_canvas.mpl_connect('button_press_event', self.button_press_event)
            self.image_canvas.mpl_connect('button_release_event', self.button_release_event)
            self.image_canvas.mpl_connect('scroll_event', self.wheel_zoom)
            matplotlib.widgets.Cursor(self.image.axis, horizOn=True, vertOn=True, useblit=True)

    def update_image_canvas(self, redraw=True, use_transformed=False):
        """
        Updates the canvas for the image
        :param redraw: Boolean specifying whether to redraw the underlying image
        :param use_transformed: Boolean specifying whether to redraw using the transformed image
        :return: None
        """
        if redraw is True:
            self.image.draw_image(use_transformed)
        self.image_canvas.draw()

    def transform_point(self, point, scale, offset_x, offset_y):
        [h, w, _] = self.image.image.shape
        y_val = max(0, min(h, int(point[0] * scale + offset_y)))
        x_val = max(0, min(w, int(point[1] * scale + offset_x)))
        return [y_val, x_val]

    def move_selections(self, scale, offset_x, offset_y):
        if self.fgnd_polygon is not None:
            for point in range(len(self.fgnd_polygon)):
                self.fgnd_polygon[point] = self.transform_point(self.fgnd_polygon[point], scale, offset_x, offset_y)
        if self.bkgd_polygon is not None:
            for point in range(len(self.bkgd_polygon)):
                self.bkgd_polygon[point] = self.transform_point(self.bkgd_polygon[point], scale, offset_x, offset_y)
        if self.group_polygons is not None:
            for polygon in range(len(self.group_polygons)):
                for point in range(len(self.group_polygons[polygon])):
                    self.group_polygons[polygon][point] = self.transform_point(self.group_polygons[polygon][point],
                                                                                scale, offset_x, offset_y)

    def key_press_event(self, event):
        redraw = False
        if event.key in self.key_modifiers:
            if event.key == 'up':
                self.move_selections(1.0, 0, 0 - self.shift_scale_factor * CC_REGION_SHIFT_PIXELS)
                redraw = True
            elif event.key == 'down':
                self.move_selections(1.0, 0, self.shift_scale_factor * CC_REGION_SHIFT_PIXELS)
                redraw = True
            elif event.key == 'left':
                self.move_selections(1.0, 0 - self.shift_scale_factor * CC_REGION_SHIFT_PIXELS, 0)
                redraw = True
            elif event.key == 'right':
                self.move_selections(1.0, self.shift_scale_factor * CC_REGION_SHIFT_PIXELS, 0)
                redraw = True
            elif event.key == '+':
                self.move_selections(1 + self.shift_scale_factor * (CC_REGION_SCALE - 1), 0, 0)
                redraw = True
            elif event.key == '-':
                self.move_selections(1.0 / (1 + self.shift_scale_factor * (CC_REGION_SCALE - 1)), 0, 0)
                redraw = True
            elif event.key == 'shift':
                self.fgnd_new_press = True
                self.fgnd_set = False
                self.key_modifiers[event.key] = True
            elif event.key == 'control':
                self.bkgd_new_press = True
                self.bkgd_set = False
                self.key_modifiers[event.key] = True
            elif event.key == 'super':
                self.key_modifiers[event.key] = True
            # For normal keys, system sends press & release events immediately, so we only use press event as toggle
            elif event.key == '1':
                self.shift_scale_factor *= 0.75
                if self.shift_scale_factor < CC_SMALL_SCALE:
                    self.shift_scale_factor = CC_SMALL_SCALE
            elif event.key == '2':
                self.shift_scale_factor *= 1.25
                if self.shift_scale_factor > CC_LARGE_SCALE:
                    self.shift_scale_factor = CC_LARGE_SCALE
            elif event.key == 'r':
                if self.key_modifiers['r'] is False:
                    self.fgnd_new_press = True
                    self.fgnd_set = False
                    self.key_modifiers[event.key] = True
                else:
                    if len(self.fgnd_polygon) > 1:
                        self.fgnd_polygon.append(self.fgnd_polygon[0])
                        self.fgnd_set = True
                        redraw = True
                    else:
                        self.fgnd_polygon = []
                        self.fgnd_new_press = True
                    self.key_modifiers[event.key] = False
            elif event.key == 'e':
                if self.key_modifiers['e'] is False:
                    self.bkgd_new_press = True
                    self.bkgd_set = False
                    self.key_modifiers[event.key] = True
                else:
                    if len(self.bkgd_polygon) > 1:
                        self.bkgd_polygon.append(self.bkgd_polygon[0])
                        self.bkgd_set = True
                        redraw = True
                    else:
                        self.bkgd_polygon = []
                        self.bkgd_new_press = True
                    self.key_modifiers[event.key] = False
            elif event.key == 'y':
                if self.key_modifiers['y'] is False:
                    self.group_new_press = True
                    self.group_set = False
                    self.key_modifiers[event.key] = True
                    if self.group_polygons is None:
                        self.group_polygons = []
                else:
                    if len(self.group_polygons[-1]) > 1:
                        self.group_polygons[-1].append(self.group_polygons[-1][0])
                        self.group_set = True
                        redraw = True
                    else:
                        self.group_polygons.pop()
                        self.group_new_press = True
                    self.key_modifiers[event.key] = False
        if redraw is True:
            self.update_image_annotations()

    def key_release_event(self, event):
        redraw = False
        if event.key in self.key_modifiers:
            if event.key == 'shift':
                self.key_modifiers[event.key] = False
                if len(self.fgnd_polygon) > 1:
                    self.fgnd_polygon.append(self.fgnd_polygon[0])
                    self.fgnd_set = True
                    redraw = True
                else:
                    self.fgnd_polygon = []
                    self.fgnd_new_press = True
            elif event.key == 'control':
                self.key_modifiers[event.key] = False
                if len(self.bkgd_polygon) > 1:
                    self.bkgd_polygon.append(self.bkgd_polygon[0])
                    self.bkgd_set = True
                    redraw = True
                else:
                    self.bkgd_polygon = []
                    self.bkgd_new_press = True
            elif event.key == 'super':
                self.key_modifiers[event.key] = False
            # Since release immediately follows press for p, o, and g, we ignore the release event.
        if redraw is True:
            self.update_image_annotations()

    def button_press_event(self, event):
        """
        Handle button presses
        :param event: the event corresponding to the button press
        :return: None
        """
        if isinstance(event.ydata, type(None)) or isinstance(event.xdata, type(None)):
            print("Click was outside image boundaries")
        else:
            iy, ix = int(round(float(event.ydata))), int(round(float(event.xdata)))
            self.click_type = None
            if event.button == 1:
                if self.key_modifiers['super'] is True:
                    direction = 'up'
                    if self.key_modifiers['control'] is True:
                        direction = 'down'
                    self.do_zoom(direction, event.ydata, event.xdata)
                elif self.key_modifiers['control'] is True or self.key_modifiers['e'] is True:
                    self.click_type = CC_CLICK_TYPES['bkgd']
                elif self.key_modifiers['shift'] is True or self.key_modifiers['r'] is True:
                    self.click_type = CC_CLICK_TYPES['fgnd']
                elif self.key_modifiers['y'] is True:
                    self.click_type = CC_CLICK_TYPES['group']
                else:
                    self.click_type = CC_CLICK_TYPES['corner']
            elif event.button == 2:
                self.click_type = CC_CLICK_TYPES['fgnd']
            elif event.button == 3:
                self.click_type = CC_CLICK_TYPES['bkgd']
            if self.click_type == CC_CLICK_TYPES['corner']:
                self.current_corner = None
                for corner in range(len(self.corners)):
                    if self.is_selected_corner(iy, ix, corner):
                        self.current_corner = corner
            elif self.click_type == CC_CLICK_TYPES['fgnd']:
                if self.fgnd_new_press is True:
                    self.fgnd_polygon = []
                    self.fgnd_new_press = False
                self.fgnd_polygon.append((iy, ix))
            elif self.click_type == CC_CLICK_TYPES['bkgd']:
                if self.bkgd_new_press is True:
                    self.bkgd_polygon = []
                    self.bkgd_new_press = False
                self.bkgd_polygon.append((iy, ix))
            elif self.click_type == CC_CLICK_TYPES['group']:
                if self.group_new_press is True:
                    self.group_polygons.append([])
                    self.group_new_press = False
                self.group_polygons[-1].append((iy, ix))

    def button_release_event(self, event):
        """
        Handle button releases
        :param event: the event corresponding to the button release
        :return: None
        """
        if isinstance(event.ydata, type(None)) or isinstance(event.xdata, type(None)):
            print("Click was outside image boundaries")
        else:
            iy, ix = int(round(float(event.ydata))), int(round(float(event.xdata)))
            if self.click_type == CC_CLICK_TYPES['corner']:
                if self.current_corner is not None:
                    self.update_selected_corner(iy, ix, self.current_corner)
            self.update_image_annotations()

    def wheel_zoom(self, event):
        self.do_zoom(event.button, event.ydata, event.xdata)

    def do_zoom(self, direction, ydata, xdata):
        if xdata is None or ydata is None:
            return False
        if direction == 'up':
            scale_factor = self.zoom_out_scale
        elif direction == 'down':
            scale_factor = self.zoom_in_scale
        else:
            print("Unexpected zoom direction {}".format(direction))
            return False
        cur_xlim = self.image.axis.get_xlim()
        cur_ylim = self.image.axis.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0]) * 0.5
        cur_yrange = (cur_ylim[0] - cur_ylim[1]) * 0.5
        [h, w, _] = self.image.image.shape
        xmin = max(xdata - cur_xrange * scale_factor, 0)
        xmax = min(xdata + cur_xrange * scale_factor, w - 1)
        ymin = max(ydata - cur_yrange * scale_factor, 0)
        ymax = min(ydata + cur_yrange * scale_factor, h - 1)
        self.image.axis.set_xlim([xmin, xmax])
        self.image.axis.set_ylim([ymax, ymin])
        self.update_image_annotations()

    def full_color_analysis(self):
        self.measure_colors()
        references = self.medians
        headers = CC_COMMON_HEADERS + CC_COLOR_NAMES
        extension = "_full_color_analyzed" + self.outname_addition + ".csv"
        self.output_analysis(references, CC_REF_COLORS, headers, extension)

    def color_analysis(self):
        self.full_color_analysis()

    def output_analysis(self, references, swatch_reference, headers, extension):
        foreground = self.squares.foreground_colors[0]
        background = self.squares.background_colors[0]
        fgnd_ys, fgnd_xs = zip(*self.fgnd_polygon)
        bkgd_ys, bkgd_xs = zip(*self.bkgd_polygon)
        lead_in = [min(fgnd_ys), min(fgnd_xs), max(fgnd_ys), max(fgnd_xs),
                    min(bkgd_ys), min(bkgd_xs), max(bkgd_ys), max(bkgd_xs)]
        self.cal_analysis.clear_values()
        self.cal_analysis.add_entry(foreground, background, references, lead_in, headers)
        if self.group_polygons is not None and len(self.group_polygons) > 1 and \
                self.squares.group_colors is not None:
            for pair in range(0, int(len(self.group_polygons) / 2) * 2, 2):
                foreground = self.squares.group_colors[pair][0]
                background = self.squares.group_colors[pair + 1][0]
                fgnd_ys, fgnd_xs = zip(*self.group_polygons[pair])
                bkgd_ys, bkgd_xs = zip(*self.group_polygons[pair + 1])
                lead_in = [min(fgnd_ys), min(fgnd_xs), max(fgnd_ys), max(fgnd_xs),
                            min(bkgd_ys), min(bkgd_xs), max(bkgd_ys), max(bkgd_xs)]
                self.cal_analysis.add_entry(foreground, background, references, lead_in, headers)
        self.cal_analysis.compute_colors(reference=swatch_reference)
        output_filename = self.image_filename + extension
        self.analysis_output_file = output_filename
        if Path(output_filename).is_file():
            output_filename = tkinter.filedialog.asksaveasfilename(
                title="Output file already exists. Choose output file", initialfile=os.path.basename(output_filename))
        self.cal_analysis.write_full_results(output_filename)
        print("Wrote {}".format(output_filename))

    def reset_key_modifiers(self):
        for keymod in self.key_modifiers:
            self.key_modifiers[keymod] = False
        if self.persistent_groups.get() == 0:
            self.group_polygons = None

    def output_all(self):
        self.measure_colors()
        self.output_reference_target_color_values()
        self.color_analysis()


if __name__ == "__main__":
    root = Tk()
    width, height = int(root.winfo_screenwidth() * 7 / 8), int(root.winfo_screenheight() * 3 / 4)
    root.geometry("%dx%d+%d+%d" % (width, height, 100, height / 4 + 50))
    color_checker = ColorChecker(root)
    root.mainloop()
