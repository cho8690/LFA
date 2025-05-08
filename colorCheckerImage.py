from matplotlib.figure import Figure
import imageio
import skimage.io
from skimage.draw import polygon
from PIL import Image
import rawpy
from colorCheckerConstants import *
import cv2

_author_ = "A.J. Aranyosi"
_project_ = "colorChecker"
_file_ = "colorCheckerImage.py"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>

Class for reading in an image of the X-Rite ColorChecker along with targets whose color is to be measured.
"""


class ColorCheckerImage(object):
    """
    A class for collecting and storing image data and metadata for images
    containing a colorChecker reference target.
    """

    def __init__(self, parent):
        self.parent = parent
        self.image = None  # The image itself
        self.lab_image = None  # The image transformed into Lab color space
        self.annotated_image = None  # The image with annotations added
        self.f = None  # Figure for drawing the image
        self.axis = None  # Axis for the image display
        self.linewidth = None
        self.is_raw_image = False
        self.image_white_balance = []

    @staticmethod
    def write_image(filename, data, is_lab=True):
        """
        Write an image to a file
        :param filename: the filename to write the image to
        :param data: the image data to write
        :param is_lab: a Boolean specifying whether the image data is in Lab color space (vs. RGB)
        :return: None
        """
        if is_lab:
            return skimage.io.imsave(filename, skimage.color.lab2rgb(data))
        else:
            return skimage.io.imsave(filename, data)

    def set_annotated_image(self, image):
        """
        Set the value of a transformed image (e.g. drawing boxes over squares)
        :param image: the transformed image to be assigned to ColorCheckerImage.
        :return: None
        """
        self.annotated_image = None  # for garbage collection
        self.annotated_image = np.copy(image)

    def load_image(self, filename):
        """
        Load an image from a file
        :param filename: The filename containing an image
        :return: None
        """
        self.image = imageio.imread(filename)
        self.reshape_image_and_make_lab()
        self.is_raw_image = False

    def reshape_image_and_make_lab(self):
        """
        Throw out the Alpha channel and convert image from RGB to CIELAB
        Conversion math is documented at:
        https://docs.opencv.org/3.4/de/d25/imgproc_color_conversions.html#color_convert_rgb_lab
        :return: None
        """
        if self.image.shape[2] > 3:  # Throw out the A in RGBA images
            self.image = self.image[:,:,0:3]
        # Note that pixels need to be 32-bit floats in the range [0, 1] per color for this conversion
        self.lab_image = cv2.cvtColor(src=np.array(self.image * 1.0 / 255, dtype=np.float32), code=cv2.COLOR_RGB2LAB)
        [h, w, _] = self.image.shape
        if self.axis is not None:
            self.axis.set_xlim([0, w - 1])
            self.axis.set_ylim([h - 1, 0])

    def load_raw_image(self, filename, no_auto_bright=True, use_camera_wb=True, use_auto_wb=False,
                       use_explicit_wb=None, bright=1.0, user_flip=0):
        """
        Load a raw image from a file
        :param filename: The filename containing a raw image
        :param no_auto_bright: Boolean of whether to auto-brighten the image
        :param use_camera_wb: Boolean of whether to use the camera white balance
        :param use_auto_wb: Boolean of whether to use auto-white balance
        :param use_explicit_wb: If not none, a 4-ple of numbers to pass as WB scale
        :param bright: brightness scale; defaults to 1.0
        :param user_flip: image rotation; 0=none, 3=180, 5=90CCW, 6=90CW. Default is to use orientation from RAW image
        :return: None
        Note: WB scale should be of the form [x, 1.0, y, 0.0] where x and y are the multipliers for red and blue
        """
        with rawpy.imread(filename) as raw:
            try:
                self.image_white_balance = raw.camera_whitebalance[:]
            except rawpy._rawpy.LibRawIOError:
                print(f"Can't get image white balance.")
                self.image_white_balance = []
            if use_explicit_wb not in (None, False):
                # Ignore use_camera_wb and use_auto_wb
                self.image_white_balance = use_explicit_wb[:]
                if no_auto_bright is False:
                    self.image = raw.postprocess(no_auto_bright=no_auto_bright, user_wb=use_explicit_wb,
                                                 output_color=rawpy.ColorSpace.raw,
                                                 demosaic_algorithm=rawpy.DemosaicAlgorithm.VNG, user_flip=user_flip)
                else:
                    self.image = raw.postprocess(no_auto_bright=no_auto_bright, user_wb=use_explicit_wb,
                                                 output_color=rawpy.ColorSpace.raw,
                                                 demosaic_algorithm=rawpy.DemosaicAlgorithm.VNG,
                                                 bright=bright, user_flip=user_flip)
            else:
                if no_auto_bright is False:
                    self.image = raw.postprocess(no_auto_bright=no_auto_bright, use_camera_wb=use_camera_wb,
                                                 use_auto_wb=use_auto_wb,
                                                 output_color=rawpy.ColorSpace.raw,
                                                 demosaic_algorithm=rawpy.DemosaicAlgorithm.VNG, user_flip=user_flip)
                else:
                    self.image = raw.postprocess(no_auto_bright=no_auto_bright, use_camera_wb=use_camera_wb,
                                                 use_auto_wb=use_auto_wb,
                                                 output_color=rawpy.ColorSpace.raw,
                                                 demosaic_algorithm=rawpy.DemosaicAlgorithm.VNG,
                                                 bright=bright, user_flip=user_flip)
        self.reshape_image_and_make_lab()
        self.is_raw_image = True

    def draw_image(self, use_annotated=False):
        """
        Make a plot of an image on the screen
        :param use_annotated: a Boolean specifying whether to use the annotated image
        :return: None
        """
        cur_xlim = cur_ylim = None
        if use_annotated is True:
            image = self.annotated_image
        else:
            image = self.image
        if image is not None:
            if self.f is None:
                self.f = Figure(frameon=False)
            if self.axis is not None:
                cur_xlim = self.axis.get_xlim()
                cur_ylim = self.axis.get_ylim()
            self.f.clf()
            self.axis = self.f.add_axes([0, 0, 1, 1])  # Fill the frame with the image
            if cur_xlim is not None and cur_ylim is not None:
                self.axis.set_xlim(cur_xlim)
                self.axis.set_ylim(cur_ylim)
            self.axis.imshow(image, interpolation='nearest')  # nearest interpolator so we don't blur edges
            self.axis.axis('off')  # Don't show axes for images
            for item in [self.f, self.axis]:
                item.patch.set_visible(False)
            return True
        return False

    def calc_linewidth(self):
        """
        Calculate the linewidth for drawing boxes based on the size and zoom of the image
        :return: the linewidth to use for drawing boxes
        """
        [h, w, _] = self.image.shape
        width_range = min(h, w)
        if self.axis is not None:
            try:
                cur_xlim = self.axis.get_xlim()
                cur_ylim = self.axis.get_ylim()
                width_range = min(abs(cur_xlim[1] - cur_xlim[0]), abs(cur_ylim[1] - cur_ylim[0]))
            except AttributeError:
                pass
        linewidth = int(width_range / CC_VIEW_SCALE)
        if linewidth < CC_MIN_LINE_WIDTH:
            linewidth = CC_MIN_LINE_WIDTH
        if linewidth > CC_MAX_LINE_WIDTH:
            linewidth = CC_MAX_LINE_WIDTH
        return linewidth

    @staticmethod
    def constrain_bounds(coords, shape):
        """
        Constrain coordinates to be within the bounds of the image
        """
        return min(max(int(coords[0]), 0), shape[0] - 1), min(max(int(coords[1]), 0), shape[1] - 1)

    def draw_boxes(self, boxes, color=None, linewidth=None, reset_annotated=False):
        """
        Draws boxes in blue where each color swatch of the Color Checker is located based on the selected corners.
        :param boxes: a tuple of box coordinates, each of which contains a 4-ple of (x,y) coordinates
        :param color: the colors to use for drawing. Should be a tuple of 2 RGB values. Auto-set if color=None
         (note that only the first color is used)
        :param linewidth: the width of the lines to draw. Auto-computed if linewidth=None
        :param reset_annotated: a Boolean specifying whether to reset the annotated image before drawing
        :return: None
        """
        if self.image is None:
            return False
        if reset_annotated is True or self.annotated_image is None:
            self.set_annotated_image(self.image)
        if linewidth is None:
            linewidth = self.calc_linewidth()
        if color is None:
            fg = CC_BOX_FG
        else:
            try:
                fg = color[0]
            except IndexError:
                fg = CC_BOX_FG
        for box in boxes:
            self.draw_line(box[0], box[1], linewidth, fg, self.annotated_image)
            self.draw_line(box[1], box[3], linewidth, fg, self.annotated_image)
            self.draw_line(box[2], box[3], linewidth, fg, self.annotated_image)
            self.draw_line(box[0], box[2], linewidth, fg, self.annotated_image)

    def draw_polygon(self, polygon_tuple, color=None, linewidth=None, close_polygon=False, reset_annotated=False):
        """
        Draw lines corresponding to a multi-point polygon
        :param polygon_tuple: a tuple of (y, x) points representing vertices of the polygon. If the polygon is open
         (i.e., the last point is not identical to the first), it gets closed automatically.
        :param color: the colors to use for drawing. Should be a tuple of 2 RGB values. Auto-set if color=None
         (note that only the first color is used)
        :param linewidth: the width of the lines to draw. Auto-computed if linewidth=None
        :param close_polygon: whether to close the polygon if it's not already closed.
        :param reset_annotated: a Boolean specifying whether to reset the annotated image before drawing
        :return: None
        """
        if self.image is None:
            return False
        if reset_annotated is True or self.annotated_image is None:
            self.set_annotated_image(self.image)
        if linewidth is None:
            linewidth = self.calc_linewidth()
        if color is None:
            fg = CC_BOX_FG
        else:
            try:
                fg = color[0]
            except IndexError:
                fg = CC_BOX_FG
        if polygon_tuple[-1] != polygon_tuple[0] and close_polygon == True:
            polygon_tuple.append(polygon_tuple[0])
        for points in zip(polygon_tuple[:-1], polygon_tuple[1:]):
            self.draw_line(points[0], points[1], linewidth, fg, self.annotated_image)

    def draw_corners(self, corners, linewidth=None, reset_annotated=True, fg=CC_CORNER_FG, bg=CC_CORNER_BG):
        """
        Draws pluses in red over white where corners are for selecting Color Checker bounds.
        :param corners: the corners of the ColorChecker (where the fiducials are)
        :param linewidth: the width of the lines to draw. Auto-computed if linewidth = None
        :param reset_annotated: whether to reset the annotated image to the plain one; default = True
        :param fg: foreground color of line, a tuple of 3 integers
        :param bg: background color of line, a tuple of 3 integers
        :return: None
        """
        if self.image is None:
            return False
        if reset_annotated is True or self.annotated_image is None:
            self.set_annotated_image(self.image)
        if linewidth is None:
            linewidth = self.calc_linewidth()
        for corner in corners:
            # Draw fgd
            for direction in range(2):
                start_loc = corner[:]
                start_loc[direction] -= linewidth * CC_LINE_LENGTH_SCALE / 2
                end_loc = corner[:]
                end_loc[direction] += linewidth * CC_LINE_LENGTH_SCALE / 2
                self.draw_line(start_loc, end_loc, linewidth, fg, self.annotated_image)

    def draw_single_width_line(self, start_loc, end_loc, color, image):
        """
        Draws a single-width line on self.annotated_image from start_loc to end_loc.
        Note that for some reason cv2.line doesn't work here, so we're using a slower method.
        :param start_loc: an (y,x) tuple of the starting pixel of the line
        :param end_loc: an (y,x) tuple of the ending pixel of the line
        :param color: an (B,G,R) tuple of the color to use to draw the line. Each value is 0-255.
        :param image: the image to draw the line onto
        :return: None
        """
        # Let's try a stupider way of doing it
        int_start = self.constrain_bounds(start_loc, image.shape)
        int_end = self.constrain_bounds(end_loc, image.shape)
        step = 1
        ratio = 0
        if abs(end_loc[0] - start_loc[0]) > abs(end_loc[1] - start_loc[1]):  # king-stepping along columns (y axis)
            if end_loc[0] >= start_loc[0]:
                from_loc = int_start
                to_loc = int_end
            else:
                from_loc = int_end
                to_loc = int_start
            if to_loc[1] == from_loc[1]:  # pure vertical
                step = 0
            else:  # Diagonal down to left or rignt
                ratio = 1.0 * (to_loc[0] - from_loc[0]) / (to_loc[1] - from_loc[1])
                if ratio < 0:  # Diagonal down to left
                    step = -1
                    ratio = 0 - ratio
            col = from_loc[1]
            for row in range(from_loc[0], to_loc[0] + 1):
                image[row, col] = color
                if col != to_loc[1]:
                    new_ratio = abs(1.0 * (to_loc[0] - row) / (to_loc[1] - col))
                    if new_ratio < ratio:
                        col += step
                        image[row, col] = color
        else:  # king-stepping along rows (x axis)
            if end_loc[1] >= start_loc[1]:
                from_loc = int_start
                to_loc = int_end
            else:
                from_loc = int_end
                to_loc = int_start
            if to_loc[0] == from_loc[0]:  # pure horizontal
                step = 0
            else:  # Diagonal up or down to right
                ratio = 1.0 * (to_loc[1] - from_loc[1]) / (to_loc[0] - from_loc[0])
                if ratio < 0:  # Diagonal up to right
                    step = -1
                    ratio = 0 - ratio
            row = from_loc[0]
            for col in range(from_loc[1], to_loc[1] + 1):
                image[row, col] = color
                if row != to_loc[0]:
                    new_ratio = abs(1.0 * (to_loc[1] - col) / (to_loc[0] - row))
                    if new_ratio < ratio:
                        row += step
                        image[row, col] = color

    def draw_line(self, start_loc, end_loc, linewidth, color, image):
        """
        Draws lines on self.annotated_image from start_loc to end_loc, with a width of linewidth.
        Note that for some reason cv2.line doesn't work here, so we're using a slower method.
        :param start_loc: an (x,y) tuple of the starting pixel of the line
        :param end_loc: an (x,y) tuple of the ending pixel of the line
        :param linewidth: width of the line, in pixels
        :param color: an (R,G,B) tuple of the color to use to draw the line. Each value is 0-255.
        :param image: the image to draw the line onto
        :return: None
        """
        # cv2.line(image, start_loc, end_loc, color, linewidth)
        half_width = int(linewidth / 2)
        for i in range(0, int(linewidth)):
            self.draw_single_width_line((start_loc[0] + i - half_width, start_loc[1]),
                                   (end_loc[0] + i - half_width, end_loc[1]), color, image)
            self.draw_single_width_line((start_loc[0], start_loc[1] + i - half_width),
                                   (end_loc[0], end_loc[1] + i - half_width), color, image)

    def get_quadrilateral_stats(self, box, use_lab=True):
        """
        Get the median, mean, stdev, and count of pixels inside the box with corners specified by box
        :param box: a list of two (x,y) coordinates of opposite corners of a closed box within which to measure.
        :param use_lab: a Boolean specifying whether to use Lab color space
        :return: a tuple containing the medians, means, standard deviations, and counts of pixels inside the box
        for each channel
        """
        r = [box[0][0], box[0][0], box[-1][0], box[-1][0], box[0][0]]
        c = [box[0][1], box[-1][1], box[-1][1], box[0][1], box[0][1]]
        return self.get_polygon_stats(r, c, use_lab)

    def get_contour_stats(self, contour, use_lab=True):
        """
        Get the median, mean, stdev, and count of pixels inside the box with corners specified by box
        :param contour: a list of the (x,y) coordinates of points along the contour within which to measure.
            Note that code will add the starting point to the end of the list to close the contour
        :param use_lab: a Boolean specifying whether to use Lab color space
        :return: a tuple containing the medians, means, standard deviations, and counts of pixels inside the box
        for each channel
        """
        r = [x[0] for x in contour]
        c = [x[1] for x in contour]
        r.append(contour[0][0])
        c.append(contour[0][1])
        return self.get_polygon_stats(r, c, use_lab)

    def get_polygon_stats(self, r, c, use_lab=True):
        """
        Get the median, mean, stdev, and count of pixels inside the polygon with vertices specified by r and c
        :param r: a list of the row coordinates of the vertices of the polygon
        :param c: a list of the column coordinates of the vertices of the polygon
        :param use_lab: a Boolean specifying whether to use Lab color space
        :return: a tuple containing the medians, means, standard deviations, and counts of pixels inside the polygon
        for each channel
        """
        if use_lab is True:
            image = self.lab_image
        else:
            image = self.image
        if len(r) < 3 or len(c) < 3 or len(r) != len(c):
            print("Wrong number of rows and columns! Got {} rows, {} columns, need at least 3 of each with equal counts".format(r, c))
            median = [0, 0, 0]
            mean = [0, 0, 0]
            stdev = [0, 0, 0]
            n = 0
            min_val = [0, 0, 0]
            max_val = [0, 0, 0]
            percentiles = [0, 0, 0, 0]
        else:
            rr, cc = polygon(r, c, shape=image.shape)
            median = []
            mean = []
            stdev = []
            min_val = []
            max_val = []
            percentiles = []
            for i in range(3):
                median.append(np.median(image[rr, cc, i]))
                mean.append(np.mean(image[rr, cc, i]))
                stdev.append(np.std(image[rr, cc, i]))
                min_val.append(np.min(image[rr, cc, i]))
                max_val.append(np.max(image[rr, cc, i]))
                percentiles.append(np.percentile(image[rr, cc, i], [10, 25, 75, 90]))
            n = len(rr)
        return median, mean, stdev, n, min_val, max_val, percentiles

    def get_subimage_quadrilateral(self, box, use_lab=True):
        """
        Get a subimage containing the pixels within the bounding box as a numpy array
        :param box: a list of (x,y) coordinates of the corners of a closed polygon containing the pixels to return.
        :param use_lab: a Boolean specifying whether to use Lab color space
        :return: a numpy array containing the subimage
        """
        r = [box[0][0], box[0][0], box[-1][0], box[-1][0], box[0][0]]
        c = [box[0][1], box[-1][1], box[-1][1], box[0][1], box[0][1]]
        return self.get_subimage_polygon(r, c, use_lab)

    def get_subimage_polygon(self, r, c, use_lab=True):
        if use_lab is True:
            image = self.lab_image
        else:
            image = self.image
        if len(r) < 3 or len(c) < 3 or len(r) != len(c):
            return [], [], np.asarray([[[]]])
        rr, cc = polygon(r, c, shape=image.shape)
        return rr, cc, np.copy(image[rr, cc, :])
