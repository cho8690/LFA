from colorCheckerImage import ColorCheckerImage
from colorCheckerConstants import *
import math

_author_ = "A.J. Aranyosi"
_project_ = "sweat-imaging"
_file_ = "SweatLineAnalysis.py"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>
"""


class ColorCheckerSquares(object):
    """
    A class for analyzing image values along specific lines.
    This class will be used for orienting and scaling images of sweat patches.
    """

    def __init__(self, parent):
        self.parent = parent
        self.image = parent.image
        self.corner_color = (CC_BRIGHT, CC_DIM, CC_DIM)
        self.corner_bkgd_color = (CC_DIM, CC_BRIGHT, CC_BRIGHT)
        self.box_color = (CC_BRIGHT, CC_BRIGHT, CC_DIM)
        self.box_bkgd_color = (CC_DIM, CC_DIM, CC_BRIGHT)
        self.analysis_coords = None
        self.corners = None
        self.w_coord = 0
        self.h_coord = 1
        self.ws = None
        self.hs = None
        self.measured_colors = None
        self.foreground_colors = None
        self.background_colors = None
        self.group_colors = None

    def set_image(self, image):
        """
        Set the image to be used for analysis.
        :param image: a ColorCheckerImage() to be used for analysis
        :return: False on error, None otherwise
        """
        if not isinstance(image, ColorCheckerImage):
            return False
        self.image = image

    def set_corners(self, corners):
        """
        Sets the corners of the ColorChecker region to analyze.
        :param corners: The corners to use. Should be a tuple of 4 (y,x) coordinates
        :return: False on error, None otherwise
        """
        if self.image is None or corners is None or len(corners) != 4:
            return False
        self.corners = []
        for corner in corners:
            if not hasattr(corner, "__iter__") or len(corner) != 2:
                return False
        # Sort corners so that we get top left, top right, bot left, bot right
        ysort_corners = sorted(corners, key=lambda x: x[0])
        for corner in sorted(ysort_corners[0:2], key=lambda x: x[1]):
            if 0 <= corner[0] < self.image.image.shape[0] and 0 <= corner[1] < self.image.image.shape[1]:
                self.corners.append(corner)
            else:
                print("Corner {} not inside image ({})!".format(corner, self.image.image.shape))
                self.corners = None
                return False
        for corner in sorted(ysort_corners[2:], key=lambda x: x[1]):
            if 0 <= corner[0] < self.image.image.shape[0] and 0 <= corner[1] < self.image.image.shape[1]:
                self.corners.append(corner)
            else:
                print("Corner {} not inside image (part2) ({})!".format(corner, self.image.image.shape))
                self.corners = None
                return False
        if self.corners[2][0] - self.corners[0][0] > self.corners[1][1] - self.corners[0][1]:
            self.w_coord = 0
            self.h_coord = 1
            self.ws = [self.corners[0][0], self.corners[2][0], self.corners[1][0], self.corners[3][0]]
            self.hs = [self.corners[0][1], self.corners[1][1], self.corners[2][1], self.corners[3][1]]
        else:
            self.w_coord = 1
            self.h_coord = 0
            self.ws = [self.corners[0][1], self.corners[1][1], self.corners[2][1], self.corners[3][1]]
            self.hs = [self.corners[0][0], self.corners[2][0], self.corners[1][0], self.corners[3][0]]

    @staticmethod
    def interpolate(x1, x2, frac):
        return x1 + (x2 - x1) * frac

    def get_xy_coords(self, wh_fraction):
        """
        Returns the image coordinates corresponding to a location specified by fractional
        width and height within the ColorChecker. Uses self.ws and self.hs as the width and
        height coordinates of the ColorChecker quadrilateral. Uses self.w_coord and self.h_coord
        to translate these back to image coordinates.
        Method: interpolate width at top and bottom, and height at left and right.
        Then interpolate between top and bottom width based on fractional height, and
        between left and right height based on fractional width.
        :param wh_fraction: a (w,h) pair corresponding to the fractional width and height
            within the ColorChecker; e.g. (0.1, 0.3) would be 10% of the width and 30%
            of the height from the top left corner.
        :return: an (y,x) pair corresponding to the pixel location specified, or False on error.
        """
        if self.image is None:
            print("Image {} is not set".format(self.image))
            return False
        if self.corners is None:
            print("Corners not set! Resetting.")
            self.corners = [[self.image.image.shape[0] * b, self.image.image.shape[1] * a]
                            for (a, b) in CC_DEFAULT_CORNERS]
        xy = [0, 0]
        w_top = self.interpolate(self.ws[0], self.ws[1], wh_fraction[0])
        w_bottom = self.interpolate(self.ws[2], self.ws[3], wh_fraction[0])
        h_left = self.interpolate(self.hs[0], self.hs[1], wh_fraction[1])
        h_right = self.interpolate(self.hs[2], self.hs[3], wh_fraction[1])
        xy[self.w_coord] = self.interpolate(w_top, w_bottom, wh_fraction[1])
        xy[self.h_coord] = self.interpolate(h_left, h_right, wh_fraction[0])
        return xy

    def get_analysis_coords(self):
        """
        Sets self.analysis_coords to a tuple of tuples. The number of tuples corresponds
        to the number of squares (SQUARE_NX * SQUARE_NY); within each tuple, four (x,y)
        pairs specify the coordinates of the quadrilaterals to use for analyzing colors.
        :return: False on error, None otherwise
        """
        self.analysis_coords = []
        wh_fraction = [0, 0]
        square_size = [1.0 / (CC_SQUARE_FRAC * CC_RECT_FRAC), 1.0 / (CC_SQUARE_FRAC * CC_RECT_FRAC) * CC_ASPECT_RATIO]
        wh_fraction[0] = 1.0 / CC_BORDER_FRAC + 1.0 / (CC_SQUARE_FRAC / (1 - 1.0 / CC_RECT_FRAC) * 2)
        for _ in range(CC_SQUARE_NX):
            wh_fraction[1] = 1.0 / CC_BORDER_FRAC * CC_ASPECT_RATIO + \
                             1.0 / (CC_SQUARE_FRAC / (1 - 1.0 / CC_RECT_FRAC) * 2) * CC_ASPECT_RATIO
            for _ in range(CC_SQUARE_NY):
                xy1 = self.get_xy_coords(wh_fraction)
                xy2 = self.get_xy_coords((wh_fraction[0] + square_size[0], wh_fraction[1]))
                xy3 = self.get_xy_coords((wh_fraction[0], wh_fraction[1] + square_size[1]))
                xy4 = self.get_xy_coords([a + b for (a, b) in zip(wh_fraction, square_size)])
                if xy1 is False or xy2 is False or xy3 is False or xy4 is False:
                    print("xy1-4 {} {} {} {} couldn't be set".format(xy1, xy2, xy3, xy4))
                    return False
                self.analysis_coords.append([xy1, xy2, xy3, xy4])
                wh_fraction[1] += 1.0 / CC_SQUARE_FRAC * CC_ASPECT_RATIO + 1.0 / CC_SPACING_FRAC * CC_ASPECT_RATIO
            wh_fraction[0] += 1.0 / CC_SQUARE_FRAC + 1.0 / CC_SPACING_FRAC

    def measure_colors(self, use_lab=True):
        if self.image is None or self.analysis_coords is None:
            return False
        self.measured_colors = []
        for i in range(CC_SQUARE_NX):
            self.measured_colors.append([])
            for j in range(CC_SQUARE_NY):
                self.measured_colors[-1].append(
                    self.image.get_quadrilateral_stats(self.analysis_coords[i * CC_SQUARE_NY + j],
                                                       use_lab=use_lab))

    def measure_foreground(self, foreground=None, background=None, use_lab=True):
        if self.image is None:
            return False
        if foreground is None and hasattr(foreground, "__iter__") and len(foreground) > 2:
            foreground = self.parent.fgnd_polygon
        if background is None and hasattr(background, "__iter__") and len(background) > 2:
            background = self.parent.bkgd_polygon
        r, c = zip(*foreground)
        self.foreground_colors = self.image.get_polygon_stats(r, c, use_lab=use_lab)
        r, c = zip(*background)
        self.background_colors = self.image.get_polygon_stats(r, c, use_lab=use_lab)

    @staticmethod
    def distance(point1, point2):
        diffs = [x - y for (x, y) in zip(point1, point2)]
        return math.sqrt(diffs[0] * diffs[0] + diffs[1] * diffs[1])

    def sort_points(self, region):
        distances = [self.distance(x, y) for x, y in zip(region[:-1], region[1:])]
        if distances[0] + distances[2] > distances[1] + distances[3]:
            mm_scale = (distances[1] + distances[3]) / 2
            if region[0][1] > region[1][1]:  # 1st & 4th point are rightmost
                start_side = (region[0], region[3])
                end_side = (region[1], region[2])
            else:
                start_side = (region[1], region[2])
                end_side = (region[0], region[3])
        else:
            mm_scale = (distances[0] + distances[2]) / 2
            if region[0][1] > region[3][1]:  # 1st & 2nd point are rightmost
                start_side = (region[0], region[1])
                end_side = (region[2], region[3])
            else:
                start_side = (region[2], region[3])
                end_side = (region[0], region[1])
        start_point = min(start_side, key=lambda x: x[0])
        end_point = min(end_side, key=lambda x: x[0])
        start_bottom = max(start_side, key=lambda x: x[0])
        end_bottom = max(end_side, key=lambda x: x[0])
        angle = math.atan2(end_point[0] - start_point[0], start_point[1] - end_point[1])
        bottom_angle = math.atan2(end_bottom[0] - start_bottom[0], start_bottom[1] - end_bottom[1])
        return start_point, end_point, start_bottom, end_bottom, angle, bottom_angle, mm_scale

    def measure_group(self, group=None, use_lab=True):
        if self.image is None or group is None:
            return False
        self.group_colors = []
        for region in group:
            if len(region) > 2:
                r, c = zip(*region)
                self.group_colors.append(self.image.get_polygon_stats(r, c, use_lab=use_lab))
            else:
                self.group_colors.append(0)

    def measure_polygon(self, polygon=None, use_lab=True):
        if self.image is None or polygon is None:
            return False
        if len(polygon) > 2:
            r, c = zip(*polygon)
            return self.image.get_polygon_stats(r, c, use_lab=use_lab)
        else:
            print("Can't measure a region with fewer than 3 points")
            return False
