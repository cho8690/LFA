import numpy as np

_author_ = "A.J. Aranyosi"
_project_ = "colorChecker"
_file_ = "colorCheckerConstants.py"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>

Constants for the ColorChecker project.
"""

# Colors and constants for displaying analysis regions on image in colorChecker.py
CC_BRIGHT = 255
CC_MODERATE = 128
CC_DIM = 64
CC_OFF = 0
CC_BOX_FG = (CC_DIM, CC_DIM, CC_BRIGHT)
CC_FGND_FG = (CC_BRIGHT, CC_DIM, CC_MODERATE)
CC_FGND_BG = (CC_BRIGHT, CC_MODERATE, CC_BRIGHT)
CC_BKGD_FG = (CC_MODERATE, CC_BRIGHT, CC_DIM)
CC_BKGD_BG = (CC_BRIGHT, CC_BRIGHT, CC_DIM)
CC_CORNER_FG = (CC_DIM, CC_BRIGHT, CC_MODERATE)
CC_CORNER_BG = (CC_BRIGHT, CC_BRIGHT, CC_BRIGHT)
CC_PH_CORNER_FG = (CC_DIM, CC_MODERATE, CC_MODERATE)
CC_PH_CORNER_BG = (CC_BRIGHT, CC_DIM, CC_MODERATE)
CC_THRESH_COLOR = (CC_MODERATE, CC_DIM, CC_BRIGHT)
CC_CONTOUR_COLOR = (CC_BRIGHT, CC_DIM, CC_DIM)
CC_GROUP_COLOR = (CC_BRIGHT, CC_OFF, CC_DIM)
CC_VIEW_SCALE = 1024  # scale for determining width of lines in pixels based on image size. This is a hack.
CC_MAX_LINE_WIDTH = 8
CC_MIN_LINE_WIDTH = 1
CC_LINE_LENGTH_SCALE = 10  # scale for determining length of lines for corner relative to line width.
CC_BKGD_THICK_SCALE = 1.25  # Scale of background line thickness to foreground
CC_BKGD_OFFSET_SCALE = 1.1  # Number of channel widths for shifting auto-analysis region from channel to background
CC_REGION_SHIFT_PIXELS = 10  # Number of pixels to shift regions by with arrow keys
CC_REGION_SCALE = 1.03  # Scale value for regions - multiply to scale up, divide to scale down
CC_LARGE_SCALE = 10  # Scale factor for region shifting when holding down shift key
CC_SMALL_SCALE = 0.1  # Scale factor for region shifting when holding down ctrl key
CC_NORMAL_SCALE = 1  # Scale factor for region shifting when not holding down shift or ctrl
CC_ZOOM_IN_SCALE = 1.6  # Scale factor for zooming in
CC_ZOOM_OUT_SCALE = 0.75  # Scale factor for zooming out

# Geometric parameters for the X-Rite Color Checker - do not change!
CC_SQUARE_NX = 6  # Number of squares in x dir
CC_SQUARE_NY = 4  # Number of squares in y dir
CC_ASPECT_RATIO = 1.5  # Aspect ratio of y vs x, in
CC_BORDER_FRAC = 60  # Fraction of width of each border
CC_SQUARE_FRAC = 7.5  # Fraction of width of each square
CC_SPACING_FRAC = 30  # Fraction of width of each inter-square space

# Constants for size of color checker analysis squares and tolerance of click-position accuracy
CC_RECT_FRAC = 2  # Fraction of square to use for analysis quadrilateral
CC_TOLERANCE_FRACTION = 10  # Tolerance of clicking on the correct corner, as a fraction of image dimensions

# General constants
CC_COLOR_CHECKER_PICKLE_FILE = "colorChecker.pickle"
CC_DEFAULT_CORNERS = [(0.68, 0.283), (0.68, 0.692), (0.867, 0.283), (0.867, 0.692)] # relative to image scale
CC_RAW_FILE_EXTENSIONS = ["CR2", "CR3", "DNG", "RAF", "RAW", "CRW", "NEF", "X3F"]
CC_CLICK_TYPES = {
    'corner': 0,
    'fgnd': 1,
    'bkgd': 2,
    'group': 3,
}

# Standard color values/names for the X-Rite Color Checker
CC_REF_COLORS = np.asarray([
    [37.54, 14.37, 14.92],
    [62.73, 35.83, 56.5],
    [28.37, 15.42, -49.8],
    [95.19, -1.03, 2.93],
    [64.66, 19.27, 17.5],
    [39.43, 10.75, -45.17],
    [54.38, -39.72, 32.27],
    [81.29, -0.57, 0.44],
    [49.32, -3.82, -22.54],
    [50.57, 48.64, 16.67],
    [42.43, 51.05, 28.62],
    [66.89, -0.75, -0.06],
    [43.46, -12.74, 22.72],
    [30.1, 22.54, -20.87],
    [81.8, 2.67, 80.41],
    [50.76, -0.13, 0.14],
    [54.94, 9.61, -24.79],
    [71.77, -24.13, 58.19],
    [50.63, 51.28, -14.12],
    [35.63, -0.46, -0.48],
    [70.48, -32.26, -0.37],
    [71.51, 18.24, 67.37],
    [49.57, -29.71, -28.32],
    [20.64, 0.07, -0.46]
])

CC_COLOR_NAMES = [
    "A1 (brown)",
    "A2 (orange)",
    "A3 (blue)",
    "A4 (white)",
    "B1 (peach)",
    "B2 (medium_blue)",
    "B3 (green)",
    "B4 (very_light_gray)",
    "C1 (light_blue)",
    "C2 (light_red)",
    "C3 (red)",
    "C4 (med-light_gray)",
    "D1 (forest_green)",
    "D2 (purple)",
    "D3 (yellow)",
    "D4 (med-dark_gray)",
    "E1 (fuschia)",
    "E2 (light_green)",
    "E3 (pinkish-purple)",
    "E4 (dark_gray)",
    "F1 (turquoise)",
    "F2 (light_orange)",
    "F3 (aqua)",
    "F4 (black)"
]

CC_COMMON_HEADERS = [
    "Min fgnd y",
    "Min fgnd x",
    "Max fgnd y",
    "Max fgnd x",
    "Min bkgd y",
    "Min bkgd x",
    "Max bkgd y",
    "Max bkgd x",
    "Color Channel",
    "Foreground",
    "Background",
]
