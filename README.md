# colorChecker

## Introduction

The colorChecker.py program is designed for measuring and correcting the colors of user-selected targets in an image. The correction is done by measuring the colors of reference targets with known color values. These colors are used to derive an affine transformation between the image color space and the known color space.  This transformation has 12 unknowns, so we need to measure a minimum of four known targets. By using an X-Rite ColorChecker as a reference target, we get 24 targets for a total of 72 measurements, allowing us to perform a least-squares fit to find the best affine transform. All measurements are done in CIELAB color space.

In general, the program measures the median color of a foreground target, a background target, and any number of additional "group" target pairs,
along with the median colors of the ColorChecker swatches. It uses the known reference colors of the ColorChecker (provided by X-Rite based on
spectophotometric measurements) to define a 3D affine transform in CIELAB color space to map the measured colors to the reference color space.
It then reports the foreground color, background color, and difference in both the measured and reference color space for each of the three
L*, a*, and b* channels in CIELAB color space.

A word of caution: colorChecker was developed for internal use, so little effort has been put into making the user interface consistent, reliable, or easy to use. If you run into problems, sometimes you might have to reload the image or even restart the program.

## Requirements

See requirements.txt for the specific packages used for development. However, generally you'll need:

* Python 3.6.x or higher
* imageio
* matplotlib
* numpy
* opencv 3.4.x or higher (4.x is fine). Note that if using conda, do not install conda's opencv package. Instead, use pip install opencv-python.
* pandas
* pathlib
* rawpy
* scikit-image
* scipy
* tkinter

## Usage Overview

This section walks you through the steps required to load an image and measure the output.

1. When you start the app, check that the following settings are available (in order from left to right):

| Option                           | Setting      | Notes                                           |
| -------------------------------- | ------------ | ----------------------------------------------- |
| **Image**                        | N/A          | Load an image to analyze                        |
| **Persistent Groups**            | Checked      | Whether to retain defined groups between images |
| **Clear Last Group**             | N/A          | Clears the last group polygon drawn             |
| **Clear Groups**                 | N/A          | Clears all group polygons drawn                 |
| **Output All**                   | N/A          | Writes measured colors to a file once defined   |

2. Click *Image* to choose an image to load. Use raw images whenever possible to bypass any color corrections
applied by the camera. Processed images such as jpegs may affect the accuracy of measurements.
3. You should see a 4x6 grid of boxes with + symbols outside the corners of the grid. These are used to define
the locations of the color swatches on the ColorChecker. You can drag and drop the + symbols to position the
boxes over the swatches on the ColorChecker. Note that you don't need to keep the orientation of the + symbols
as they are, e.g. you can drag the top-left + down to the lower-right. For fine positioning, simply click
the left mouse button near a + symbol over the pixel you'd like to move it to.
4. Use the mouse wheel to zoom in on your foreground target.
    * The zoomed image is centered on the current mouse location.
    * If you don't have a mouse wheel, you can zoom by holding down Super (Fn on Mac) and left-clicking.
    * Similarly, hold down Super, then Ctrl, then left-click to zoom out. Release Ctrl before releasing Super!
4. Next you'll draw a polygon to define the region that represents your foreground target.
    * Press *r* to toggle into foreground mode.
    * Click on each vertex of your foreground target; it doesn't matter if you go clockwise or counter-clockwise,
    but try to draw simple (non-crossing) polygons. You don't need to click back to your first point.
    * Press *r* again to close the polygon.
    Alternately, you can hold down Shift while clicking on the points; some may find this easier, others more cumbersome.
5. Repeat this process to define your background target, but use *e* instead of *r* (or Ctrl instead of Shift).
6. Repeat this process using *y* instead of *r* or *e* to define additional groups of targets. Note that groups
need to be defined in pairs since the program measures color differences between targets. You can define as many
groups as you like.
6. You can use the arrow keys to fine-tune the positions of your targets, but be aware that all targets will
move together. You can also use *+* and *-* to scale the target sizes, although this will likely move your
targets to a different part of the image.
7. If you need to re-draw any of your targets, just use *r* or *e* again as appropriate. You can also use the
buttons on top to erase the most recently drawn group target, or to erase all group targets.
8. Click *Output All* to generate the results, which will be written into .csv files in the same folder as your
original image. 
    * The image_filename_colors.csv file contains the measured target colors (including the colors
    measured from the ColorChecker). In this file each row corresponds to one target; the foreground target is in
    the first row after the header, the background target in the second, and any group targets following. The ColorChecker
    swatch colors are in the remaining rows. 
     * The image_filename_full_color_analyzed.csv file contains 3 rows for each pair of targets, corresponding
     to the measured L*, a*, and b* color values for the targets. Each row includes the measured colors of the 
     targets, the measured colors of the ColorChecker, the computed affine transform, the corrected colors of the
     targets, the corrected colors of the ColorChecker, and the difference between the foreground and background
     target colors in both measured and corrected color spaces.
9. When you exit the program, the coordinates of the foreground and background target polygons and of the ColorChecker
reference will be saved and loaded again the next time the program is run. The coordinates of group targets are not saved.

## Command Summary

| Command       | Alternate     | Action                                | Notes                     |
| :-----------: | :-----------: | ------------------------------------- | ------------------------- |
| *Super-click* | *Mouse wheel* | Zoom in/out                           | Use Ctrl to toggle in/out |
| *1*           |               | Decrease region shift size            |                           |
| *2*           |               | Increase region shift size            |                           |
| *Up*          |               | Move all regions up                   |                           |
| *Down*        |               | Move all regions down                 |                           |
| *Left*        |               | Move all regions left                 |                           |
| *Right*       |               | Move all regions right                |                           |
| *+*           |               | Increase size of all regions          |                           |
| *-*           |               | Decrease size of all regions          |                           |
| *e*           | *Ctrl-click*  | Define default background region      |                           |
| *r*           | *Shift-click* | Define default channel region         |                           |
| *y*           |               | Define group member region            |                           |

## Tips and Tricks

* If you need to press multiple modifier keys for a command (e.g. Super-Ctrl to zoom out), be sure to press them in the specified order and release them in the opposite order, or the app will get confused. If that happens, press and release them in the correct order to fix it.
If all else fails, quit the program and start again.
