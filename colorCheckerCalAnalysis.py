import random
import math
import matplotlib
matplotlib.use("TkAgg")  # Need this to avoid iOS crash
import numpy as np
from colorCheckerConstants import CC_COMMON_HEADERS, CC_COLOR_NAMES, CC_REF_COLORS
from colorCalc import colorCalc
from pathlib import Path
import sys
from numbers import Number
import csv

_author_ = "A.J. Aranyosi"
_project_ = "colorChecker"
_file_ = "colorCheckerCalAnalysis.py"

__data_file__ = "calibration-image-values.csv"
__output_file__ = "calibration-image-values-corrected.csv"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>
"""

def compute_affine_distance(affine1, affine2):
    return math.sqrt(sum([sum([(q - w) * (q - w) for (q, w) in zip(x, y)]) for (x, y) in zip(affine1, affine2)]))


class ColorCheckerCalAnalysis(object):
    """
    TODO: insert class doc
    """

    def __init__(self, parent):
        """
        Constructor for colorCheckerCalAnalysis
        """
        self.color_calc = colorCalc(self)
        self.parent = parent
        self.cal_values = []
        self.headers = None
        self.fgnd_entry = None
        self.color_vals = []
        self.references = None

    def load_values(self, filename=__data_file__):
        with open(filename) as file:
            self.headers = file.readline().strip().split(",")
            self.fgnd_entry = self.headers.index("Foreground")
            l_vals = None
            a_vals = None
            for line in file:
                entries = line.strip().split(",")
                if entries[self.fgnd_entry - 1] == "L*":
                    self.cal_values.append({})
                    self.cal_values[-1]["lead-in"] = []
                    for i in range(self.fgnd_entry - 1):
                        self.cal_values[-1]["lead-in"].append(entries[i])
                    l_vals = [float(x) for x in entries[self.fgnd_entry:]]
                elif entries[self.fgnd_entry - 1] == "a*":
                    a_vals = [float(x) for x in entries[self.fgnd_entry:]]
                elif entries[self.fgnd_entry - 1] == "b*":
                    b_vals = [float(x) for x in entries[self.fgnd_entry:]]
                    self.cal_values[-1]["foreground"] = [l_vals.pop(0), a_vals.pop(0), b_vals.pop(0)]
                    self.cal_values[-1]["background"] = [l_vals.pop(0), a_vals.pop(0), b_vals.pop(0)]
                    self.cal_values[-1]["medians"] = [[x, y, z] for (x, y, z) in zip(l_vals, a_vals, b_vals)]
                else:
                    print("Invalid image channel type {} in input file".format(entries[self.fgnd_entry - 1]))

    def load_from_colors_file(self, filename=__data_file__, nrefs=4):
        l_vals = []
        a_vals = []
        b_vals = []
        x_min = []
        y_min = []
        x_max = []
        y_max = []
        with open(filename) as file:
            reader = csv.DictReader(file)
            for row in reader:
                l_vals.append(float(row["L*"]))
                a_vals.append(float(row["a*"]))
                b_vals.append(float(row["b*"]))
                x_min.append(float(row["x1"]))
                y_min.append(float(row["y1"]))
                x_max.append(float(row["x2"]))
                y_max.append(float(row["y2"]))
        self.headers = CC_COMMON_HEADERS + CC_COLOR_NAMES
        self.fgnd_entry = self.headers.index("Foreground")
        self.references = [[x, y, z] for (x, y, z) in zip(l_vals[-nrefs:], a_vals[-nrefs:], b_vals[-nrefs:])]
        for pair in range(0, len(l_vals) - nrefs, 2):
            channel = [l_vals[pair], a_vals[pair], b_vals[pair]]
            background = [l_vals[pair + 1], a_vals[pair + 1], b_vals[pair + 1]]
            lead_in = [x_min[pair], y_min[pair], x_max[pair], y_max[pair], x_min[pair + 1], y_min[pair + 1],
                        x_max[pair + 1], y_max[pair + 1]]
            self.add_entry(foreground=channel, background=background, references=self.references, lead_in=lead_in)

    def add_entry(self, foreground=None, background=None, references=None, lead_in=(), headers=None):
        if headers is not None:
            self.headers = headers
            self.fgnd_entry = self.headers.index("Foreground")
        if foreground is None or not hasattr(foreground, '__iter__') or len(foreground) != 3:
            print("foreground colors {} not valid".format(foreground))
            return False
        if background is None or not hasattr(background, '__iter__') or len(background) != 3:
            print("background colors {} not valid".format(background))
            return False
        if references is None or not hasattr(references, '__iter__') or not (len(references) == 4 or len(references) == 24):
            print("reference colors {} not valid".format(references))
            return False
        for reference in references:
            if reference is None or not hasattr(reference, '__iter__') or len(reference) != 3:
                print("reference color {} not valid".format(reference))
                return False
        self.cal_values.append({})
        self.cal_values[-1]["lead-in"] = []
        for entry in lead_in:
            self.cal_values[-1]["lead-in"].append(entry)
        self.cal_values[-1]["foreground"] = foreground
        self.cal_values[-1]["background"] = background
        self.cal_values[-1]["medians"] = references

    def clear_values(self):
        self.cal_values = []

    def arrays_equal(self, a, b):
        if not hasattr(a, '__iter__'):
            if not hasattr(b, '__iter__'):
                if a != b:
                    return False
                return True
            return False
        if not hasattr(b, '__iter__'):
            return False
        if len(a) != len(b):
            return False
        for i in range(len(a)):
            if not self.arrays_equal(a[i], b[i]):
                return False
        return True

    def compute_colors(self, reference=None):
        if self.cal_values is None or len(self.cal_values) == 0:
            return False
        for cal in self.cal_values:
            if "foreground" in cal:
                self.color_calc.set_foregrounds([cal["foreground"]])
            else:
                continue
            if "background" in cal:
                self.color_calc.set_backgrounds([cal["background"]])
            else:
                continue
            if "medians" in cal:
                self.color_calc.set_cal_values(cal["medians"])
            else:
                continue
            if reference is None:
                reference = CC_REF_COLORS
            self.color_calc.compute_affine(reference=reference)
            self.color_calc.apply_affine()
            self.color_calc.measure_colors()
            cal["affine"] = self.color_calc.affine.a_matrix.T.tolist()
            cal["corrected_foreground"] = self.color_calc.corrected_foregrounds[0][:]
            cal["corrected_background"] = self.color_calc.corrected_backgrounds[0][:]
            cal["corrected_medians"] = self.color_calc.corrected_cal_values[:]
            cal["diffs"] = self.color_calc.differences[0][:]
            cal["corr_diffs"] = self.color_calc.corrected_differences[0][:]

    def write_results(self, filename=__output_file__):
        with open(filename, 'w') as file:
            if self.headers is not None:
                for header in self.headers:
                    file.write("{},".format(header))
            file.write("Affine L* scale,Affine a* scale,Affine b* scale,Affine constant,")
            for header in self.headers[self.fgnd_entry:]:
                file.write("{} corrected,".format(header))
            file.write("Difference,Corrected Difference,")
            file.write("\n")
            for cal in self.cal_values:
                # Sequentially write L*, a*, b* data
                colors = {0: "L*", 1: "a*", 2: "b*"}
                for color in colors:
                    for i in range(self.fgnd_entry - 1):
                        file.write("{},".format(cal["lead-in"][i]))
                    file.write("{},".format(colors[color]))
                    file.write("{},{},".format(cal["foreground"][color], cal["background"][color]))
                    for i in range(len(cal["medians"])):
                        file.write("{},".format(cal["medians"][i][color]))
                    for i in range(4):
                        file.write("{},".format(cal["affine"][color][i]))
                    file.write("{},{},".format(cal["corrected_foreground"][color],
                                               cal["corrected_background"][color]))
                    for i in range(len(cal["corrected_medians"])):
                        file.write("{},".format(cal["corrected_medians"][i][color]))
                    file.write("{},{},".format(cal["diffs"][color], cal["corr_diffs"][color]))
                    file.write("\n")

    def write_full_results(self, filename=__output_file__):
        with open(filename, 'w') as file:
            if self.headers is not None:
                for header in self.headers:
                    file.write("{},".format(header))
            file.write("Affine L* scale,Affine a* scale,Affine b* scale,Affine constant,")
            for header in self.headers[self.fgnd_entry:]:
                file.write("{} corrected,".format(header))
            file.write("Difference,Corrected Difference,")
            file.write("\n")
            for cal in self.cal_values:
                # Sequentially write L*, a*, b* data
                colors = {0: "L*", 1: "a*", 2: "b*"}
                for color in colors:
                    for i in range(self.fgnd_entry - 1):
                        file.write("{},".format(cal["lead-in"][i]))
                    file.write("{},".format(colors[color]))
                    file.write("{},{},".format(cal["foreground"][color], cal["background"][color]))
                    for i in range(len(cal["medians"])):
                        file.write("{},".format(cal["medians"][i][color]))
                    for i in range(4):
                        file.write("{},".format(cal["affine"][color][i]))
                    file.write("{},{},".format(cal["corrected_foreground"][color],
                                               cal["corrected_background"][color]))
                    for i in range(len(cal["corrected_medians"])):
                        file.write("{},".format(cal["corrected_medians"][i][color]))
                    file.write("{},{},".format(cal["diffs"][color], cal["corr_diffs"][color]))
                    file.write("\n")


if __name__ == "__main__":
    import tkinter.filedialog
    colorCheckerCalAnalysis = ColorCheckerCalAnalysis(None)
    if len(sys.argv) > 1:
        for input_filename in sys.argv[1:]:
            colorCheckerCalAnalysis.load_from_colors_file(input_filename, nrefs=4)
            colorCheckerCalAnalysis.compute_colors(reference=CC_FOUR_COLOR_VALUES_V2_PATCH_WHITE)
            output_filename = input_filename[:-4] + "_analyzed.csv"
            if Path(output_filename).is_file():
                output_filename = tkinter.filedialog.asksaveasfilename(
                    title="Output file already exists. Choose output file")
            colorCheckerCalAnalysis.write_results(output_filename)
            colorCheckerCalAnalysis.clear_values()
    else:
        input_filename = tkinter.filedialog.askopenfilename()
        while input_filename:
            colorCheckerCalAnalysis.load_values(input_filename)
            colorCheckerCalAnalysis.compute_colors()
            output_filename = input_filename[:-4] + "_analyzed.csv"
            if Path(output_filename).is_file():
                output_filename = tkinter.filedialog.asksaveasfilename(title = "Output file already exists. Choose output file")
            colorCheckerCalAnalysis.write_results(output_filename)
            input_filename = tkinter.filedialog.askopenfilename()
            colorCheckerCalAnalysis.clear_values()
