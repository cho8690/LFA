from colorCheckerImage import ColorCheckerImage
from colorCheckerAffine import ColorCheckerAffine
from colorCheckerConstants import *
import tkinter.filedialog
import csv
import math

_author_ = "A.J. Aranyosi"
_project_ = "colorChecker"
_file_ = "colorCalc.py"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>
"""


class colorCalc(object):
    """
    Class for calculating color values and correcting them using a color checker
    """

    def __init__(self, parent, image=None, affine=None, ab_affine=None):
        """
        Constructor for colorCalc
        """
        self.parent = parent
        self.image = image
        if image is None:
            self.image = ColorCheckerImage(self)
        self.affine = affine
        if affine is None:
            self.affine = ColorCheckerAffine(self)
        self.ab_affine = ab_affine
        if ab_affine is None:
            self.ab_affine = ColorCheckerAffine(self)
        self.cal_values = []
        self.foregrounds = []
        self.backgrounds = []
        self.differences = []
        self.corrected_data = None
        self.corrected_cal_values = []
        self.corrected_foregrounds = []
        self.corrected_backgrounds = []
        self.corrected_differences = []
        self.headers = None
        self.filename = None

    def read_csv(self, header_lines=1, filename=None):
        values = [0, 0, 0]
        if filename is None:
            self.filename = tkinter.filedialog.askopenfilename()
            if self.filename[-3:].lower() != "csv":
                print("Error: this program needs a CSV file as input.")
                return False
        with open(self.filename, newline='') as csvfile:
            file_reader = csv.reader(csvfile, delimiter=',')
            for _ in range(header_lines):
                self.headers = next(file_reader, None)
            for row in file_reader:
                values[(file_reader.line_num + 1) % 3] = float(row[-1])
                if 3 < file_reader.line_num <= 73 and file_reader.line_num % 3 == 1:
                    self.cal_values.append(values[:])
                if 73 < file_reader.line_num and file_reader.line_num % 3 == 1:
                    if file_reader.line_num % 6 == 1:
                        self.backgrounds.append(values[:])
                    else:
                        self.foregrounds.append(values[:])

    def set_foregrounds(self, values):
        self.foregrounds = values[:]

    def set_backgrounds(self, values):
        self.backgrounds = values[:]

    def set_affine(self, affine):
        self.affine = affine

    def set_cal_values(self, values):
        self.cal_values = values

    def compute_affine(self, reference=CC_REF_COLORS):
        if self.cal_values is None or len(self.cal_values) == 0:
            return False
        self.affine.get_affine(data=self.cal_values, reference=reference)

    def apply_affine(self, data=None):
        if self.cal_values is not None and len(self.cal_values) > 0:
            self.corrected_cal_values = self.affine.apply_affine(data=self.cal_values).tolist()
        if self.foregrounds is not None and len(self.foregrounds) > 0:
            self.corrected_foregrounds = self.affine.apply_affine(data=self.foregrounds).tolist()
        if self.backgrounds is not None and len(self.backgrounds) > 0:
            self.corrected_backgrounds = self.affine.apply_affine(data=self.backgrounds).tolist()
        if data is not None and len(data) > 0:
            self.corrected_data = self.affine.apply_affine(data=data).tolist()

    def measure_colors(self):
        self.differences = [[x - y for (x, y) in zip(t, r)] for (t, r) in zip(self.foregrounds, self.backgrounds)]
        self.corrected_differences = [[x - y for (x, y) in zip(t, r)] for (t, r) in
                                      zip(self.corrected_foregrounds, self.corrected_backgrounds)]

    def write_results(self, output_filename=None):
        if output_filename is None:
            output_filename = self.filename[:-4] + "_analyzed" + self.filename[-4:]
        with open(output_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['Color', 'L*', 'a*', 'b*', 'Corrected L*', 'Corrected a*', 'Corrected b*'])
            for i in range(len(self.cal_values)):
                writer.writerow([CC_COLOR_NAMES[i]] + self.cal_values[i] + self.corrected_cal_values[i])
            for i in range(4):
                writer.writerow(["Affine row {}".format(i)] + self.affine.a_matrix[i, :].tolist())
            writer.writerow(["Foreground Analysis"])
            for i in range(len(self.foregrounds)):
                writer.writerow(["Foreground"] + self.foregrounds[i] + self.corrected_foregrounds[i])
                writer.writerow(["Background"] + self.backgrounds[i] + self.corrected_backgrounds[i])
                writer.writerow(["Difference"] + self.differences[i] + self.corrected_differences[i])


if __name__ == "__main__":
    colorCalc = colorCalc(None)
    colorCalc.read_csv()
    colorCalc.compute_affine()
    colorCalc.apply_affine()
    colorCalc.write_results()
