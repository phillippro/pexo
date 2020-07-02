import matplotlib.pyplot as plt
import matplotlib
import numpy as np


class PDF(object):
    def __init__(self, path, size=(200,200)):
        self.path = path
        self.size = size

    def _repr_html_(self):
        return '<iframe src={0} width={1[0]} height={1[1]}></iframe>'.format(self.path, self.size)

    def _repr_latex_(self):
        return r'\includegraphics[width=1.0\textwidth]{{{0}}}'.format(self.path)


class Table(object):
    def __init__(self, path, table=True, take_all=False, take_top=4, take_bottom=1, header=0, delimiter=" ", fontsize=10):
        self.path = path
        self.table = table
        self.take_all = take_all
        self.take_top = take_top
        self.take_bottom = take_bottom
        self.header = header
        self.delimiter = delimiter
        self.fontsize = fontsize


    def _read_table(self):
        with open(self.path, "r") as f:
            lines = f.readlines()
            if self.take_all:
                return lines

            n = len(lines)
            top_lines    = lines[0:self.take_top+1]       if self.take_top >= 0    else []
            bottom_lines = lines[n-self.take_bottom-1:-1] if self.take_bottom >= 0 else []
            skipped = n - self.take_top - self.take_bottom
            return top_lines + ["..."] + bottom_lines


    def _repr_html_(self):
        rows = self._read_table()
        html_rows = []
        
        for i in range(len(rows)):
            row = rows[i]
            element_type = "th" if self.header == i else "td"

            html_row_elements = ['<{} style="text-align: left; min-width:100px">{}</{}>'.format(element_type, x, element_type) for x in row.split(self.delimiter)]
            html_row = "".join(html_row_elements)
            html_rows.append(html_row)
        
        table_body = ['<tr>{}</tr>'.format(x) for x in html_rows]
        
        return '<table style="text-align: left; font-size: {}pt">{}</table>'.format(self.fontsize, "".join(table_body))


class Plot(object):
    def __init__(self, path, figsize=(12, 8), factors=(1, 1)):
        self.figsize = figsize
        self.data = np.genfromtxt(path, delimiter=" ", names=True)
        self.colours = ["#031926", "#95190C", "#FB8B24", "#5B9FAE", "#05C793", "#32746D"]
        self.nplots = 0
        self.factors = factors

        plt.style.use('seaborn-whitegrid')
        plt.rcParams.update({'font.size': 15})
        self.figure, self.axis = plt.subplots(1, 1, figsize=self.figsize)


    def _repr_html_(self):
        return ""


    def _get_column(self, col):
        if isinstance(col, str):
            return self.data[col]
        elif len(col) == 1:
            return self.data[col[0]]
        elif len(col) == 2:
            return self.data[col[0]] + self.data[col[1]]
        return []


    def _update_legend(self):
        handles, labels = self.axis.get_legend_handles_labels()
        if len(handles) > 0:
            self.axis.legend(handles, labels, ncol=5, loc="upper left", bbox_to_anchor=(0.0, 1.08), frameon=False)


    def add_plot(self, xcol, ycol, colour=None, yxdiff=False, legend="", xlabel="", ylabel="", factors=None, line="none", marker=None):
        x = self._get_column(xcol)
        y = self._get_column(ycol)

        if (line is None or line.lower() == "none") and (marker is None or marker.lower() == "none"):
            line = "-"
            marker = None

        if colour is None:
            colour = self.colours[ self.nplots % len(self.colours) ]
            self.nplots += 1
        
        if factors is None:
            factors = self.factors

        if yxdiff:
            y = y - x
        
        x *= factors[0]
        y *= factors[1]

        self.axis.plot(x, y, color=colour, lw=2, label=legend, alpha=0.7, ls=line, marker=marker)

        self.axis.set_xlabel(xlabel)
        self.axis.set_ylabel(ylabel)
        self._update_legend()

        return self
