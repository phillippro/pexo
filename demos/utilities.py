
import os
from subprocess import Popen, PIPE, call, check_output

import matplotlib.pyplot as plt
import matplotlib

from numpy import genfromtxt, array
from rpy2.robjects import r
from rpy2.rinterface import NARealType


class Struct(object):
    def __init__(self, dictionary):
        """
        Simple helper class to "convert" dictionaries into property containers.

        Usge:
        In [1]: d = {'foo': 34, 'bar': 19}
        In [2]: s = Struct(d)
        In [3]: s
                <__main__.Struct at 0x7fda49cf11c0>
        In [4]: s.foo
                34
        In [5]: s.bar
                19
        In [6]: s.bar = 23
        In [7]: s.bar
                23
        In [8]: s.woohooo = 23
        In [9]: s.woohooo
                23
        """
        self.__dict__.update(dictionary)

        for k, v in dictionary.items():
            if isinstance(v, dict):
                self.__dict__[k] = Struct(v)


    def __setattr__(self, name, value):
        self.__dict__[name] = value


    @property
    def dictionary(self):
        return self.__dict__


class PDF(object):
    def __init__(self, path, size=(200,200)):
        self.path = path
        self.size = size

    def _repr_html_(self):
        return '<iframe src={0} width={1[0]} height={1[1]}></iframe>'.format(self.path, self.size)

    def _repr_latex_(self):
        return r'\includegraphics[width=1.0\textwidth]{{{0}}}'.format(self.path)


class Table(object):
    def __init__(self, path=None, lines=None, table=True, take_all=False, take_top=4, take_bottom=1, header=0, left_header=-1, delimiter=" ", fontsize=10, digits=3):
        if lines is not None:
            self.lines = lines
        else:
            self.path = path
            with open(self.path, "r") as f:
                self.lines = f.readlines()

        self.table = table
        self.take_all = take_all
        self.take_top = take_top
        self.take_bottom = take_bottom
        self.header = header
        self.left_header = left_header
        self.delimiter = delimiter
        self.fontsize = fontsize
        self.digits = digits


    def _read_table(self):
        if self.take_all:
            return self.lines

        n = len(self.lines)
        top_lines    = self.lines[0:self.take_top+1]       if self.take_top >= 0    else []
        bottom_lines = self.lines[n-self.take_bottom-1:-1] if self.take_bottom >= 0 else []
        skipped = n - self.take_top - self.take_bottom
        return top_lines + ["..."] + bottom_lines


    def _repr_html_(self):
        def _is_float(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        rows = self._read_table()
        html_rows = []
        
        for i in range(len(rows)):
            row = rows[i]
            element_type = "th" if self.header == i else "td"

            # apply the float precision
            row_items = row.replace("\n", "").split(self.delimiter)
            for i, item in enumerate(row_items):
                if _is_float(item):
                    split = item.split(".")
                    if len(split) == 2:
                        split[1] = split[1][0:self.digits]
                    row_items[i] = ".".join(split)

            html_row_elements = ['<{} style="text-align: right; padding: 7px;">{}</{}>'.format(element_type, x, element_type) for x in row_items]
            if self.left_header > -1:
                html_row_elements[self.left_header] = html_row_elements[self.left_header].replace('style="', 'style="font-weight: bold; ')
            html_row = "".join(html_row_elements)
            html_rows.append(html_row)
        
        table_body = ['<tr>{}</tr>'.format(x) for x in html_rows]
        
        return '<table style="text-align: left; font-size: {}pt">{}</table>'.format(self.fontsize, "".join(table_body))


class Plot(object):
    def __init__(self, path=None, data=None, figsize=(12, 8), factors=(1, 1), plot_instance=None):
        if data is None:
            self.data = genfromtxt(path, delimiter=" ", names=True)
        else:
            self.data = data

        self.figsize = figsize
        self.colours = ["#031926", "#95190C", "#FB8B24", "#5B9FAE", "#05C793", "#32746D"]
        self.nplots = 0
        self.factors = factors

        plt.style.use('seaborn-whitegrid')
        plt.rcParams.update({'font.size': 15})
        if plot_instance is None:
            self.figure, self.axis = plt.subplots(1, 1, figsize=self.figsize)
        else:
            self.figure, self.axis = plot_instance.figure, plot_instance.axis


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

        x, y = (array(t) for t in zip(*sorted(zip(x, y))))

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


class EmulationOutput(object):
    """
    Read PEXO emulation output file (.txt) from the specified `path`.
    """
    def __init__(self, path):
        if not os.path.isfile(path):
            errormessage = "PEXO output not found in the specified path: {}".format(path)
            raise FileNotFoundError(errormessage)
        
        self.path = path
        self.contents = genfromtxt(path, names=True)


class FitOutput(object):
    """
    Read PEXO fit output file (.Robj) from the specified `path`.
    """
    def __init__(self, path):
        if not os.path.isfile(path):
            errormessage = "PEXO output not found in the specified path: {}".format(path)
            raise FileNotFoundError(errormessage)

        self.path = path
        self.contents = r.load(path)


    @property
    def parstat(self):
        r_obj = r['ParStat']
        obj  = array(r_obj).T
        colnames = r.colnames(r_obj)

        properties = ["xopt", "x1per", "x99per", "x10per", "x90per", "xminus", "xplus", "mode", "mean", "sd", "skewness", "kurtosis"]
        values = [dict(zip(properties, obj[i])) for i in range(len(colnames))]

        matrix = [ [""] + properties, ]
        for i, param in enumerate(colnames):
            temp = [param,]
            for prop in properties:
                temp.append(values[i][prop])
            matrix.append(temp)

        return matrix


    @property
    def data(self):
        def _fix_nones_list(column_list):
            return [None if isinstance(x, NARealType) else x for x in column_list]

        r_obj     = r['Data']
        colnames  = r.colnames(r_obj)
        table = [_fix_nones_list(x) for x in r_obj]

        return Struct(dict(zip(colnames, table)))


    @property
    def model(self):
        def _fix_nones_list(column_list):
            return [None if isinstance(x, NARealType) else x for x in column_list]

        r_obj     = r['model']
        colnames  = r.colnames(r_obj)
        table = [_fix_nones_list(x) for x in r_obj]

        return Struct(dict(zip(colnames, table)))
