class PDF(object):
  def __init__(self, path, size=(200,200)):
    self.path = path
    self.size = size

  def _repr_html_(self):
    return '<iframe src={0} width={1[0]} height={1[1]}></iframe>'.format(self.path, self.size)

  def _repr_latex_(self):
    return r'\includegraphics[width=1.0\textwidth]{{{0}}}'.format(self.path)


class Txt(object):
  def __init__(self, path, table=True, take_top=5, take_bottom=3, header=0, delimiter=" "):
    self.path = path
    self.table = table
    self.take_top = take_top
    self.take_bottom = take_bottom
    self.header = header
    self.delimiter = delimiter
      
  def _read_table(self):
    with open(self.path, "r") as f:
      lines = f.readlines()
      n = len(lines)
        
      top_lines    = lines[0:self.take_top]       if self.take_top >= 0    else []
      bottom_lines = lines[n-self.take_bottom:-1] if self.take_bottom >= 0 else []
        
      return top_lines + ["..."] + bottom_lines

  def _repr_html_(self):
    rows = self._read_table()
    html_rows = []
    
    for i in range(len(rows)):
      row = rows[i]
      element_type = "th" if self.header == i else "td"
      html_row = ['<{}  style="text-align: left;">{}</{}>'.format(element_type, x, element_type) for x in row.split(self.delimiter)]
      html_rows.append("".join(html_row))
    
    table_body = ['<tr>{}</tr>'.format(x) for x in html_rows]
    
    return '<table style="text-align: left;">{}</table>'.format("".join(table_body))
