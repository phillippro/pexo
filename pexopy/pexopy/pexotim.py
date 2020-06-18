from .uniquefilename import UniqueFilename
from .settings import tim_storage
import numbers
import os


class PexoTim(object):
   """
   PEXOpy .tim file handler.
   The constructor's `tim` argument is either a list/array of floats ot tuples of floats (1-part or 2-part JD[UTC]), or a path to a .tim file.

   If it's a path, the file is parsed and validated.

   If it's a list/array, it is validated and saved to a file.

   Either way, the class instance has both the path (<PexoTim.tim_path>) and the dictionary (<PexoTim.tim>).
   """
   def __init__(self, tim):
      self._storage = tim_storage

      if isinstance(tim, type(self)):
         self.data = tim.data
         self.path = tim.path
      elif isinstance(tim, str): # this is a path to a file
         self.data = self._parse_tim(tim)
         self.path = tim
      elif isinstance(tim, list) or 'array' in str(type(tim)): # this is a list or array
         self.data = tim
         self.path = self._generate_tim(tim)
      else:
         raise ValueError("PexoTim's argument `tim` should be a list of numbers, a list of tuples of numbers, or a path to a .tim file")

   
   def _parse_tim(self, tim_path):
      # read the .tim file into an array
      tim = []
      with open(tim_path) as f:
         for line in f:
            s = line.split()
            if len(s) == 1:
               tim.append(s[0])
            elif len(s) == 2:
               tim.append(tuple(s))
            else:
               raise ValueError(f"Error while parsing the .tim file: {tim_path}. Must only have one or two columns, 1-part or 2-part JD[UTC].")
         
      return tim
   

   def _generate_tim(self, tim):
      contents = ""
      try:
         for jd in tim:
            if isinstance(jd, numbers.Number):
               contents += f"{jd}\n"
            elif len(jd) == 2:
               contents += f"{jd[0]} {jd[1]}\n"
            else:
               raise TypeError()
      except TypeError:
         raise ValueError("`tim` argument should be a list of numbers, a list of tuples of numbers, or a path to a .tim file")
   
      filename = UniqueFilename(contents, append=".tim") # TODO : generate a unique name for this file
      path = os.path.join(self._storage, filename)
               
      with open(path, "w") as f:
         f.write(contents)

      return path