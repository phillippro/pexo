from numpy import genfromtxt
from .settings import out_storage
import os


class PexoOut(object):
   """
   Read PEXO output file from the specified `path`.
   """
   def __init__(self):
      self._storage = out_storage


   def read(self, path):
      if not os.path.isfile(path):
         raise FileNotFoundError(f"PEXO output not found in the specified path: {path}")

      return genfromtxt(path, names=True)