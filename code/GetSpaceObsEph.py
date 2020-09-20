import sys
import numpy as np
from astropy.io import ascii
from astroquery.jplhorizons import Horizons

#script = sys.argv[0]
#code = sys.argv[1]
#filename = sys.argv[2]

code = -48
filename = 'test.tim'
epochs = np.loadtxt(filename)
obj = Horizons(id=code, location='@ssb',id_type='id',epochs=epochs)
print(obj.uri)
vec = obj.vectors()

ascii.write(vec[['x','y','z','vx','vy','vz']], '../observatories/space.csv', format='csv', fast_writer=False)
