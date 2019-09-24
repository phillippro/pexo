import sys
from astropy.io import ascii
import numpy as np
from astroquery.jplhorizons import Horizons
def main():
    script = sys.argv[0]
    code = sys.argv[1]
    filename = sys.argv[2]
    epochs = np.loadtxt(filename, delimiter=',')
    obj = Horizons(id=code, location='@ssb',id_type='id',epochs=epochs)
    print(obj.uri)
    vec = obj.vectors()
    ascii.write(vec[['x','y','z','vx','vy','vz']], '../observatories/space.csv', format='csv', fast_writer=False)

if __name__ == '__main__':
   main()
