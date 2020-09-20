import sys
import string
import numpy as np

import astropy.units as u
from astropy.table import Column, Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord

from astroquery.simbad import Simbad
from astroquery.esasky import ESASky
from astroquery.gaia import Gaia


star = sys.argv[1]
#star = 'HD128620'
#tables = Gaia.load_tables(only_names=True)
##get data from simbad
#print('star: ' + star)
#check all votable simbad fields: Simbad.list_votable_fields()

Simbad.add_votable_fields('pm', 'plx','rv_value','rvz_error','flux(V)', 'flux_error(V)')
result_table = Simbad.query_object(star)

ra0 = [float(x) for x in result_table['RA'][0].split()]
ra1 = (ra0[0] + ra0[1]/60 + ra0[2]/3600) * 360/24 #deg

dec0 = [float(x) for x in result_table['DEC'][0].split()]
dec1 = dec0[0] + dec0[1]/60 + dec0[2]/3600 #deg

pmra  = result_table['PMRA'][0]
pmdec = result_table['PMDEC'][0]

##convert to Gaia epoch
#dec = dec1 + 15.5*pmdec/3600000 #deg
#ra = ra1 + (15.5*pmra/3600000)/np.cos((dec + dec1)*np.pi/360) #deg
dec = dec1
ra = ra1

rv  = result_table['RV_VALUE'][0]
erv = result_table['RVZ_ERROR'][0]
plx = result_table['PLX_VALUE'][0]

##get data from Gaia
##get source id
r = Simbad.query_objectids(star)
g = [x for x in r if 'Gaia DR2' in str(x)]
h = [x for x in r if 'HIP' in str(x)]

if len(h) == 0:
    sys.exit("No Hipparcos source found to match the target!")

h1 = h[0][0].split()[-1]
hid = int(h1)

if len(g) == 0:
    query = "SELECT * FROM gaiadr2.gaia_source " + \
        "WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{},{},0.1))=1 ".format(ra, dec) + \
        "AND gaiadr2.gaia_source.parallax/{}<1.1 ".format(plx) + \
        "AND gaiadr2.gaia_source.parallax/{}>0.9;".format(plx)
else:
    g1 = str(g).split()[-1]
    gid = int(g1.split(']')[0])
    query = "select * from gaiadr2.gaia_source where source_id={}".format(gid)

job = Gaia.launch_job(query=query)
out = job.get_results()
N = np.shape(out)[0]
if N == 0:
    N = 1

hh = Column(name='HIP', data=[hid]*N)

#ra = Column(name='ra', data=[ra1]*N)
#dec = Column(name='dec', data=[dec1]*N)
#pmra = Column(name='pmra', data=[pmra]*N)
#pmdec = Column(name='pmdec', data=[pmdec]*N)
#plx = Column(name='parallax', data=[plx]*N)
#rv = Column(name='rv', data=[rv]*N)
#erv = Column(name='erv', data=[erv]*N)

fout = star + '_gaia_hip.csv'

v  = Column(name='RV',  data=[rv]*N)
ev = Column(name='eRV', data=[erv]*N)
if len(out) > 0:
    out.add_columns([hh, v, ev])
    ascii.write(out, fout, format='csv', overwrite=True)
else:
    out = np.array([hid, ra, dec, plx, pmra, pmdec, rv, erv]).transpose()
    ascii.write(out, fout, 
        format='csv',
        overwrite=True,
        names=['HIP','ra','dec','parallax','pmra','pmdec','radial_velocity','radial_velocity_error']
    )

