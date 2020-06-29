import sys
import string
import astropy.units as u
import numpy as np
from astropy.table import Column
from astropy.table import Table
from astroquery.simbad import Simbad
from astropy.io import ascii
from astroquery.esasky import ESASky
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
star = sys.argv[1]
#star = 'GJ534'
#tables = Gaia.load_tables(only_names=True)
##get data from simbad
print('star: '+star)
#check all votable simbad fields: Simbad.list_votable_fields()
Simbad.add_votable_fields('pm', 'plx','rv_value','rvz_error','flux(V)', 'flux_error(V)')
result_table = Simbad.query_object(star)
ra0 = [float(x) for x in result_table['RA'][0].split()]
ra1 = (ra0[0]+ra0[1]/60+ra0[2]/3600)*360/24#deg                                 
print('ra1=%f3'%ra1)

dec0 = [float(x) for x in result_table['DEC'][0].split()]
dec1 = dec0[0]+dec0[1]/60+dec0[2]/3600#deg                                                                            
print('dec1=%f3'%dec1)

pmra = result_table['PMRA'][0]
pmdec = result_table['PMDEC'][0]
##convert to Gaia epoch                                                                                               
#dec = dec1 + 15.5*pmdec/3600000#deg                                                                                   
dec = dec1
#ra = ra1 + (15.5*pmra/3600000)/np.cos((dec+dec1)*np.pi/360)#deg                                                      
ra = ra1 
rv = result_table['RV_VALUE'][0]
erv = result_table['RVZ_ERROR'][0]
plx = result_table['PLX_VALUE'][0]

##get data from Gaia
##get source id
r = Simbad.query_objectids(star)
g = [x for x in r if 'Gaia DR2' in str(x)]
h = [x for x in r if 'HIP' in str(x)]

if len(h)==0:
    sys.exit("No Hipparcos source found to match the target!")

h1 = h[0][0].split()[-1]
hid = int(h1)
print('HIP'+str(hid))

if len(g)==0:
    query = "SELECT *  FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',"+str(ra)+","+str(dec)+","+"0.1))=1 AND gaiadr2.gaia_source.parallax/"+str(plx)+"<1.1 AND gaiadr2.gaia_source.parallax/"+str(plx)+">0.9;"
else:
    g1 = str(g).split()[-1]
    gid = int(g1.split(']')[0])
    query = "select * from gaiadr2.gaia_source where source_id="+str(gid)

#print(query)
job = Gaia.launch_job(query=query)
out = job.get_results()
N = np.shape(out)[0]
hh = Column(name='HIP', data=[hid]*N)
v = Column(name='RV', data=[rv]*N)
ev = Column(name='eRV', data=[erv]*N)
out.add_columns([hh,v,ev])
ascii.write(out, 'gaia_hip.csv',format='csv',overwrite=True)
