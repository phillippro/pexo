from sgp4.api import Satrec
import numpy as np
from astropy.io import ascii
np.set_printoptions(precision=2)
def minmaxloc(num_list):
    return np.argmin(num_list), np.argmax(num_list)

f = open('hipparcos-tle.txt', 'r')
st = f.read().split('\n')
jds = np.genfromtxt("hipparcos_tle_jd.txt", dtype=None)
Njd = len(jds)
s = st[0:(Njd*2-1):2]
t = st[1:(Njd*2):2]
Nsim = 100000
jdsim = np.linspace(min(jds),max(jds),num=Nsim)
tmp=np.modf(jdsim)
jdi = tmp[1]
jdf = tmp[0]
es = np.zeros(Nsim)
rs = np.zeros((Nsim,3))
vs = np.zeros((Nsim,3))
for k in range(0,Nsim):
#    index = np.minmaxloc(abs(jds-jdsim[k]))[0]
    if k<(Nsim-1):
        index = np.where(jds-jdsim[k]>0)[0][0]-1
    else:
        index = Njd-1
    satellite = Satrec.twoline2rv(s[index], t[index])
    e, r, v = satellite.sgp4(jdi[k],jdf[k])
    es[k] = e
    rs[k][:] = r
    vs[k][:] = v

out = np.concatenate((jdsim,rs,vs),axis=1)
np.savetxt('hipparcos_state.txt', out, delimiter=' ')
