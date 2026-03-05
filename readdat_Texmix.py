import numpy as np                          
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

##Iteration for each method
dat1 = ascii.read('tkin_iter1.dat')
ld1  = len(dat1)
dat2 = ascii.read('tbg_iter1.dat')
ld2  = len(dat2)
dat3 = ascii.read('mix_iter1.dat')
ld3  = len(dat3)

r1 = dat1['col1'][0:ld1]
r2 = dat2['col1'][0:ld2]
r3 = dat3['col1'][0:ld3]

d1 = dat1['col2'][0:ld1]
d2 = dat2['col2'][0:ld2]
d3 = dat3['col2'][0:ld3]

fig1 = figure(1, figsize=(8, 4))#, layout='constrained')
plt.grid(True, linestyle='--')
plt.plot(r2[0:ld2], d2[0:ld2],  linewidth=2.0, label='$T_{bg}$')
plt.plot(r1[0:ld1], d1[0:ld1], linewidth=2.0, label='$T_{kin}$')
plt.plot(r3[0:ld3], d3[0:ld3],  linewidth=2.0, label='$T_{bg};T_{kin}$')
plt.yscale('log', nonpositive='clip')
plt.ylabel('$R_{diff}$', fontsize=10)
plt.xlabel('Nr. Iterations', fontsize=10)
plt.ylim([10e-7, 0])
plt.legend()
plt.savefig('comp_iter_tkin_hco_var1.png')
fig1.show()

##Tex for each method
#number of transistion -hco++nl = 20
nl = 20
nj = 16

dat_tbg  = ascii.read('tbg1.out')
dat_tkin = ascii.read('tkin1.out')
dat_mix  =  ascii.read('mix1.out')

lb = len(dat_tbg)
lt = len(dat_tkin)
lm = len(dat_mix)

db = {}
dk = {}
dm = {}
dm1={}
dm2={}
for j in range(nj):
  db[j] = dat_tbg['col3'][j:lb:nl]
  dk[j] = dat_tkin['col3'][j:lt:nl]
  dm[j] = dat_mix['col5'][j:lm:nl]
  dm1[j] = dat_mix['col3'][j:lm:nl]
  dm2[j] = dat_mix['col4'][j:lm:nl]
  
n1 = int(lb/nl)
n2 = int(lt/nl)
n3 = int(lm/nl)
xb = np.arange(0,n1,1)
xt = np.arange(0,n2,1)
xm = np.arange(0,n3,1)

#plt.show()
fig2 = figure(2, figsize=(8, 4))
ax = plt.gca()
plt.plot(xb,db[0], label="$J=1-0$")
plt.plot(xb,db[4], label="$J=5-4$")
plt.plot(xb,db[6], label="$J=7-6$")
plt.plot(xb,db[11], label="$J=12-11$")
plt.grid(True, linestyle='--')
plt.legend()
plt.ylabel('$T_{ex}$ [K]', fontsize=10)
xl = 'Nr. Iterations'
plt.xlabel(xl, fontsize=10)
#plt.title("J=3-2; f={:.2f}; Optical depth={:.3E}".format(f[2][0], df[2][ldf-1]), fontsize=10)
fig2.show()

fig3 = figure(3, figsize=(8, 4))
ax = plt.gca()
plt.plot(xt,dk[0], label="$J=1-0$")
plt.plot(xt,dk[4], label="$J=5-4$")
plt.plot(xt,dk[6], label="$J=7-6$")
plt.plot(xt,dk[11], label="$J=12-11$")
plt.grid(True, linestyle='--')
plt.legend()
plt.ylabel('$T_{ex}$ [K]', fontsize=10)
xl = 'Nr. Iterations'
plt.xlabel(xl, fontsize=10)
#plt.title("J=3-2; f={:.2f}; Optical depth={:.3E}".format(f[2][0], df[2][ldf-1]), fontsize=10)
fig3.show()

fig3, axs = plt.subplots(2, 2)
ax = plt.gca()
axs[0,0].plot(xm,dm1[0], label="$t_{bg}:J=1-0$")
axs[0,0].plot(xm,dm2[0], label="$t_{kin}:J=1-0$")
axs[0,0].plot(xm,dm[0], label="$t_{bg}; t_{kin}:J=1-0$")
axs[0,0].grid(True, linestyle='--')
axs[0,0].legend()
axs[0,0].set_ylabel('$T_{ex}$ [K]', fontsize=8)
xl = 'Nr. Iterations'
axs[0,0].set_xlabel(xl, fontsize=8)
xl = 'Nr. Iterations'
axs[0,1].set_xlabel(xl, fontsize=8)
axs[0,1].plot(xm,dm1[4], label="$t_{bg}:J=5-4$")
axs[0,1].plot(xm,dm2[4], label="$t_{kin}:J=5-4$")
axs[0,1].plot(xm,dm[4], label="$t_{bg}; t_{kin}:J=5-4$")
axs[0,1].grid(True, linestyle='--')
axs[0,1].legend()
axs[0,1].set_ylabel('$T_{ex}$ [K]', fontsize=8)
xl = 'Nr. Iterations'
axs[0,1].set_xlabel(xl, fontsize=8)
axs[1,0].plot(xm,dm1[6], label="$t_{bg}:J=7-6$")
axs[1,0].plot(xm,dm2[6], label="$t_{kin}:J=7-6$")
axs[1,0].plot(xm,dm[6], label="$t_{bg}; t_{kin}:J=7-6$")
axs[1,0].grid(True, linestyle='--')
axs[1,0].legend()
axs[1,0].set_ylabel('$T_{ex}$ [K]', fontsize=8)
xl = 'Nr. Iterations'
axs[1,0].set_xlabel(xl, fontsize=8)
#axs[1,1] = plt.gca()
axs[1,1].plot(xm,dm1[11], label="$t_{bg}:J=12-11$")
axs[1,1].plot(xm,dm2[11], label="$t_{kin}:J=12-11$")
axs[1,1].plot(xm,dm[11], label="$t_{bg}; t_{kin}:J=12-11$")
axs[1,1].grid(True, linestyle='--')
axs[1,1].legend()
axs[1,1].set_ylabel('$T_{ex}$ [K]', fontsize=8)
xl = 'Nr. Iterations'
axs[1,1].set_xlabel(xl, fontsize=8)
fig3.tight_layout(rect=[0, 0, 1, 0.95])
#plt.title("J=3-2; f={:.2f}; Optical depth={:.3E}".format(f[2][0], df[2][ldf-1]), fontsize=10)
fig3.show()
