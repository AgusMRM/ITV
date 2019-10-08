import numpy as np
import matplotlib.pylab as plt

t=np.loadtxt('halosprop_2rvir_2.dat')
xbox=403.896
ybox=459.8882
zbox=440.9021  

x=408.2054  
y=457.7778
z=441.53868 

x=x-xbox+250
y=y-ybox+250
z=z-zbox+250
d   = np.sqrt((t[:,0]-x)**2+(t[:,1]-y)**2+(t[:,2]-z)**2)
gas = t[:,5]*1.81729961E-02
dm  = t[:,6]*(9.32880491E-02)
st  = t[:,8]
sf  = t[:,9]
spn = t[:,10]
hsml= t[:,12]
rvr = t[:,3]*1e-3    #para tenerlo en Mpc
rho_gs = t[:,11]*(1.81729961E-02)
rho_dm = dm/((4./3.)*np.pi*(2*rvr)**3)

up = np.where((gas/dm >= 0.1) & (hsml>0))
under = np.where((gas/dm < 0.1) & (hsml>0))

plt.plot(t[up,5])
plt.savefig(fname='/home/arodriguez/tesis/R1198/R_gdm.pdf',format='pdf')
