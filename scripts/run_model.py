import numpy as np
from numpy import pi, cos, sin, cosh, tanh
import matplotlib.pyplot as plt
from TwoDimBoussinesq import spectral_model as model

sech = lambda x: 1./np.cosh(x)

reload(model)

Fr =  0.1
Fr2 = Fr**2
Reb = 100

# the model object
m = model.Boussinesq2d(Lx=2.*np.pi,nx=128, tmax = 10, dt = 0.001, \
        use_fftw=True,ntd=4,Fr=Fr,use_filter=True,tsave=100,
        twrite=100,nu=Fr2/Reb,sig=1.e5,kf=25,
        save_snapshots=True,
        tsave_snapshots=1,
        ext_forc=True)

# run the model and plot some figs
plt.rcParams['image.cmap'] = 'RdBu'

plt.ion()

# forcing q
m0 = 3.
A0 = 1/2.
ep = pi/2.
zmin,zmax = pi-ep,pi+ep

# initial conditions
sig = 1.e-7
qi = sig*np.random.randn(m.nx,m.nx)
bi = sig*np.random.randn(m.nx,m.nx)

qi = qi
m.set_q(qi)
m.set_b(bi)

# forcing in the vorticity equation
#fq = ( (m0**3) / (Reb/Fr2))*cos(m0*m.y)
fq = cos(m0*m.z)
m.set_forcing(fq)

m.run()

#for snapshot in m.run_with_snapshots(tsnapstart=0, tsnapint=100*m.dt):
#
#    plt.clf()
#    p1 = plt.contourf(m.q,np.linspace(-30,30,15))
#    plt.clim([-30., 30.])
#    plt.title('t='+str(m.t))
#
#    plt.xticks([])
#    plt.yticks([])
#
#    plt.pause(0.001)
#
#    plt.draw()
#
#plt.show()
#plt.ion()
