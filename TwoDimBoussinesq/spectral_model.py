from __future__ import division
import numpy as np
from numpy import pi, exp, sqrt, cos, sin

import os
import h5py

try:
    import pyfftw
    pyfftw.interfaces.cache.enable()
except ImportError:
    pass

class Boussinesq2d(object):

    """ 2D Boussinesq code in vorticity-streamfunction/buoyancy formulation """

    def __init__(
            # grid parameters
            self,
            nx = 128,
            nz=None,
            Lx=2*pi,                    
            Lz=None,                    
            # physical parameters
            nu = 0.,
            Fr   = .1,
            R = 0.99,
            sig = 1.,
            ext_forc=False,
            kf = 4.,
            # timestepping parameters
            dt=.0025,               # numerical timestep
            twrite=100,             # interval for cfl and ke printout (in timesteps)
            tmax=100.,              # total time of integration
            filt=True,              # spectral filter flag
            use_fftw=True,
            ntd = 1,                # number of threads for fftw
            # filter or dealiasing (False for 2/3 rule)
            use_filter=True,
            # saving parameters
            tsave=100,              # interval to save (in timesteps)
            save_snapshots=True,
            overwrite=True,
            tsave_snapshots=100,
            path = 'output/'):

        if nz is None: nz = nx
        if Lz is None: Lz = Lx
       
        # initialize parameters

        # domain
        self.nx = nx
        self.nz = nz
        self.Lx = Lx
        self.Lz = Lz
        self.dx = Lx/nx
        self.dz = Lz/nz
        self.x = np.arange(0.,Lx,self.dx)
        self.z = np.arange(0.,Lz,self.dz)
        self.x,self.z = np.meshgrid(self.x,self.z)

        # physical
        self.nu = nu
        self.Fr   = Fr
        self.Fr2 = Fr**2
        self.sig = 1.e6
        self.R = R
        self.sig = sig
        self.kf = kf
        
        # ext forcing
        self.ext_forc = ext_forc

        # time related variables
        self.nmax = int(np.ceil(tmax/dt))        
        self.dt = dt
        self.twrite = twrite
        self.tmax = tmax
        self.t = 0.
        self.ndt = 0
        self.tsave = tsave
        self.tsave_snapshots=tsave_snapshots        
        self.nsave_max = int(np.ceil(self.nmax/tsave))
        self.save_snapshots = save_snapshots
        self.nsave = 0 
        self.overwrite=True
        self.snapshots_path = path

        # fourier settings
        self._init_kxky()
        self.kappa2 = self.k**2 + self.l**2
        self.kappa = sqrt(self.kappa2)
        self.fnz = self.kappa2 != 0
        self.kappa2i = np.zeros_like(self.kappa2)   # inversion not defined at kappa=0
        self.kappa2i[self.fnz] = self.kappa2[self.fnz]**-1
         
        # exponential filter or dealising
        self.use_filter = use_filter        
        self._init_filter()

        # fftw
        self.use_fftw = use_fftw
        self.ntd = ntd
        
        # allocate variables
        self._allocate_variables()

        # DFT
        self._initialize_fft()

        # initialize time-stepper
        self._init_rk3w()

        # initialize fno
        if self.save_snapshots:
            self._init_save_snapshots(self.snapshots_path)

        # setup file
        self._save_setup()

    def run(self):
        """ step forward until tmax """

        while(self.t < self.tmax): 
           
            if not(self.ext_forc):
                self._random_forcing()

            self._stepforward()

            if (self.ndt%self.twrite == 0.):
                self._printout()
            if self.save_snapshots:
                self._save_snapshots()
                
            self.t += self.dt
            self.ndt += 1

        self._save_diagnostics()

    def _stepforward(self):

        """ march the system forward using a RK3W-theta scheme """
 
        self.nl1h = -self.jacobian() + self.kj*self.bh/self.Fr2 + self.fh
        self.nl1h_b = -self.jacobian_b() - self.kj*self.ph 
        self.qh = (self.L1*self.qh + self.c1*self.dt*self.nl1h).copy()
        self.qh = self.filt*self.qh
        self.bh = (self.L1*self.bh + self.c1*self.dt*self.nl1h_b).copy()
        self.bh = self.filt*self.bh

        self.nl2h = self.nl1h.copy()
        self.nl1h = -self.jacobian() + self.kj*self.bh/self.Fr2 + self.fh
        self.nl2h_b = self.nl1h_b.copy()
        self.nl1h_b = -self.jacobian_b() - self.kj*self.ph
        self.qh = (self.L2*self.qh + self.c2*self.dt*self.nl1h +\
                self.d1*self.dt*self.nl2h).copy()
        self.qh = self.filt*self.qh
        self.bh = (self.L2*self.bh + self.c2*self.dt*self.nl1h_b +\
                self.d1*self.dt*self.nl2h_b).copy()
        self.bh = self.filt*self.bh

        self.nl2h = self.nl1h.copy()
        self.nl1h = -self.jacobian()  + self.kj*self.bh/self.Fr2 + self.fh
        self.nl2h_b = self.nl1h_b.copy()
        self.nl1h_b = -self.jacobian_b() - self.kj*self.ph
        self.qh = (self.L3*self.qh + self.c3*self.dt*self.nl1h +\
                self.d2*self.dt*self.nl2h).copy()
        self.qh = self.filt*self.qh
        self.bh = (self.L3*self.bh + self.c3*self.dt*self.nl1h_b +\
                self.d2*self.dt*self.nl2h_b).copy()
        self.bh = self.filt*self.bh

    def _allocate_variables(self):
        """ Allocate variables in memory """

        dtype_real = np.dtype('float64')            
        dtype_cplx = np.dtype('complex128')
        shape_real = (self.nz, self.nx)
        shape_cplx = (self.nz, self.nx/2+1)
            
        # vorticity
        self.q  = np.zeros(shape_real, dtype_real)
        self.qh = np.zeros(shape_cplx, dtype_cplx)
        # streamfunction
        self.p  = np.zeros(shape_real, dtype_real)
        self.ph = np.zeros(shape_cplx, dtype_cplx)
        # buoyancy
        self.b = np.zeros(shape_real, dtype_real)
        self.bh = np.zeros(shape_cplx, dtype_cplx)
        # velocity
        self.u = np.zeros(shape_real, dtype_real)
        self.v = np.zeros(shape_real, dtype_real)
        # nonlinear-terms
        self.nl1h = np.zeros(shape_cplx, dtype_cplx)
        self.nl2h = np.zeros(shape_cplx, dtype_cplx)
        self.nl1h_b = np.zeros(shape_cplx, dtype_cplx)
        self.nl2h_b = np.zeros(shape_cplx, dtype_cplx)

        # amplitude of forcing
        self.A  = np.zeros(shape_cplx, dtype_real)

    def _initialize_fft(self):
        # set up fft functions for use later
        if self.use_fftw: 
            self.fft2 = (lambda x :
                    pyfftw.interfaces.numpy_fft.rfft2(x, threads=self.ntd,\
                            planner_effort='FFTW_ESTIMATE'))
            self.ifft2 = (lambda x :
                    pyfftw.interfaces.numpy_fft.irfft2(x, threads=self.ntd,\
                            planner_effort='FFTW_ESTIMATE'))
        else:
            self.fft2 =  (lambda x : np.fft.rfft2(x))
            self.ifft2 = (lambda x : np.fft.irfft2(x))

    def _init_kxky(self):
        """ Calculate wavenumbers """

        self.dl = 2.*pi/self.Lz
        self.dk = 2.*pi/self.Lx
        self.ll = self.dl*np.append( np.arange(0.,self.nz/2),
                np.arange(-self.nz/2,0.) )
        self.kk = self.dk*np.arange(0.,self.nx/2+1)
        self.k,self.l = np.meshgrid(self.kk,self.ll)
        self.kj = 1j*self.k
        self.lj = 1j*self.l

    def _invert(self):
        """ Compute streamfunction from vorticity """
        self.ph = -self.kappa2i*self.qh

    def set_q(self,q):
        """ Initialize vorticity """ 

        self.q = q
        self.qh = self.fft2(self.q)
        self._invert()
        self.ph = self.filt * self.ph

    def set_b(self,b):
        """ Initialize buoyancy """ 

        self.b = b
        self.bh = self.fft2(self.b)
        self.bh = self.filt * self.bh

    def set_forcing(self,f):
        """ Initialize forcing """
        if self.ext_forc:
            self.f = f
            self.fh = self.fft2(self.f)
            self.fh = self.filt * self.fh
        else:
            self.fhp = 0.
            self.A = self.sig*np.exp(-((self.kappa-self.kf)/2)**4)
            self.A[:,0] = 0.
            self.A[0,:] = 0.
            self._random_forcing()

    def _random_forcing(self):
        """ Random forcing """
        th = 2*pi*np.random.rand(self.nx,self.nz/2+1)
        self.fh = self.A*(1.-self.R)*np.exp(1j*th) + self.R*self.fhp
        self.fhp = self.fh.copy() 

    def jacobian(self):
        """ Compute the Jacobian in conservative form """

        self._invert()
        self.ph = self.filt*self.ph
        self.q = self.ifft2(self.qh)
        self.u = self.ifft2(-self.lj*self.ph) 
        self.v = self.ifft2( self.kj*self.ph)

        jach = self.kj*self.fft2(self.u*self.q) +\
                self.lj*self.fft2(self.v*self.q)

        return jach

    def jacobian_b(self):
        """ Compute the Jacobian between psi and b in conservative form """

        self.b = self.ifft2(self.bh)
        jach = self.kj*self.fft2(self.u*self.b) +\
                self.lj*self.fft2(self.v*self.b)

        return jach

    def _printout(self):
        """ Print model status """
        if (self.ndt%self.twrite == 0.):        
            self.ke = self._calc_ke()
            self.pe = self._calc_pe()
            self.ens = self._calc_ens()
            self.ani = self._calc_anisotropy()
            self.cfl = self._calc_cfl()
            print "t= %e, cfl= %e, ke= %e, pe= %e, ens= %e, ani=%e" %(self.t, self.cfl, 
                self.ke, self.pe, self.ens, self.ani)
            assert self.cfl<1., "CFL condition violated"
 
    def _init_save_snapshots(self,path):

        self.fno = path

        if not os.path.isdir(self.fno):
            os.makedirs(self.fno)
            os.makedirs(self.fno+"/snapshots/")

    def _file_exist(self, fno):
        if os.path.exists(fno):
            if self.overwrite:
                os.remove(fno)
            else: raise IOError("File exists: {0}".format(fno))

    def _save_setup(self,):

        """Save setup  """

        fno = self.fno + 'setup.h5'

        self._file_exist(fno)

        h5file = h5py.File(fno, 'w')
        
        h5file.create_dataset("grid/nx", data=(self.nx),dtype=int)
        h5file.create_dataset("grid/nz", data=(self.nz),dtype=int)
        h5file.create_dataset("grid/x", data=(self.x))
        h5file.create_dataset("grid/z", data=(self.z))
        h5file.create_dataset("grid/kappa", data=self.kappa)

        h5file.close()

    def _save_snapshots(self, fields=['q','b','u','v']):

        """ Save snapshots of fields """

        fno = self.fno + '/snapshots/{:015.0f}'.format(self.t)+'.h5'

        self._file_exist(fno)

        h5file = h5py.File(fno, 'w')

        for field in fields:
            h5file.create_dataset(field, data=eval("self."+field))

        h5file.close()


    def _save_diagnostics(self, diagnostics=['t','ke','pe','ens']):

        """ Save diagnostics """

        fno = self.fno + 'diagnostics.h5'

        self._file_exist(fno)

        h5file = h5py.File(fno, 'w')

        for diagnostic in diagnostics:
            h5file.create_dataset(diagnostic, data=eval("self."+diagnostic))

        h5file.close()


    # step forward
    def _init_rk3w(self):
        """Initialize time-stepper stuff"""

        self.a1, self.a2, self.a3 = 29./96., -3./40., 1./6.
        self.b1, self.b2, self.b3 = 37./160., 5./24., 1./6.
        self.c1, self.c2, self.c3 = 8./15., 5./12., 3./4.
        self.d1, self.d2 = -17./60., -5./12.

        self.Lin = -self.nu*self.kappa2*self.dt
        self.L1 = ( (1. + self.a1*self.Lin)/(1. - self.b1*self.Lin) )     
        self.L2 = ( (1. + self.a2*self.Lin)/(1. - self.b2*self.Lin) )
        self.L3 = ( (1. + self.a2*self.Lin)/(1. - self.b3*self.Lin) )

    def _init_filter(self):
        """ Set spectral filter """

        if self.use_filter:
            cphi=0.65*pi
            cphi = 0.715*pi;
            wvx=sqrt((self.k*self.dx)**2.+(self.l*self.dz)**2.)
            self.filt = exp(-23.6*(wvx-cphi)**4.)
            self.filt[wvx<=cphi] = 1. 
        else:
            # if not use exponential filter,
            #   then dealias using 2/3 rule
            self.filt = np.ones_like(self.kappa2)
            self.filt[self.nx/3:2*self.nx/3,:] = 0.
            self.filt[:,self.nz/3:] = 0.

    # some diagnostics
    def _calc_cfl(self):
        return np.abs(
            np.hstack([self.u, self.v])).max()*self.dt/self.dx

    def _calc_ke(self):
        ke = .5*self.spec_var(self.kappa*self.ph)
        return ke.sum()

    def _calc_pe(self):
        pe = .5*self.spec_var(self.bh/self.Fr)
        return pe.sum()

    def _calc_ens(self):
        ens = .5*self.spec_var(self.kappa2*self.ph)
        return ens.sum()

    def _calc_anisotropy(self):
        return (self.u**2).mean()/(self.v**2).mean()

    def spec_var(self,ph):
        """ compute variance of p from Fourier coefficients ph """
        var_dens = 2. * np.abs(ph)**2 / (self.nx*self.nz)**2 
        # only half of coefs [0] and [nx/2+1] due to symmetry in real fft2
        var_dens[:,0],var_dens[:,-1] = var_dens[:,0]/2.,var_dens[:,-1]/2.
        return var_dens.sum()
 
