import numpy as np
from scipy.special import factorial, lpmn
import quadpy

from mayavi import mlab

class sphbasis:
    '''
    Class for constructing spherical harmonics (Ylms), their gradients
    and related magnetic field 'basis vectorfunctions'
    (Psilms for r**l components, Philms for r**(-l) components).
    
    Uses notations and definitions by Plattner and Simons (2014; https://arxiv.org/pdf/1306.3201.pdf)
    and the same normalization conventions.
    
    Integration over a surface of unit sphere is 
    used as the inner product <C,D> = int C dot D dOmega.
    
    Has also lot of functions for spherical <-> cartesian transformations.
    
    TODO: mu0 might be missing!!!
    '''  
    
    def __init__(self,Np):
        '''
        Initialises the sphbasis object.
        
        Parameters:
            
            Np: single, number of points along theta and phi in spherical coordinates
            
        Returns:
            
            sphbasis object with fields:
                sp: Nx3 array - spherical coordinates of the points [r, theta, phi]
                p: Nx3 array - cartesian coordinates
                qp: Mx3 array - cartesian quadrature points
                sqp: Mx3 array - spherical quadrature points
            
        '''
        
        theta = np.linspace(0.01, np.pi-0.01, Np)
        phi = np.linspace(0, 2*np.pi, Np)
        phi, theta = np.meshgrid(phi, theta)
        phi = phi.flatten()
        theta = theta.flatten()
        
        self.sp = np.zeros((theta.shape[0],3))
        self.sp[:,0] = 1     
        self.sp[:,1] = theta
        self.sp[:,2]= phi
        
        self.p = self.spherical2cartesian(self.sp)
        self.Np = Np
        
        self.initqp()
            
    def initqp(self):
        '''
        Initialises quadrature points on the sphere.
        
        Default points are McLaren(10) so that we avoid singularities.
        '''
        
        self.qp = quadpy.sphere.McLaren(10)
        sp = self.cartesian2spherical(self.qp.points)
        self.sqp = sp
        
    def cartesian2spherical(self,p):
        '''
        Maps cartesian coordinates to spherical.
        
        Parameters:
            
            p: Nx3 array - cartesian coordinates
            
        Returns:
            
            sp: Nx3 array - spherical coordinates [r, theta, phi]
        
        '''
        
        r = np.sqrt(p[:,0]**2+p[:,1]**2+p[:,2]**2)
        theta = np.arctan2(np.sqrt(p[:,0]**2+p[:,1]**2),p[:,2])
        phi = np.arctan2(p[:,1],p[:,0])
        
        sp = np.array((r,theta,phi)).T
        
        return sp
        
    def spherical2cartesian(self,sp):
        '''
        Maps spherical coordinates to cartesian.
        
        Parameters:
            
            sp: Nx3 array - spherical coordinates [r, theta, phi]
            
        Returns:
            
            p: Nx3 array - cartesian croodinates
        '''
        
        X = sp[:,0] * np.sin(sp[:,1]) * np.cos(sp[:,2])
        Y = sp[:,0] * np.sin(sp[:,1]) * np.sin(sp[:,2])
        Z = sp[:,0] * np.cos(sp[:,1]) 
        
        p = np.array((X,Y,Z)).T
        return p
    
    
    def Rotmatrix(self,sp):
        '''
        Constructs rotation matrix from cartesian coordinates to spherical.
        
        Parameters:
            
            sp: Nx3 array - spherical coordinates [r, theta, phi]
            
        Returns:
            
            vmat: 3x3 array - rotation matrix from cartesian to spherical.
        '''
        
        vmat = np.zeros((3,3))     
        vmat[0,0] = np.sin(sp[1])*np.cos(sp[2])
        vmat[0,1] = np.sin(sp[1])*np.sin(sp[2])
        vmat[0,2] = np.cos(sp[1])
        vmat[1,0] = np.cos(sp[1])*np.cos(sp[2])   
        vmat[1,1] = np.cos(sp[1])*np.sin(sp[2])   
        vmat[1,2] = -1*np.sin(sp[1])
        vmat[2,0] = -1*np.sin(sp[2])
        vmat[2,1] = np.cos(sp[2])
        vmat[2,2] = 0
        
        return vmat

    def cartvec2sph(self,sp,vec):
        '''
        Transforms cartesian vector to spherical coordinates.
        
        Parameters:
            
            sp: Nx3 array - spherical coordinates [r, theta, phi]
            vec: Nx3 array - vector in cartesian coordinates [e_x, e_y, e_z]
            
        Returns:
            
            svec: Nx3 array - vector in spherical coordinates [e_r, e_theta, e_phi]
        
        '''
        
        svec = np.zeros(vec.shape)
        for i in range(sp.shape[0]):
            vmat = self.Rotmatrix(sp[i])
            svec[i] = vmat@vec[i]
        return svec
    
    def sphvec2cart(self,sp,vec):
        '''
        Transforms cartesian vector to spherical coordinates.
        
        Parameters:
            
            sp: Nx3 array - spherical coordinates [r, theta, phi]
            vec: Nx3 array - vector in spherical coordinates [e_r, e_theta, e_phi]
            
        Returns:
            
            svec: Nx3 array - vector in cartesian coordinates [e_x, e_y, e_z]
        
        '''        
        
        svec = np.zeros(vec.shape)
        for i in range(sp.shape[0]):
            vmat = self.Rotmatrix(sp[i])
            svec[i] = np.transpose(vmat)@vec[i]
            
        return svec
    
    def lpmn_em(self,l,m,x):
        '''
        Computes associated Legendre function (Plm) of the first kind of order m and degree l.
        
        Parameters:
            
            l: single - degree of Plm
            m: single - order of Plm
            x: Nx1 array - evaluation points
            
        Returns:
            
            lp: Nx1 array - Plm at x
            
        '''
        
        lp = np.zeros(x.shape)
        for i in range(x.shape[0]):
            a, b = lpmn(m,l,x[i])
            lp[i] = a[np.abs(m),l]
        return lp
    
    def derlpmn_em(self,l,m,x):
        '''
        Computes derivative of associated Legendre function (Plm) of the first kind of order m and degree l
        with respect to the argument x.
        
        Parameters:
            
            l: single - degree of Plm
            m: single - order of Plm
            x: Nx1 array - evaluation points
            
        Returns:
            
            derlp: Nx1 array - dPlm/dx at x
            
        '''
        
        derlp = np.zeros(x.shape)
        for i in range(x.shape[0]):
            a, b = lpmn(m,l,x[i])
            derlp[i] = b[np.abs(m),l]
        return derlp
    
    def xlm(self, l, m, theta):
        '''
        Xlm-function used in the definition of spherical harmonics (Ylm).
        Follows notation of Plattner and Simons (2014);
        see Eqs. 1--3 in https://arxiv.org/pdf/1306.3201.pdf.
        
        Parameters:
            
            l: single - degree of Xlm
            m: single - order of Xlm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            
        Returns:
            
            Xlm: Nx1 array - Xlm at theta
            
        '''
        
        xlm = ((2*l+1)/(4*np.pi))**.5*(factorial(l-m)/factorial(l+m))**.5
        xlm *= self.lpmn_em(l,m,np.cos(theta))
        return xlm
        
    def ylm(self, l, m, theta, phi):
        '''
        Real spherical harmonics as defined by Plattner and Simons (2014);
        see Eqs. 1--3 in https://arxiv.org/pdf/1306.3201.pdf.
        
        Parameters:
            
            l: single - degree of Ylm
            m: single - order of Ylm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            Ylm: Nx1 array - Ylm at (theta,phi)
        
        '''
        
        if m < 0:
            ylm = np.sqrt(2)*self.xlm(l,m, theta)*np.cos(m*phi)
        elif m == 0:
            ylm = self.xlm(l,0,theta)
        elif m > 0:
            ylm = np.sqrt(2)*self.xlm(l,m, theta)*np.sin(m*phi)
            
        return ylm
    
    
    def derxlm(self, l, m, theta):
        '''
        Derivative of Xlm with respect to theta.
        
        Parameters:
            
            l: single - degree of Xlm
            m: single - order of Xlm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            
        Returns:
            
            derxlm: Nx1 array - dXlm/dtheta at theta
        
        '''
        
        derxlm = ((2*l+1)/(4*np.pi))**.5*(factorial(l-m)/factorial(l+m))**.5
        derxlm *= self.derlpmn_em(l,m,np.cos(theta))
        derxlm *= -1*np.sin(theta) #this comes from dXlm(cos(theta))/dtheta = dXlm(cos(theta))/dcos(theta)*(-sin(theta))
        return derxlm
    
    def sinxlm(self, l, m, theta):
        '''
        Computes m/(sin(theta))*Xlm.
        
        Parameters:
            
            l: single - degree of Xlm
            m: single - order of Xlm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            
        Returns:
            
            sinxlm: Nx1 array - m/(sin(theta))*Xlm at theta
        
        '''        
        
        sinxlm = m/(np.sin(theta))*self.xlm(l,m, theta)
        return sinxlm    
    
    def dthylm(self, l, m, theta, phi):
        '''
        Derivative of Ylm with respect to theta dYlm/dtheta.
        
        Parameters:
            
            l: single - degree of Ylm
            m: single - order of Ylm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            dthylm: Nx1 array - dYlm/dtheta at (theta, phi).
        
        '''    
        
        if m < 0:
            dthylm = np.sqrt(2)*self.derxlm(l,m, theta)*np.cos(m*phi)
        elif m == 0:
            dthylm = self.derxlm(l,0,theta)
        elif m > 0:
            dthylm = np.sqrt(2)*self.derxlm(l,m, theta)*np.sin(m*phi)
        return dthylm
    
    
    def dphiylm(self, l, m, theta, phi):
        '''
        Derivative of Ylm with respect to phi dYlm/dphi.
        
        Parameters:
            
            l: single - degree of Ylm
            m: single - order of Ylm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            dphiylm: Nx1 array - dYlm/dphi at (theta, phi).
        
        '''           
        
        if m < 0:
            dphiylm = -np.sqrt(2)*np.sin(m*phi)*self.sinxlm(l,m,theta)
        elif m == 0:
            dphiylm = 0
        elif m > 0:
            dphiylm = np.sqrt(2)*np.cos(m*phi)*self.sinxlm(l,m,theta)
        return dphiylm
    
    
    def Plm(self, l, m, theta, phi):
        '''
        Plm vector function (see Eq. 18 Plattner and Simons (2014)).
        
        Parameters:
            
            l: single - degree of Plm
            m: single - order of Plm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            Plm: Nx3 array - Plm at (theta, phi).
        
        '''         
        
        Plm = np.zeros((theta.shape[0],3))
        Plm[:,0] = self.ylm(l,m,theta,phi)
        return Plm

    def Blm(self, l, m, theta, phi):
        '''
        Blm vector function (see Eq. 19 Plattner and Simons (2014)).
        
        Parameters:
            
            l: single - degree of Plm
            m: single - order of Plm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            Blm: Nx3 array - Blm at (theta, phi).
        
        '''
        
        Blm = np.zeros((theta.shape[0],3))
        
        Blm[:,1] = self.dthylm(l,m,theta,phi)
        
        Blm[:,2] = self.dphiylm(l,m,theta,phi)
        
        Blm *= 1/np.sqrt(l*(l+1))
        return Blm
    
    def Psilm(self,l,m,theta,phi):
        '''
        Vector basis function (Psilm) for r**l component of the magnetic field.
        Normalization <Psilm,Psikn> = delta_lk,mn.
        
        Parameters:
            
            l: single - degree of Psilm
            m: single - order of Psilm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            Psilm: Nx3 array - Psilm at (theta, phi).
        
        '''        
        
        Psilm = l*self.Plm(l,m,theta,phi) + np.sqrt(l*(l+1))*self.Blm(l,m,theta,phi)
        Psilm *= 1/np.sqrt(2*l**2+l)
        return Psilm

    def Philm(self,l,m,theta,phi):
        '''
        Vector basis function (Philm) for r**(-l) component of the magnetic field.
        Normalization <Philm,Phikn> = delta_lk,mn.
        
        Parameters:
            
            l: single - degree of Philm
            m: single - order of Philm
            theta: Nx1 array - evaluation points, theta at spherical coordinates
            phi: Nx1 array - evaluation points, phi at spherical coordinates
            
        Returns:
            
            Philm: Nx3 array - Philm at (theta, phi).
        
        '''        
                
        Philm = -1*(l+1)*self.Plm(l,m,theta,phi) + np.sqrt(l*(l+1))*self.Blm(l,m,theta,phi)
        Philm *= 1/np.sqrt((l+1)*(2*l+1))
        return Philm
    
    def innerproduct(self,fun1,fun2):
        '''
        Inner product of vector functions fun1 and fun2.
        Defined as integration over a surface of unit sphere <C,D> = int C dot D dOmega.
        Quadrature rule defined in qp is used.
        
        Parameters:
            
            fun1: Nx3 array - vector function
            fun2: Nx3 array - vector function
            
        Returns:
            
            dotp: single - inner product of fun1 and fun2
        
        '''
        
        dotp = np.sum(self.qp.weights*np.sum(fun1*fun2,axis=1))*4*np.pi
        return dotp
    
    def avsphspectra(self, fun, lmax):
        '''
        Calculate the l,m-spectra (over r**l-terms) of vector function defined in quadrature points
        using the inner product.
        
        Parameters:
            
            fun: Nx3 array - vector function computed at quadrature points self.sqp
            lmax: single - maximum degree l for which the spectra is computed
            
        Returns:
            
            coeffs: lmax*(lmax+2)x1 array - spectral coefficients
        
        '''        
        
        coeffs = []

        for l in range(1,lmax+1):
            for m in range(-1*l,l+1):
                Psilm = self.Psilm(l,m,self.sqp[:,1],self.sqp[:,2])
                ctemp = self.innerproduct(fun, Psilm)
                ctemp /= (self.sqp[0,0]**(l-1)*np.sqrt(2*l**2 +l)) # we use this normalization
#                ctemp /= (self.sqp[0,0]**(l-1))
                coeffs.append(ctemp)
        
        coeffs = np.array(coeffs)
        return coeffs

    def bvsphspectra(self, fun, lmax):
        '''
        Calculate the l,m-spectra (over r**(-l)-terms) of vector function defined in quadrature points
        using the inner product.
        
        Parameters:
            
            fun: Nx3 array - vector function computed at quadrature points self.sqp
            lmax: single - maximum degree l for which the spectra is computed
            
        Returns:
            
            coeffs: lmax*(lmax+2)x1 array - spectral coefficients
        
        '''
        
        coeffs = []

        for l in range(1,lmax+1):
            for m in range(-1*l,l+1):
                Philm = self.Philm(l,m,self.sqp[:,1],self.sqp[:,2])
                ctemp = self.innerproduct(fun, Philm)
                ctemp /= (self.sqp[0,0]**(l-1)*np.sqrt((l+1)*(2*l+1))) # we use this normalization
#                ctemp /= (self.sqp[0,0]**(l-1))
                coeffs.append(ctemp)
        
        coeffs = np.array(coeffs)
        return coeffs
    
    def potential(self,p, acoeffs, bcoeffs, lmax):
        '''
        Computes magnetic scalar potential from the sph coefficients.
        Ignores the 'DC' component l=0.

        Parameters:
            
            p: Nx3 array - coordinates in which the potential is computed
            acoeffs: lmax*(lmax+2)x1 array - spectral coefficients of r**l terms
            bcoeffs: lmax*(lmax+2)x1 array - spectral coefficients of r**(-l) terms
            lmax: single - maximum degree l which is used in computing
            
        Returns:
            
            pot: Nx1 array - magnetic scalar potential at p
        
        '''
        
        pot = np.zeros(p.shape[0])
        
        sp = self.cartesian2spherical(p)
        
        lind = 0
        for l in range(1,lmax+1):
            for m in range(-1*l,l+1):
                ylm = self.ylm(l, m, sp[:,1], sp[:,2])
                pot += (acoeffs[lind]*sp[:,0]**l + bcoeffs[lind]*sp[:,0]**(-1*l-1))*ylm
                lind += 1
        return pot
    
        
class sphfittools:
    '''
    Class for fitting spherical harmonics basis functions to measured magnetic field data.
    Can be used to calculate the l,m-spectra of measured data and reconstructing field with
    the spherical harmonics representation.
    '''       
    
    def fitSpectra(sph, coords, Bmeas, lmax):
        '''
        Fits spherical harmonics representation (r**l) to measured data.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            coords: Nx3x3 array - measurement coordinates, each measured field direction
            in the third dimension: e.g. coords[:,:,2] gives the coordinates of measured z-components.
            Bmeas: Nx3 array - the measured field values along different directions (x,y,z)
            lmax: single - maximum degree l for which the fit is done
            
        Returns:
            
            coeffs: lmax*(lmax+2)x1 array - the unnormalized coefficients
            coeffs2: lmax*(lmax+2)x1 array - the 'properly' normalized squared coefficients
            nrmse: single - normalized rms error in percents between the data and fit
        
        '''        
        
        Nmeas = coords.shape[0]
        A = np.zeros((3*Nmeas, lmax*(lmax+2))) #initialize the fitting matrix
        
        #loop over the components
        for e in range(3):
            p = coords[:,:,e]
            sp = sph.cartesian2spherical(p)
            
            lind = 0
            for l in range(1,lmax+1):
                for m in range(-1*l,l+1):
                    Psilm = sph.Psilm(l,m, sp[:,1],sp[:,2])
                    Psilm *= np.sqrt(2*l**2 + l)
                    Psilm[:,0] *= sp[:,0]**(l-1)
                    Psilm[:,1] *= sp[:,0]**(l-1)
                    Psilm[:,2] *= sp[:,0]**(l-1)
                    Psilm = sph.sphvec2cart(sp, Psilm)
                    A[e*Nmeas:(e+1)*Nmeas, lind] = Psilm[:,e]
                    lind += 1
        print("Condition number = %f" % (np.linalg.cond(A))) #print the condition number of A
        
        coeffs = np.linalg.pinv(A)@Bmeas.T.flatten() #compute coefficients using pseudoinverse of A
        
        #the following script calculates the normalized coefficients
        #the coefficients are normalized so that squared norm of magnetic field integrates to 1 over the measurement volume
        coeffs2 = np.zeros(coeffs.shape)
        lind = 0
        Rmax = np.max(sp[:,0])
        for l in range(1,lmax+1):
                for m in range(-1*l,l+1):
#                    coeffs2[lind] = coeffs[lind]*np.sqrt(2*l**2 + l)
                    temp = (2*l**2 + l)*Rmax**(2*l-1)/(2*l-1)
                    coeffs2[lind] = coeffs[lind]**2*temp
                    lind += 1
        Breco = A@coeffs
        
        nrmse = np.sqrt(np.mean((Bmeas.T.flatten()-Breco)**2))/np.max(np.abs(Bmeas.T.flatten()))*100
        print("Normalized RMS error = %f%%" % (nrmse)) #print the normalized rms error
        
        return coeffs, coeffs2, nrmse
    
    def reconstructB(sph, p, coeffs,lmax):
        '''
        Reconstructs the magnetic field using the spherical harmonics coefficients.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            p: Nx3 array - coordinates where B is reconstructed
            coeffs: lmax*(lmax+2)x1 array - the unnormalized l,m-coefficients
            lmax: single - maximum degree l of the fit
            
        Returns:
            
            B: Nx3 array - reconstructed magnetic field at p
        
        '''          
        
        B = np.zeros(p.shape)
        sp = sph.cartesian2spherical(p)
        idx = 0
        for l in range(1,lmax+1):
                for m in range(-1*l,l+1):
                    Psilm = sph.Psilm(l,m, sp[:,1],sp[:,2])
                    Psilm *= np.sqrt(2*l**2 + l)
                    Psilm[:,0] *= sp[:,0]**(l-1)
                    Psilm[:,1] *= sp[:,0]**(l-1)
                    Psilm[:,2] *= sp[:,0]**(l-1)
                    Psilm *= coeffs[idx]
                    Psilm = sph.sphvec2cart(sp, Psilm)
                    B += Psilm
                    idx += 1
        return B
    
class plotsph:
    '''
    Class for visualization of spherical harmonics and basis vector functions.
    '''
    
    def plotYlms(sph, lmax):
        '''
        Plots real spherical harmonics up to lmax.
        Inspired by https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            lmax: single - maximum degree l
        
        '''
        
        theta = np.reshape(sph.sp[:,1], (sph.Np, sph.Np))
        phi = np.reshape(sph.sp[:,2], (sph.Np, sph.Np))
        r = 0.3
        x=r*np.sin(theta)*np.cos(phi)
        y=r*np.sin(theta)*np.sin(phi)
        z=r*np.cos(theta)
        
        
        for l in range(1, lmax+1):
            for m in range(l):
                ylm = sph.ylm(l,m,theta.flatten(),phi.flatten())
                ylm = np.reshape(ylm, (sph.Np, sph.Np))

                mlab.mesh(x - m, y - l, z, scalars=ylm, colormap='bwr')
                ylm /= ylm.max()
                mlab.mesh(ylm * x - m, ylm * y - l, ylm * z + 1.3,
                          scalars=np.abs(ylm), colormap='Spectral')
        
        mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))
        

    def plotYlm(sph, l,m):
        '''
        Plots real spherical harmonics of order m and degree l.
        Inspired by https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html
        
        Parameters:
            
            sph: spherical harmonics analysis object
            l: single - degree l
            m: single - order m
        
        '''
        
        theta = np.reshape(sph.sp[:,1], (sph.Np, sph.Np))
        phi = np.reshape(sph.sp[:,2], (sph.Np, sph.Np))
        r = 0.6
        x=r*np.sin(theta)*np.cos(phi)
        y=r*np.sin(theta)*np.sin(phi)
        z=r*np.cos(theta)
        
        ylm = sph.ylm(l,m,theta.flatten(),phi.flatten())
        ylm = np.reshape(ylm, (sph.Np, sph.Np))

        mlab.mesh(x - m, y - l, z, scalars=ylm, colormap='bwr')

        ylm /= ylm.max()
        mlab.mesh(ylm * x - m, ylm * y - l, ylm * z + 1.3,
                  scalars=np.abs(ylm), colormap='Spectral')
    
            
    def plotPsilm(sph,l, m):
        '''
        Plots magnetic field basis function 'Psilm' (r**l) over a sphere.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            l: single - degree l
            m: single - order m
            
        Returns:
            
            obj: mayavi object
            
        '''
        
        Psilm = sph.Psilm(l,m, sph.sp[:,1], sph.sp[:,2])
        Psilm = sph.sphvec2cart(sph.sp, Psilm)
        obj = mlab.quiver3d(sph.p[:,0],sph.p[:,1],sph.p[:,2], Psilm[:,0], Psilm[:,1], Psilm[:,2])
        obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))        
        return obj
    
    def plotBPsilm_volume(sph,l, m, lim, Np, offset):
        '''
        Plots magnetic field basis function 'Psilm' (r**l) over a 3D volume.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            l: single - degree l
            m: single - order m
            lim: single - limits for coordinates, e.g., xmin = -lim, xmax = lim
            Np: single - number of points along different coordinates
            offset: 1x3 array - offset of the volume in which Psilm is plotted
            
        Returns:
            
            obj: mayavi object   
            
        '''
        
        x, y, z = np.meshgrid(np.linspace(-lim+offset[0],lim+offset[0],Np),np.linspace(-lim+offset[1],lim+offset[1],Np),np.linspace(-lim+offset[2],lim+offset[2],Np))
        
        p = np.array((x.flatten(), y.flatten(), z.flatten())).T
        sp = sph.cartesian2spherical(p)
        
        Psilm = sph.Psilm(l,m, sp[:,1], sp[:,2])
        Psilm *= np.sqrt(2*l**2 +l)
        Psilm[:,0] *= sp[:,0]**(l-1)
        Psilm[:,1] *= sp[:,0]**(l-1)
        Psilm[:,2] *= sp[:,0]**(l-1)
        
        Psilm = sph.sphvec2cart(sp, Psilm)
        obj = mlab.quiver3d(p[:,0],p[:,1],p[:,2], Psilm[:,0], Psilm[:,1], Psilm[:,2])
        obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))        
        return obj
    
    

    def plotPhilm(sph,l, m):
        '''
        Plots magnetic field basis function 'Philm' (r**(-l)) over a sphere.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            l: single - degree l
            m: single - order m

        Returns:
            
            obj: mayavi object
            
        '''
        
        Philm = sph.Philm(l,m, sph.sp[:,1], sph.sp[:,2])
        Philm = sph.sphvec2cart(sph.sp, Philm)
        obj = mlab.quiver3d(sph.p[:,0],sph.p[:,1],sph.p[:,2], Philm[:,0], Philm[:,1], Philm[:,2])
        obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))        
        return obj
    
    def plotBPhilm_volume(sph,l, m, lim, Np, offset):
        '''
        Plots magnetic field basis function 'Philm' (r**(-l)) over a 3D volume.
        
        Parameters:
            
            sph: spherical harmonics analysis object
            l: single - degree l
            m: single - order m
            lim: single - limits for coordinates, e.g., xmin = -lim, xmax = lim
            Np: single - number of points along different coordinates
            offset: 1x3 array - offset of the volume in which Philm is plotted
            
        Returns:
            
            obj: mayavi object   
            
        '''
        
        x, y, z = np.meshgrid(np.linspace(-lim+offset[0],lim+offset[0],Np),np.linspace(-lim+offset[1],lim+offset[1],Np),np.linspace(-lim+offset[2],lim+offset[2],Np))
        
        p = np.array((x.flatten(), y.flatten(), z.flatten())).T
        sp = sph.cartesian2spherical(p)
        
        Philm = sph.Philm(l,m, sp[:,1], sp[:,2])
        Philm *= np.sqrt((l+1)*(2*l+1))
        
        Philm[:,0] *= sp[:,0]**(-1*(l+2))
        Philm[:,1] *= sp[:,0]**(-1*(l+2))
        Philm[:,2] *= sp[:,0]**(-1*(l+2))
        
        Philm = sph.sphvec2cart(sp, Philm)
        obj = mlab.quiver3d(p[:,0],p[:,1],p[:,2], Philm[:,0], Philm[:,1], Philm[:,2])
        obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))        
        return obj
        
    

if __name__ == '__main__':
    
    sph = sphbasis(20) 
    plotsph = plotsph
    
    obj = plotsph.plotYlms(sph,4)
    obj = plotsph.plotYlm(sph,3,3)
    
#    obj = plotsph.plotPsilm(sph,5,2)
#    obj = plotsph.plotPhilm(sph,2,0)
     
    offset = np.array((0, 0, 2))
    mlab.figure()
    obj = plotsph.plotBPhilm_volume(sph,5,0, 1, 10,offset)
    
    mlab.figure()
    obj = plotsph.plotBPsilm_volume(sph,5,0, 1, 10,offset)
    
    
    
    
    Psilm1 = sph.Psilm(1,0, sph.sqp[:,1], sph.sqp[:,2])
    Psilm2 = sph.Psilm(7,0, sph.sqp[:,1], sph.sqp[:,2])
    
#    print(sph.innerproduct(Psilm1,Psilm2))
#    
#    Philm1 = sph.Philm(1,0, sph.sqp[:,1], sph.sqp[:,2])
#    Philm2 = sph.Philm(7,0, sph.sqp[:,1], sph.sqp[:,2])
#    
#    print(sph.innerproduct(Philm2,Philm2))
    
#    
#    
#    B = np.zeros(sph.sqp.shape)
#    #B[:,0] = 1
#    B[:,2] = sph.qp.points[:,0]/np.max(sph.qp.points[:,0])
#    B += 0.1*sci.random.randn(B.shape[0], B.shape[1])
#    
#    B = sph.cartvec2sph(sph.sqp,B)
#    
#    coeffs = sph.vsphspectra(B, 7)
#    
#    plt.figure()
#    plt.semilogy(coeffs**2)
#    
#    obj = plotsph.plotYlm(sph,5,3)
    
#    Np = 10
#    lim = 3
#    x, y, z = np.meshgrid(np.linspace(-lim,lim,Np),np.linspace(-lim,lim,Np),np.linspace(-lim,lim,Np))
#        
#    
##    p = np.array((x.flatten(), y.flatten(), z.flatten())).T
#    
#    p = np.array((x.flatten(), y.flatten(), np.zeros(y.flatten().shape))).T
#    lmax = 2
#    acoeffs = np.zeros(lmax*(lmax+2))
#    bcoeffs = np.zeros(lmax*(lmax+2))
#    acoeffs[7] = 1
##    bcoeffs[2] = 1
#    
#    pot = sph.potential(p, acoeffs,bcoeffs, lmax)
#    
#    pot = np.reshape(pot, x.shape)
#    
#    mlab.mesh(x[:,:,0], y[:,:,0] ,z[:,:,0] , scalars=pot[:,:,0], colormap='Spectral')
#    
#    coords = np.zeros((p.shape[0],p.shape[1],3))
#    coords[:,:,0] = p
#    coords[:,:,1] = p
#    coords[:,:,2] = p
#    
#    B = np.zeros((coords.shape[0],coords.shape[1]))
#    B[:,2] = p[:,0]/np.max(p[:,0])
#    B[:,1] = 0.3
#    B += 0.4*sci.random.randn(B.shape[0], B.shape[1])
#    
#    lmax = 5
#    coeffs,coeffs2, mse = sphfittools.fitSpectra(sph,coords, B, lmax)
#    
#    plt.figure()
#    plt.semilogy(coeffs**2,'.')
#    
    