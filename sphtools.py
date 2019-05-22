import numpy as np
#import numdifftools as nd
from scipy.special import lpmv, factorial, lpmn

import matplotlib.pyplot as plt

import quadpy
import scipy as sci

from mayavi import mlab

class sphbasis:
    
    """
    
    """    
    def __init__(self,Np):
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
        self.qp = quadpy.sphere.McLaren(10)
        sp = self.cartesian2spherical(self.qp.points)
        self.sqp = sp
        
    def cartesian2spherical(self,p):
        
        r = np.sqrt(p[:,0]**2+p[:,1]**2+p[:,2]**2)
        theta = np.arctan2(np.sqrt(p[:,0]**2+p[:,1]**2),p[:,2])
        phi = np.arctan2(p[:,1],p[:,0])
        
        sp = np.array((r,theta,phi)).T
        
        return sp
        
    def spherical2cartesian(self,sp):
        
        X = sp[:,0] * np.sin(sp[:,1]) * np.cos(sp[:,2])
        Y = sp[:,0] * np.sin(sp[:,1]) * np.sin(sp[:,2])
        Z = sp[:,0] * np.cos(sp[:,1]) 
        
        p = np.array((X,Y,Z)).T
        return p
    
    
    def Rotmatrix(self,sp):
        vmat = np.zeros((3,3))        
        temp = np.sin(sp[1])*np.cos(sp[2])
        vmat[0,0] = temp
        
        temp = np.sin(sp[1])*np.sin(sp[2])
        vmat[0,1] = temp
        
        temp = np.cos(sp[1])
        vmat[0,2] = temp
        
        temp = np.cos(sp[1])*np.cos(sp[2])
        vmat[1,0] = temp   
        
        temp = np.cos(sp[1])*np.sin(sp[2])
        vmat[1,1] = temp   
        
        temp = -1*np.sin(sp[1])
        vmat[1,2] = temp
        
        temp = -1*np.sin(sp[2])
        vmat[2,0] = temp
        
        temp = np.cos(sp[2])
        vmat[2,1] = temp
    
        vmat[2,2] = 0
        
        return vmat

    def cartvec2sph(self,sp,vec):
        svec = np.zeros(vec.shape)
        for i in range(sp.shape[0]):
            vmat = self.Rotmatrix(sp[i])
            svec[i] = vmat@vec[i]
        return svec
    
    def sphvec2cart(self,sp,vec):
        svec = np.zeros(vec.shape)
        for i in range(sp.shape[0]):
            vmat = self.Rotmatrix(sp[i])
            svec[i] = np.transpose(vmat)@vec[i]
            
        return svec
    
    def lpmv_em(self,l,m,x):
        lp = lpmv(m,l,x)
        return lp
    
    def lpmn_em(self,l,m,x):
        
        lp = np.zeros(x.shape)
        for i in range(x.shape[0]):
            a, b = lpmn(m,l,x[i])
            lp[i] = a[np.abs(m),l]
        return lp
    
    def derlpmn_em(self,l,m,x):
        
        derlp = np.zeros(x.shape)
        for i in range(x.shape[0]):
            a, b = lpmn(m,l,x[i])
            derlp[i] = b[np.abs(m),l]
        return derlp
    
    def xlm(self, l, m, theta):
        xlm = ((2*l+1)/(4*np.pi))**.5*(factorial(l-m)/factorial(l+m))**.5
        xlm *= self.lpmn_em(l,m,np.cos(theta))
        return xlm
        
    def ylm(self, l, m, theta, phi):
        if m < 0:
            ylm = np.sqrt(2)*self.xlm(l,m, theta)*np.cos(m*phi)
        elif m == 0:
            ylm = self.xlm(l,0,theta)
        elif m > 0:
            ylm = np.sqrt(2)*self.xlm(l,m, theta)*np.sin(m*phi)
            
        return ylm
    
    
    def derxlm(self, l, m, theta):
        derxlm = ((2*l+1)/(4*np.pi))**.5*(factorial(l-m)/factorial(l+m))**.5
        derxlm *= self.derlpmn_em(l,m,np.cos(theta))
        derxlm *= -1*np.sin(theta)
        return derxlm
    
    def sinxlm(self, l, m, theta):
        sinxlm = m/(np.sin(theta))*self.xlm(l,m, theta)
        return sinxlm    
    
    def dthylm(self, l, m, theta, phi):
        if m < 0:
            dthylm = np.sqrt(2)*self.derxlm(l,m, theta)*np.cos(m*phi)
        elif m == 0:
            dthylm = self.derxlm(l,0,theta)
        elif m > 0:
            dthylm = np.sqrt(2)*self.derxlm(l,m, theta)*np.sin(m*phi)
        return dthylm
    
    
    def dphiylm(self, l, m, theta, phi):
        if m < 0:
            dphiylm = -np.sqrt(2)*np.sin(m*phi)*self.sinxlm(l,m,theta)
        elif m == 0:
            dphiylm = 0
        elif m > 0:
            dphiylm = np.sqrt(2)*np.cos(m*phi)*self.sinxlm(l,m,theta)
        return dphiylm
    
    
    def Plm(self, l, m, theta, phi):
        Plm = np.zeros((theta.shape[0],3))
        Plm[:,0] = self.ylm(l,m,theta,phi)
        return Plm

    def Blm(self, l, m, theta, phi):
        Blm = np.zeros((theta.shape[0],3))
        
        Blm[:,1] = self.dthylm(l,m,theta,phi)
        
        Blm[:,2] = self.dphiylm(l,m,theta,phi)
        
        Blm *= 1/np.sqrt(l*(l+1))
        return Blm
    
    def Psilm(self,l,m,theta,phi):
        Psilm = l*self.Plm(l,m,theta,phi) + np.sqrt(l*(l+1))*self.Blm(l,m,theta,phi)
        Psilm *= 1/np.sqrt(2*l**2+l)
        return Psilm

    
    def innerproduct(self,fun1,fun2):
        dotp = np.sum(self.qp.weights*np.sum(fun1*fun2,axis=1))*4*np.pi
        return dotp
    
    def vsphspectra(self, fun, lmax):
        coeffs = []

        for l in range(1,lmax):
            for m in range(-1*l,l+1):
                Psilm = self.Psilm(l,m,self.sqp[:,1],self.sqp[:,2])
                ctemp = self.innerproduct(fun, Psilm)
                ctemp /= (self.sqp[0,0]**(l-1)*np.sqrt(2*l**2 +l))
#                ctemp /= (self.sqp[0,0]**(l-1))
                coeffs.append(ctemp)
        
        coeffs = np.array(coeffs)
        return coeffs

        
    
class sphfittools:
        
    def fitSpectra(sph, coords, Bmeas, lmax):
        Nmeas = coords.shape[0]
        A = np.zeros((3*Nmeas, lmax*(lmax+2)))
        

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
        print("Condition number = %f" % (np.linalg.cond(A)))
        coeffs = np.linalg.pinv(A)@Bmeas.T.flatten()
        
        coeffs2 = np.zeros(coeffs.shape)
        lind = 0
        for l in range(1,lmax+1):
                for m in range(-1*l,l+1):
                    coeffs2[lind] = coeffs[lind]/np.sqrt(2*l**2 + l)
                    lind += 1
        Breco = A@coeffs
#        mse = np.mean((Bmeas.T.flatten()-Breco)**2)
        nrmse = np.sqrt(np.mean((Bmeas.T.flatten()-Breco)**2))/np.max(np.abs(Bmeas.T.flatten()))*100
        print("Normalized RMS error = %f%%" % (nrmse))
        return coeffs, coeffs2, nrmse
    
    def reconstructB(sph, p, coeffs,lmax):
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
    
    def plotYlms(sph, lmax):
        #https://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html
        theta = np.reshape(sph.sp[:,1], (sph.Np, sph.Np))
        phi = np.reshape(sph.sp[:,2], (sph.Np, sph.Np))
        r = 0.3
        x=r*np.sin(theta)*np.cos(phi)
        y=r*np.sin(theta)*np.sin(phi)
        z=r*np.cos(theta)
        
        
        for l in range(1, lmax):
            for m in range(l):
                ylm = sph.ylm(l,m,theta.flatten(),phi.flatten())
                ylm = np.reshape(ylm, (sph.Np, sph.Np))

                mlab.mesh(x - m, y - l, z, scalars=ylm, colormap='bwr')
                ylm /= ylm.max()
                mlab.mesh(ylm * x - m, ylm * y - l, ylm * z + 1.3,
                          scalars=np.abs(ylm), colormap='Spectral')
        
        mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))
        

    def plotYlm(sph, l,m):
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
        Psilm = sph.Psilm(l,m, sph.sp[:,1], sph.sp[:,2])
        Psilm = sph.sphvec2cart(sph.sp, Psilm)
        obj = mlab.quiver3d(sph.p[:,0],sph.p[:,1],sph.p[:,2], Psilm[:,0], Psilm[:,1], Psilm[:,2])
        obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))        
        return obj
    
    def plotPsilm_volume(sph,l, m, lim, Np):
        x, y, z = np.meshgrid(np.linspace(-lim,lim,Np),np.linspace(-lim,lim,Np),np.linspace(-lim,lim,Np))
        
        p = np.array((x.flatten(), y.flatten(), z.flatten())).T
        sp = sph.cartesian2spherical(p)
        
        Psilm = sph.Psilm(l,m, sp[:,1], sp[:,2])
        
        Psilm[:,0] *= sp[:,0]**(l-1)
        Psilm[:,1] *= sp[:,0]**(l-1)
        Psilm[:,2] *= sp[:,0]**(l-1)
        
        Psilm = sph.sphvec2cart(sp, Psilm)
        obj = mlab.quiver3d(p[:,0],p[:,1],p[:,2], Psilm[:,0], Psilm[:,1], Psilm[:,2])
        obj.glyph.glyph_source.glyph_source.center = np.array((0, 0, 0))        
        return obj
    
        
    
    
#https://arxiv.org/pdf/1306.3201.pdf    


if __name__ == '__main__':
    
    sph = sphbasis(40) 
    plotsph = plotsph
    
#    obj = plotsph.plotYlms(sph,8)
#    obj = plotsph.plotYlm(sph,3,3)
    
#    obj = plotsph.plotPsilm(sph,2,0)
#    obj = plotsph.plotPsilm_volume(sph,5,3, 4.34, 7)
#    
#    Psilm1 = sph.Psilm(1,0, sph.sqp[:,1], sph.sqp[:,2])
#    Psilm2 = sph.Psilm(7,0, sph.sqp[:,1], sph.sqp[:,2])
#    
#    print(sph.innerproduct(Psilm1,Psilm2))
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
    
    Np = 5
    lim = 1
    x, y, z = np.meshgrid(np.linspace(-lim+0.1,lim,Np),np.linspace(-lim+0.1,lim,Np),np.linspace(-lim+0.1,lim,Np))
        
    p = np.array((x.flatten(), y.flatten(), z.flatten())).T
    coords = np.zeros((p.shape[0],p.shape[1],3))
    coords[:,:,0] = p
    coords[:,:,1] = p
    coords[:,:,2] = p
    
    B = np.zeros((coords.shape[0],coords.shape[1]))
    B[:,2] = p[:,0]/np.max(p[:,0])
    B[:,1] = 0.3
    B += 0.4*sci.random.randn(B.shape[0], B.shape[1])
    
    lmax = 5
    coeffs,coeffs2, mse = sphfittools.fitSpectra(sph,coords, B, lmax)
    
    plt.figure()
    plt.semilogy(coeffs**2,'.')
    
    