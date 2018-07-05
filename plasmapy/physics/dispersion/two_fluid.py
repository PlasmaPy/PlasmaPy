#
#
#
# FUNCTION TO CALCULATE PHASE SPEEDS OF THE THREE BRANCHES OF TWO FLUID
# DISPERSION RELATION (e.g. STRINGER JPP 1963, ROGERS PRL 2001, 
# and BELLAN JGR 2012)
#
#                   Tulasi Nandan Parashar
#
import numpy as np
def tfps(beta=0.6, ca=1., de2=0.000545, theta=0., kk=1.):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation in degrees
      kk: Wavenumber of interest in units of kdi
      
      Output is frequencies of the roots and the phase speeds w/k
      The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron
   """
   tht = theta*pi/180.
   ct = np.cos(tht)
   st = np.sin(tht)
   tt = st/ct
   cs=np.sqrt(beta/2.)*ca
   di =1.
   caksq=np.cos(tht)**2*ca**2
   cmsq=ca**2 + cs**2
   
   pcs=np.zeros(4)
   D = 1 + kk**2*de2
   # Find out the root of the quadratic dispersion relation
   pcs[0] = 1.
   pcs[1] = -(ca**2/D + cs**2 + caksq*(1+kk**2*di**2/D)/D)*kk**2
   pcs[2] = caksq*kk**4*(ca**2/D + cs**2 + cs**2*(1+kk**2*di**2/D))/D
   pcs[3] = -(cs**2*kk**6*caksq**2)/D**2
   w2 = np.roots(pcs); w = np.sqrt(w2)
   speeds= w/kk
   return w,speeds
#
# Compute the dispersion relation w(k) vs k for a single theta value.
#
def tfst(beta=0.6, ca=1., de2=0.000545, theta=0., kmin=1e-2, kmax=10., npoints=200,wrt='n'):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation in degrees
      kkmin: Minimum Wavenumber of interest in units of kdi
      kkmax: Maximum wavenumber of interest in units of kdi

      Output is an array with 4 columns, k, w-fast, w-alf, w-slow
   """
   import matplotlib.pyplot as plt
   kmmn=np.log10(kmin)
   kmmx=np.log10(kmax)
   kk=np.logspace(kmmn,kmmx,npoints)
   warray=np.zeros((3,npoints))
   for i in range(0,npoints):
      f,s = tfps(beta, ca, de2, theta, kk[i])
      warray[:,i]=f
   plt.loglog(kk,warray[0,:], label='Fast/Magnetosonic')
   plt.loglog(kk,warray[1,:], label='Alfven/KAW')
   plt.loglog(kk,warray[2,:], label='Slow')
   plt.xlabel('$kd_i$')
   plt.ylabel('$\omega/\omega_{ci}$')
   plt.legend(loc='best',fancybox=True,framealpha=0.2)
   plt.title('Dispersion Relation for beta='+str(beta)+' and me/mi='+str(de2))
   plt.show()
   if wrt == 'y':
      ofile=open('disp'+str(theta)+'.dat','w')
      print>> ofile,'#  k', 'Acoustic', 'Alfven', 'Fast'
      for i in range(npoints):
         print>> ofile, kk[i],warray[2,i],warray[1,i],warray[0,i]
      ofile.close()

def tfev(beta=0.6, ca=1., de2=0.000545, theta=0., k=1.,aa=0.1):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation in degrees
      k: wavenumber of interest in units of kdi

      Output: Prints the two fluid eigenvector to the screen
   """
   def amp(beta,de2,k,w,theta,aa):
      th=theta*pi/180.
      bb=1-w**2*(1+de2*k**2)/(k**2*np.cos(th)**2)
      sk='sin('+str(round(k,3))+'x)'
      ck='cos('+str(round(k,3))+'x)'
      def st(a):
         return str(round(a,3))
      return 'bx=0.'+\
            '\tby = '+st(round(2*aa,3))+ck+\
            '\tbz = '+st(-np.cos(th)*bb*2*aa/w)+sk+\
            '\nux = '+st(aa*k*bb*np.sin(2*th)/(w**2-beta*k**2))+sk+\
            '\tuy = '+st(-2*aa*k*np.cos(th)/w)+ck+\
            '\tuz = '+st(2*aa*k*bb*np.cos(th)**2/w**2)+sk+\
            '\t n = '+st((k*np.cos(th)/w)*aa*k*bb*np.sin(2*th)/(w**2-beta*k**2))+sk
   f,s=tfps(beta,ca,de2,theta,k)
#     The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron
   print 'Fast/ Whistler'
   print amp(beta,de2,k,f[0],theta,aa)
   print '##################'
   print
   print 'Alfven/ KAW'
   print amp(beta,de2,k,f[1],theta,aa)
   print '##################'
   print
   print 'Slow/ Cyclotron'
   print amp(beta,de2,k,f[2],theta,aa)
