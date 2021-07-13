from __future__ import division
from __future__ import print_function
from builtins import zip
from past.utils import old_div
import numpy as np
import os
import sys
#sys.path.insert(1,os.path.abspath('..'))
import fmmfronthelper
import param

# compute fronts at a given times for a set of features
# Francesca Tesser and Wolfram Moebius
# April - August 2016 - May 2021

# four input arguments
print('Note: four input parameters, index for rho, index for gamma, the number of the realization (me) and the seed for the random number generator.')
if not len(sys.argv)==5:
   sys.exit('Wrong number of input arguments.')
rhoid=int(sys.argv[1])
gammaid=int(sys.argv[2])
me=int(sys.argv[3])
randseed=int(sys.argv[4])
print('rhoid:',rhoid)
print('gammaid:',gammaid)
print('me:',me)
print('random number seed:',randseed)
if randseed<0 or randseed>2**32-1:
  sys.exit('Invalid random number seed.')
np.random.seed(randseed)

# setting the stage, defining the domain
auxhalfls=0.5/param.invdx
X,Y=np.meshgrid(np.linspace(-param.neglengthch+auxhalfls,param.lengthch-auxhalfls,np.rint((param.neglengthch+param.lengthch)*param.invdx).astype(int)),np.linspace(auxhalfls,param.widthch-auxhalfls,np.rint(param.widthch*param.invdx).astype(int)))
# setting the stage, initial condition
phi=np.ones_like(X)
phi[:,0]=0 # this choice matters for front finding

# loop over parameters number density
for rho in param.rhos[rhoid:rhoid+1]:
  # number of features
  nspots=int(round(param.widthch*param.lengthch*rho))
  nspotsroughen=int(round(param.widthch*param.neglengthch/2.0*rho))

  # loop over parameter hotspot / obstacle intensity
  for gamma,gamma_shakeup in zip(param.gammas[gammaid:gammaid+1],param.gammas_shakeup[gammaid:gammaid+1]):
    vfeat=param.vbasic*gamma; # feature speed
    vfeat_shakeup=param.vbasic*gamma_shakeup; # feature speed for shake-up features (placed in front of ROI)
    vback=param.vbasic; # background speed

    # create folder for output
    folder="OUTPUT/data.gamma.{0:.6f}/rho.{1:.6f}".format(gamma,rho)
    if not os.path.exists(folder):
      os.makedirs(folder)
    print('rho, gamma, me:', rho, gamma, me)

    # construct speedmatrix
    # set everything to background speed
    speedmatrix=vback*np.ones_like(X)
    # axes of the features with elliptical shape
    axis1=param.radiusfeat*param.ratiofeat
    axis2=param.radiusfeat

    # shake up the front a little bit before the actual ROI - use gamma_shakeup=1 if this is not desired
    # coordinates of feature positions
    featpos=np.random.rand(nspotsroughen,2)
    featpos[:,0]=-featpos[:,0]*param.neglengthch/2.0
    featpos[:,1]=featpos[:,1]*param.widthch
    # place individual features
    for onefeat in featpos:
      speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1])**2,axis2**2)<=1] = vfeat_shakeup 
      # because of periodic boundary conditions
      speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1]-param.widthch)**2,axis2**2)<=1] = vfeat_shakeup 
      speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1]+param.widthch)**2,axis2**2)<=1] = vfeat_shakeup 

    if param.refraction==0:
      # coordinates of feature positions
      featpos=np.random.rand(nspots,2)
      # can ignore that there should be features past the length of the channel which overlap with the channel b/c we won't fit in this range to get speed 
      featpos[:,0]=featpos[:,0]*param.lengthch
      featpos[:,1]=featpos[:,1]*param.widthch
      # place individual features
      for onefeat in featpos:
        speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1])**2,axis2**2)<=1] = vfeat
        # because of periodic boundary conditions
        speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1]-param.widthch)**2,axis2**2)<=1] = vfeat
        speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1]+param.widthch)**2,axis2**2)<=1] = vfeat
    if param.refraction==1:
      # coordinates of feature positions
      featpos=np.random.rand(nspots*2,2)
      # can ignore that there should be features past the length of the channel which overlap with the channel
      featpos[:,0]=featpos[:,0]*param.lengthch
      featpos[:,1]=featpos[:,1]*2.0*param.widthch-0.5*param.widthch
      # place individual features
      for onefeat in featpos:
        # if in left part of channel (tilted boundary)
        # if in right part of channel, but random number < 0.1
        if onefeat[1]-0.5*param.widthch>onefeat[0]-0.5*param.lengthch or np.random.random_sample()<0.1:
          speedmatrix[old_div((X-onefeat[0])**2,axis1**2)+old_div((Y-onefeat[1])**2,axis2**2)<=1] = vfeat

    # and the actual work computation
    # time interval is chosen so that front is detected every time it has passed unit length (in or ourside features, whatever is faster)
    # boundary condition determined by param.refraction 
    auxphi=1.0-np.exp(-np.pi*axis1*axis2*rho)
    auxdt=(1-auxphi)/vback+auxphi/np.amax([1,vfeat])
    fmmfronthelper.fmmandfront(X, Y, phi, speedmatrix, (1-param.refraction,0), auxdt, 0, folder, me, param.dump_data==1 or (param.dump_data==2 and me==1))
