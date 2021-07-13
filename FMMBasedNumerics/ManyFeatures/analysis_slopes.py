from __future__ import print_function
from builtins import range
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import os
import sys
import param # parameters for fmmmanyfeatures.py
import paramanalysis # parameters for analysis

# set seed for random numbers (bootstrapping)
np.random.seed(42)

# auxiliary function to check for monotonicity
def monotonic(x):
    increment=np.diff(x)
    return np.all(increment >= 0)

# create output folder(s) if necessary
if not os.path.exists('OUTPUT/'+paramanalysis.anfolder):
  os.makedirs('OUTPUT/'+paramanalysis.anfolder)
pngfolder='OUTPUT/'+paramanalysis.anfolder+'/png'
if not os.path.exists(pngfolder):
  os.makedirs(pngfolder)

# results for front speed and corresponding errors
# from fits to individual trajectories...
slopes_indfit=np.full((len(param.gammas),len(param.rhos)),np.nan)
errors_indfit=slopes_indfit.copy()
# ...and from fits to ensembles of trajectories
slopes_bootstrap=slopes_indfit.copy()
errors_bootstrap=slopes_indfit.copy()

# loop over hotspot/obstacle intensity
for idxgamma,gamma in np.ndenumerate(param.gammas):
  # loop over number density of obstacles
  for idxrho,rho in np.ndenumerate(param.rhos):

    print('rho, gamma:', rho, gamma)

    # read and organize data
    inputfolder='OUTPUT/data.gamma.{0:.6f}/rho.{1:.6f}'.format(gamma,rho)
    inputfiles=glob.glob(inputfolder+'/meanfront.me*')
    numinstances=len(inputfiles)
    tlist=[] # list containing times for each realization
    lengthtlist=[] # list of number of time points for each realization
    poslist=[] # list containing mean front pos for each realization
    widthlist=[] # list containing front width for each realization
    # read all files and concatenate data
    for myfile in inputfiles:
      aux=np.loadtxt(myfile)
      if aux.size==0:
        sys.exit('Exit: encountered file without data.')
      if not monotonic(aux[:,0]):
        sys.exit('Exit: non monotonic time array!!') 
      tlist.append(aux[:,0])
      lengthtlist.append(np.alen(aux))
      poslist.append(aux[:,1])
      widthlist.append(aux[:,2])

    # crop trajectories to equal length
    tcutidx=min(lengthtlist)
    for i in range(numinstances):
      tlist[i]=tlist[i][0:tcutidx]
      poslist[i]=poslist[i][0:tcutidx]
      widthlist[i]=widthlist[i][0:tcutidx]

    # convert into numpy arrays
    tlist=np.asarray(tlist)
    poslist=np.asarray(poslist)
    widthlist=np.asarray(widthlist)
    # tlist should have identical rows - test, but use tlist for readability
    # http://stackoverflow.com/questions/26163727/how-to-test-if-all-rows-are-equal-in-a-numpy
    if not (tlist==tlist[0]).all():
      sys.exit('Exit: problem with non-synchronous trajectories.')
    # fitting range should be within time interval of trajectories analyzed
    for i in range(numinstances):
      if poslist[i,0]>paramanalysis.posstart_fit or poslist[i,-1]<paramanalysis.posend_fit:
        sys.exit('Exit: fitting range extends beyond interval for which data are available.')

    # setting up the plot
    # http://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib
    fig=plt.figure(figsize=(30,15))
    ax1=fig.add_subplot(2,4,1)
    ax2=fig.add_subplot(2,4,2)
    ax3=fig.add_subplot(2,4,3)
    ax4=fig.add_subplot(2,4,4)
    ax5=fig.add_subplot(2,4,5)
    ax6=fig.add_subplot(2,4,6)
    ax7=fig.add_subplot(2,4,7)
    ax8=fig.add_subplot(2,4,8)

    # individual fits
    indfitslopes=np.full(numinstances,np.nan)
    # loop over instances
    for i in range(numinstances):
      # crop region for fit
      my1dcropmask=np.logical_and(paramanalysis.posstart_fit<=poslist[i,:],poslist[i,:]<=paramanalysis.posend_fit)
      croppedt=tlist[i,my1dcropmask]
      croppedpos=poslist[i,my1dcropmask]
      # actual fitting
      auxres=stats.linregress(croppedt,croppedpos)
      myslope=auxres[0]
      myintercept=auxres[1]
      if (not np.isfinite(myslope)) or (not np.isfinite(myintercept)):
        sys.exit('Some problem with fit to individual trajectories.')
      indfitslopes[i]=myslope
      # plotting results
      # the position
      ax1.plot(tlist[i,:],poslist[i,:],'lightblue',alpha=0.2)
      # the fit
      ax1.plot(croppedt,croppedt*myslope+myintercept,'red',alpha=0.5)
      # the speed
      auxlocalvel=np.divide(np.gradient(poslist[i,:]),np.gradient(tlist[i,:]))
      ax3.plot(tlist[i,:],auxlocalvel,'lightblue',alpha=0.2)
      ax3.plot(croppedt,np.mean(auxlocalvel[my1dcropmask])*np.ones_like(croppedt),'red',alpha=0.2)
      # the width
      ax4.plot(tlist[i,:],widthlist[i,:],'lightblue',alpha=0.2)
      ax4.plot(croppedt,np.mean(widthlist[i,my1dcropmask])*np.ones_like(croppedt),'red',alpha=0.2)
    # output statistics
    ax2.hist(indfitslopes)
    slopes_indfit[idxgamma[0],idxrho[0]]=np.mean(indfitslopes)
    errors_indfit[idxgamma[0],idxrho[0]]=stats.sem(indfitslopes)

    # bootstrapping
    if paramanalysis.numbootstrap<=numinstances:
      print('This seems not good: paramanalysis.numbootstrap<=numinstances, just saying...', file=sys.stderr)
    bootstrapfitslopes=np.full(paramanalysis.numbootstrap,np.nan)
    # loop over ensembles of trajectories
    for j in range(paramanalysis.numbootstrap):
      # drawing observed trajectories (instances) with replacement
      obsindex=np.random.randint(0,high=numinstances,size=numinstances)
      # pool instances
      emptyauxensemble=np.full((numinstances,tcutidx),np.nan) # waste of comp time to do this here...
      auxposensemble=emptyauxensemble.copy()
      auxwidthensemble=emptyauxensemble.copy()
      auxfitregionensemble=emptyauxensemble.copy()
      for k,obs in np.ndenumerate(obsindex):
        auxposensemble[k[0],:]=poslist[obs,:]
        auxwidthensemble[k[0],:]=widthlist[obs,:]
        auxfitregionensemble[k[0],:]=np.logical_and(paramanalysis.posstart_fit<=poslist[obs,:],poslist[obs,:]<=paramanalysis.posend_fit)
      # compute mean
      meanposlist=np.mean(auxposensemble,axis=0)
      meanwidthlist=np.mean(auxwidthensemble,axis=0)
      meantlist=np.mean(tlist,axis=0) # just for convenience...
      # crop region for fit
      my1dcropmask=np.all(auxfitregionensemble,axis=0)
      croppedmeantlist=meantlist[my1dcropmask]
      croppedmeanposlist=meanposlist[my1dcropmask]
      if croppedt.shape[0]<30:
        print('This seems not good: fewer than 30 data points for fit, just saying...', file=sys.stderr)
      # fit
      auxres=stats.linregress(croppedmeantlist,croppedmeanposlist)
      myslope=auxres[0]
      myintercept=auxres[1]
      if (not np.isfinite(myslope)) or (not np.isfinite(myintercept)):
        sys.exit('Some problem with fit to ensembles of trajectories.')
      bootstrapfitslopes[j]=myslope
      # plotting results
      # ensemble mean of positions
      ax5.plot(meantlist,meanposlist,'blue',alpha=0.2)
      # fits to ensemble mean of positions
      ax5.plot(croppedmeantlist,croppedmeantlist*myslope+myintercept,'red',alpha=0.5)
      # speed in ensemble mean of positions
      auxlocalmeanvel=np.divide(np.gradient(meanposlist),np.gradient(meantlist))
      ax7.plot(meantlist,auxlocalmeanvel,'blue',alpha=0.2)
      ax7.plot(croppedmeantlist,np.mean(auxlocalmeanvel[my1dcropmask])*np.ones_like(croppedmeantlist),'red',alpha=0.2)
      # ensemble mean of widths
      ax8.plot(meantlist,meanwidthlist,'blue',alpha=0.2)
      ax8.plot(croppedmeantlist,np.mean(meanwidthlist[my1dcropmask])*np.ones_like(croppedmeantlist),'red',alpha=0.2)
    # output statistics
    ax6.hist(bootstrapfitslopes)
    slopes_bootstrap[idxgamma[0],idxrho[0]]=np.mean(bootstrapfitslopes)
    errors_bootstrap[idxgamma[0],idxrho[0]]=np.std(bootstrapfitslopes,ddof=1)

    # save figure
    pars={'bbox_inches':"tight",'dpi':400}
    fig.savefig(pngfolder+"/gamma{0:.6f}_rho{1:.6f}.png".format(gamma,rho),**pars)
    plt.close(fig)

# output into text files
# one file for each gamma
for idxgamma,gamma in np.ndenumerate(param.gammas):
  foutgamma=open('OUTPUT/'+paramanalysis.anfolder+'/rho_speed_error_gamma{0:.6f}.dat'.format(gamma),'w')
  for idxrho,rho in np.ndenumerate(param.rhos):
    print(rho, slopes_indfit[idxgamma[0],idxrho[0]], errors_indfit[idxgamma[0],idxrho[0]], slopes_bootstrap[idxgamma[0],idxrho[0]], errors_bootstrap[idxgamma[0],idxrho[0]], file=foutgamma)
  foutgamma.close()

# ... and now one file for each rho
for idxrho,rho in np.ndenumerate(param.rhos):
  foutrho=open('OUTPUT/'+paramanalysis.anfolder+'/gamma_speed_error_rho{0:.6f}.dat'.format(rho),'w')
  for idxgamma,gamma in np.ndenumerate(param.gammas):
    print(gamma, slopes_indfit[idxgamma[0],idxrho[0]], errors_indfit[idxgamma[0],idxrho[0]], slopes_bootstrap[idxgamma[0],idxrho[0]], errors_bootstrap[idxgamma[0],idxrho[0]], file=foutrho)
  foutrho.close()
