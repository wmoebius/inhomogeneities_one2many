from __future__ import division
from __future__ import print_function
from builtins import str
from past.utils import old_div
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import subprocess
import param # parameters from main simulation
import paramanalysis # parameters from analysis (for determining region of interest)

def fromh5topngandmp4():

  # read simulation data
  falldata=pd.HDFStore(inputfolder+'/'+inputfile)
  xfm=np.squeeze(falldata.xformatrix.values)
  yfm=np.squeeze(falldata.yformatrix.values)
  allfr=falldata.allfronts
  auxnumdig=len(str(int(allfr['counter'].max()))) # number of digits needed for file names
  fronts=allfr.groupby('counter')

  # region of interest depends on whether we have a uniform channel or the refraction case
  if param.refraction==1:
    xmin=np.amin(xfm)
    xmax=np.amax(xfm)
  else:
    xmin=paramanalysis.posstart_fit
    xmax=paramanalysis.posend_fit
  ymin=np.amin(yfm)
  ymax=np.amax(yfm)

  # 'crop' matrices and handle nans (see fmmfronthelper.py)
  sm=falldata.speedmatrix.values[:,np.logical_and(xfm>=xmin,xfm<=xmax)] # speed matrix (to plot environment)
  tm=falldata.nanfilledtm.values[:,np.logical_and(xfm>=xmin,xfm<=xmax)] # matrix of arrival times (to plot occupied region)
  auxtm=tm.copy()
  auxtm[np.isnan(tm)]=-1.0
  falldata.close()

  # create a figure to be used for all frames
  fig=plt.figure(figsize=[19.2/2,old_div(19.2/2*(1.0*ymax-ymin),(1.0*xmax-xmin))])
  startframe=np.nan # we will need the number of the first frame as an argument to ffmpeg later, and to set frame rate
  # loop over different frames
  for aux,front in fronts:
    valsfront=front.values
    # only relevant frames
    if xmin<=np.nanmean(valsfront[:,2]) and np.nanmean(valsfront[:,2])<=xmax:
      if np.isnan(startframe):
        startframe=int(aux)
      # add axes 
      ax=fig.add_axes([0,0,1,1])
      ax.set_xlim(xmin-0.5/param.invdx,xmax+0.5/param.invdx)
      ax.set_ylim(ymin-0.5/param.invdx,ymax+0.5/param.invdx)
      ax.set_aspect('equal')
      # the environment
      ax.imshow(np.log(old_div(sm,param.vbasic)+1.0e-17),extent=[xmin-0.5/param.invdx,xmax+0.5/param.invdx,ymin-0.5/param.invdx,ymax+0.5/param.invdx],origin='lower',vmin=-3.0,vmax=3.0,cmap=plt.cm.bwr)
      # the region occupied by the population
      auxpop=np.ones_like(sm)
      auxpop[valsfront[0,1]<auxtm]=0 # population has not reached a given lattice point, see fmmfronthelper.py
      auxpop[auxtm<0]=0 # no population where there were nans in tm
      ax.imshow(auxpop,extent=[xmin-0.5/param.invdx,xmax+0.5/param.invdx,ymin-0.5/param.invdx,ymax+0.5/param.invdx],origin='lower',cmap='Greys',alpha=0.2)
      # the front
      ax.plot(valsfront[:,2],valsfront[:,3],marker='.',color='black',markersize=0.5,linestyle='') # not sure why linestyle is required
      # save and clear figure
      fig.savefig(outputfolder+'/'+outputsubfolder+'/'+str(int(aux)).zfill(auxnumdig)+'.png',bbox_inches='tight')
      fig.clf()
      endframe=int(aux) # to set frame rate below
  plt.close()

  # create video from individual frames
  # https://hamelot.io/visualization/using-ffmpeg-to-convert-a-set-of-images-into-a-video/
  # https://stackoverflow.com/questions/20847674/ffmpeg-libx264-height-not-divisible-by-2
  # https://superuser.com/questions/326629/how-can-i-make-ffmpeg-be-quieter-less-verbose
  ffmpegcmd='ffmpeg -loglevel panic -r '+str(int(np.round((endframe-startframe)/20.0)))+' -start_number '+str(startframe)+' -i '+outputfolder+'/'+outputsubfolder+'/%0'+str(auxnumdig)+'d.png -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2:color=white" -vcodec libx264 -pix_fmt yuv420p '+outputfolder+'/'+outputvideofn
  print(ffmpegcmd)
  # https://www.journaldev.com/16140/python-system-command-os-subprocess-call
  returned_value=subprocess.call(ffmpegcmd,shell=True)
  if returned_value!=0:
    sys.exit('Problem with creating video using ffmpeg.')

#------------
# main script
#------------

# create output folder(s) if necessary
outputfolder='OUTPUT/Videos'
if not os.path.exists(outputfolder):
  os.makedirs(outputfolder)

# three input arguments
print('Note: three input parameters, index for rho, index for gamma, the number of the realization (me).')
if not len(sys.argv)==4:
   sys.exit('Wrong number of input arguments.')
rhoid=int(sys.argv[1])
gammaid=int(sys.argv[2])
me=int(sys.argv[3])
print('rhoid:',rhoid)
print('gammaid:',gammaid)
print('me:',me)

# loop over parameters number density
for rho in param.rhos[rhoid:rhoid+1]:
  # loop over parameter hotspot / obstacle intensity
  for gamma in param.gammas[gammaid:gammaid+1]:
    inputfolder='OUTPUT/data.gamma.{0:.6f}/rho.{1:.6f}'.format(gamma,rho)
    inputfilewpath=inputfolder+'/alldata.me{0:d}.h5'.format(me)
    # just input file
    inputfile=os.path.basename(inputfilewpath)
    # create output subfolder and video file names, include parameters
    outputsubfolder=inputfile
    outputsubfolder=outputsubfolder.replace('alldata','png.gamma.{0:.6f}.rho.{1:.6f}'.format(gamma,rho))
    outputsubfolder=outputsubfolder.replace('.h5','')
    outputvideofn=outputsubfolder.replace('png.','video.')+'.mp4'
    if not os.path.exists(outputfolder+'/'+outputsubfolder):
      os.makedirs(outputfolder+'/'+outputsubfolder)
    # create frames and video 
    print(inputfilewpath)
    fromh5topngandmp4()
