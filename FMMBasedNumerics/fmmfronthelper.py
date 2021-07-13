from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import numpy as np
import pandas as pd
import skfmm
import sys
#from skimage import measure

def fmmandfront(myx,myy,myphi,myspeedmatrix,myperarg,mydt,mycontourfront,myfolder,myme,mydumpdata):

  # matrix of x values, matrix of y values, matrix phi, speed matrix (may include value 0 for speed),
  # argument for PBC (such as in skfmm), time step for front finding, option for front finding (only one option, 0, is supported at the moment),
  # output folder, realization / counter for environments, option for dumping data in addition to mean front (True/False)

  # the actual FMM part
  # note that after talking with developer of skfmm decided to leave obstacles as zero speed entries
  # in myspeedmatrix instead of masking myphi
  tmatrix=skfmm.travel_time(myphi,myspeedmatrix,dx=np.asscalar(myx[0,1]-myx[0,0]),periodic=myperarg)
  # if myspeedmatrix contains 0s (or smaller (smaller equal?) machine precision),
  # than tmatrix should be masked array (and vice versa), check for that;
  # use 1e-10 instead of machine precision because that's a good enough test for our choice of parameters;
  if hasattr(tmatrix,'mask') and (not np.amin(myspeedmatrix)<1e-10):
    sys.exit('Problem with tmatrix.')
  if (not hasattr(tmatrix,'mask')) and np.amin(myspeedmatrix)<1e-10:
    sys.exit('Problem with treating obstacles 1.')
  # inaccessible lattice sites could result in masked entries in tmatrix, thought the speed is finite there;
  # but all entries with vanishingly small speed should result in masked entries in tmatrix, check for that
  if hasattr(tmatrix,'mask'):
    tmmask=tmatrix.mask
    smsmall=(myspeedmatrix<1e-10)
    if not np.all(smsmall.astype(int)<=tmmask.astype(int)):
      sys.exit('Problem with treating obstacles 2.')

  # find front and output
  myt=0.0 # keep track of time
  timestepcounter=0 # count time steps
  allfronts=[] # all fronts at all times
  meanfront=[] # mean front position (and width) as fcn. of time
  # if tmatrix is a masked array, fill with nans
  if hasattr(tmatrix,'mask'):
    nanfilledtm=tmatrix.filled(np.nan)
  else:
    nanfilledtm=tmatrix

  # when has front reached right side of domain?
  finalt=0.99*np.nanmin(nanfilledtm[:,-1])
  # if that does not give meaningful result, e.g., because
    # example: phi=np.array([[-1,1,1,1],[-1,1,1,1],[-1,1,1,1],[-1,1,1,1]])
    # example: speed=np.ones_like(phi)
    # example: speed[:,2]=0
    # example: tm=skfmm.travel_time(phi,speed)
    # example: auxtm=tm.filled(np.nan)
  # throw an error and and do not go ahead with front finding 
  if not np.isfinite(finalt):
    print('Front has not reached end of domain, no front output,', myfolder,'me',myme, file=sys.stderr)

  # use that again in order to avoid long else block
  while np.isfinite(finalt) and myt<finalt:
    ### definition of front similar to event driven model: max xpos at defined yvalues ###
    if mycontourfront==0:
      # 'advanced' thresholding to find front
      xsreached=myx.copy()
      xsnotreached=myx.copy()
      # use Boolean array indexing
      # because of problem with comparing nan in np.array with a value, auxiliary variable with nans replaced by -1
      auxtm=nanfilledtm.copy()
      auxtm[np.isnan(nanfilledtm)]=-1.0
      xsreached[np.logical_or(myt<auxtm,np.isnan(nanfilledtm))]=np.nan
      maxxsreached=np.nanmax(xsreached,axis=1)
      maxxsreached_idx=np.nanargmax(xsreached,axis=1)
      otherpoint_idx=maxxsreached_idx+1
      # now loop over y to identify front
      # write down the whole front
      # this loop might not be the most Python way to do it...
      frontx=np.full_like(maxxsreached,np.nan)
      for yj in range(0,maxxsreached.size):
        # if next lattice site in x direction is not one to consider
        if otherpoint_idx[yj]>np.argmax(myx) or np.isnan(nanfilledtm[yj,otherpoint_idx[yj]]):
          frontx[yj]=maxxsreached[yj]
        # otherwise interpolate
        else:
          otherpoint=myx[yj,otherpoint_idx[yj]]
          distls=otherpoint-maxxsreached[yj]
          fracprop=old_div((1.0*myt-nanfilledtm[yj,maxxsreached_idx[yj]]),(nanfilledtm[yj,otherpoint_idx[yj]]-nanfilledtm[yj,maxxsreached_idx[yj]]))
          frontx[yj]=maxxsreached[yj]+distls*fracprop
        allfronts.append([timestepcounter, myt, frontx[yj], myy[yj,0]])
      if np.isnan(frontx).any():
        sys.exit('Problem with front detection.')
      # ... and now mean front position and width
      # useful info: https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation, 
      #              http://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.nanstd.html
      meanfront.append([myt,np.nanmean(frontx),np.nanstd(frontx,ddof=1)])
    ### definition of the front using countours finder ###
    elif mycontourfront==1:
      sys.exit('Choice of front definition currently not supported.')
      # this is old non-function code - don't trust, use for inspiration only
      # http://scikit-image.org/docs/dev/auto_examples/edges/plot_contours.html#example-edges-plot-contours-py
      # contours=measure.find_contours(tmatrix,myt)
      # numopencontours=0
      # for thiscon in contours:
      #   thiscon=thiscon/param.invdx
      #   thiscon[:,0]=thiscon[:,0]-param.widthch
      #   xlen=np.max(thiscon[:,0])-np.min(thiscon[:,0])
      #   # comparing doubles should be save here
      #   if any(np.not_equal(thiscon[0,:],thiscon[-1,:])) and (xlen >= param.widthch):
      #     numopencontours=numopencontours+1
      #     if numopencontours > 1:
      #       sys.exit("more than one open contour")
      #     if dump_front==1 or (dump_front==2 and me==1):
      #       np.savez_compressed(folder+"/front.me{0:d}.time".format(me)+str(myt)+".npz",thiscon=thiscon)
      #     np.savetxt(fmf,np.c_[myt,np.mean(thiscon[:,1])],fmt='%.6f')
      #   # if abs(numopencontours-1)>0.1:
      #   #    sys.exit("# of closed contours != 1 - aborting...")
    else:
      sys.exit('Invalid choice for front finding.')
    myt=myt+mydt
    timestepcounter=timestepcounter+1

  # write mean front position as function of time to file 
  if np.isfinite(finalt):
    np.savetxt(myfolder+'/meanfront.me{0:d}.dat'.format(myme),np.asarray(meanfront))

  # dump data (in addition to the above)
  if mydumpdata: 
    fdd=pd.HDFStore(myfolder+'/alldata.me{0:d}.h5'.format(myme),mode='w',complevel=9,complib='bzip2')
    fdd['phi']=pd.DataFrame(myphi)
    fdd['speedmatrix']=pd.DataFrame(myspeedmatrix)
    fdd['nanfilledtm']=pd.DataFrame(nanfilledtm)
    fdd['xformatrix']=pd.DataFrame(np.asarray(myx[0,:]))
    fdd['yformatrix']=pd.DataFrame(np.asarray(myy[:,0]))
    if np.isfinite(finalt):
      fdd['allfronts']=pd.DataFrame(np.asarray(allfronts),columns=['counter','time','x','y'])
      fdd['meanfrontpos']=pd.DataFrame(np.asarray(meanfront),columns=['time','mean pos','std dev pos'])
    fdd.close()
