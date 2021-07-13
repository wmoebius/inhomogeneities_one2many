import numpy as np
import pandas as pd
import os
import sys
sys.path.insert(1,os.path.abspath('..'))
import fmmfronthelper

X=np.loadtxt('TEMP/myX.txt')
Y=np.loadtxt('TEMP/myY.txt')
speedmatrix=np.loadtxt('TEMP/myspeedmatrix.txt')
phi=np.ones_like(X)
phi[:,0]=0

fmmfronthelper.fmmandfront(X, Y, phi, speedmatrix, (1,0), 0.25, 0, 'TEMP', 1, True)

falldata=pd.HDFStore('TEMP/alldata.me1.h5')
falldata.allfronts.to_csv('TEMP/fronts.dat',sep=' ',header=False,index=False)
falldata.speedmatrix.to_csv('TEMP/speedmatrix.dat',sep=' ',header=False,index=False)
falldata.nanfilledtm.to_csv('TEMP/nanfilledtm.dat',sep=' ',header=False,index=False)
falldata.close()
