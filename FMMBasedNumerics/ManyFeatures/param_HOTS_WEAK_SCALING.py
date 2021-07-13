import numpy as np

# parameters

# setting the stage
widthch=50.0 # float - width of channel - along y
lengthch=1300.0 # float - length of channel - along x - in propagation direction
neglengthch=50.0 # float - length of channel before centers of features - along x - in propagation direction
invdx=15.0 # float - discretization parameter: number of lattice sites per unit length

# properties of individual features (hotspots/obstacles)
radiusfeat=np.sqrt(1.5) # float - radius of circular feature or in case of ellipses the semi-axis perpendicular to propagation
ratiofeat=1.0 # float - ratio: axis in the direction of propagation / axis perpendicular to propagation
vbasic=1.0 # float - 'propagation speed scale', see main function

# parameters characterizing the heterogeneous environment
# numpy arrays with more than one entry -> main program will loop
gammas=np.array([1.0005,1.0015,1.0045,1.0135,1.0405,1.1215,1.3645,2.0935]) # float - ratio: speed inside / outside of features
gammas_shakeup=np.array([1.005,1.01,1.03,1.05,1.08,1.0,1.0,1.0]) # float - ratio: speed inside / outside of features
                               # features before ROI, to shake up the front
                               # needs to be same length as array gammas
rhos=np.array([np.log(2)/np.pi/1.5]) # float - number density of features

if gammas.shape!=gammas_shakeup.shape:
      sys.exit('Error: Arrays gammas and gammas_shakeup do not have same shape.')

# illustration of metamaterial properties
# affects density of features and boundary condition
# integer please; 0 means no, 1 means yes
refraction=0

# flags for defining the front
# 0 means front as maxpos at fixed yvalues, 1 means front defined from contour finding
contour_front=0

# flags for dumping files
# for dumping this information: 1 means yes, 2 means yes if also me=1
dump_data=1
