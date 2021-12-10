import numpy as np 
import pandas as pd

#####################################################################
# Parameters: The various input params are seperated into dictionaries
######################################################################
# dims:         nx = number of cells in i direction
#               ny = number of cells in j direction
#
# solv_param:   fcyc      = the number of cycles
#               fprnt     = the interval at which the solution is printed
#               fout      = the interval at which convergence is monitored
#               ftim      = the interval at which the permissible time step
#                            is calculated (1 for multigrid calculations)
#               gprnt       controls the initial printout
#               gprnt     = 0 gives no initial printout
#               gprnt     = 1 gives a printout of the mesh
#               gprnt     = 2 gives printouts of the mesh and the initial flow
#               hprnt     = the interval used in printing the solution
#               hmesh     = the number of meshes used in the multigrid sequence
#
#               cflf      = the courant number for the time step on the fine mesh
#                           (cflf<0 selects the use of a variable local step)
#               cflim     = ?
#               bc        = optional  far field boundary conditions
#               vis2      = the coefficient for the adaptive dissipation
#               vis4      = the coefficient for the background dissipation
#               adis      = exponent for directional scaling of the dissipation
#                           (adis = 1. for isotropic dissipation)
#               hmf       = the enthalpy damping factor for the fine mesh
#               cstp(k)   = fraction of the time step at each stage
#               cdis(k)   = fraction of the dissipative terms to be replaced
#                           (no evaluation if cdis(k).eq.0.)
#               mstage    = number of stages in the integration scheme
#               smoopi    = smoothing coefficient for the i direction
#               smoopj    = smoothing coefficient for the j direction
#               ksmoop      controls the use of residual averaging:
#               ksmoop    = 0. for no residual averaging
#               ksmoop    = 1. for residual averaging at all stages
#               ksmoop    = -1. for residual averaging at alternate stages
#               vt        = local time step (1=same local step, 0=variable local step)
#               iprec     = ?
#               epsf      = ?
#               epsc      = ?
#               diag      = ?
#               cflc      = the courant number for time steps on the coarse meshes
#               hmc       = the enthalpy damping factor for the coarse meshes
#               fbc         controls the far field boundary condition
#               fbc       = 0. to freeze the far field on the coarse meshes
#               fbc       = 1. to update the far field on the coarse meshes
#               fcoll     = a relaxation factor for the collected residuals
#               fadd        controls the smoothing of the interpolated corrections
#               vis0      = the dissipative coefficient for the coarse meshes
#               lcyc        controls the multigrid cycle
#               lcyc      = 1. for a v cycle
#               lcyc      = 2. for a w cycle
#
# flo_param:    rm        = the mach number
#               al        = the angle of attack in degrees
#               alpha     = the angle of attack in radians
#               ca        = cos(alpha)
#               sa        = sin(alpha)
#               fcl         controls the option to fix the lift coefficient
#               fcl       = 0. for fixed angle of attack
#               fcl       = 1. for fixed lift coefficient
#               clt       = the target lift coefficient
#               cd0       = the parasite drag coefficient
#               re        = the reynolds number
#               prn       = the prandtl number
#               prt       = the turbulent prandtl number
#               xtran     = the transition point
#               t0        = stagnation temperature of the free-stream
#               kvis        selects the mathematical model
#               kvis      = 0 for inviscid flow
#               kvis      = 1  for laminar flow
#               kvis      = 2  for turbulent flow
#               gamma     = ratio of specific heats (C_p & C_v) for air at STP
#               rho0      = density of free-stream
#               p0        = pressure of the free-stream
#               c0        = speed of sound for the the free-stream
#               ei0       = ???
#               u0        = x-velocity for the free-stream
#               v0        = y-velocity for the free-stream
#               h0        = enthalpy for the free-stream
#               mu0       = kinematic viscosity of the free-stream
# 
# geo_param:    boundx    = ??
#               boundy    = ??
#               bunch     = ??
#               xte       = ??
#               ylim1     = ??
#               ylim2     = ??
#               ax        = ??
#               ay        = ??
#               sy        = ??
#               aplus     = ??
#               ncut      = ??
#               isym      tells us if the airfoil is symmetric about the x-axis
#               isym      = 
#               nu        = total number of points on the upper airfoil surface
#               nl        = total number of points on the lower airfoil surface
#               nn        = total number of points on the airfoil surface
#               trail     = trailing edge angle in degrees
#               slopt     = trailing edge slope
#               xsing     = x-coordinate in physical space that 
#                           becomes origin after square root mapping
#               ysing     = y-coordinate in physical space that 
#                           becomes origin after square root mapping

# in_var:       xn        =x-coordinate of airfoil geometry in physical space
#               yn        =y-coordinate of airfoil geometry in physical space
class Input:

    dim_p=[["nx","ny"]]
    solv_p=[["fcyc","fprnt","fout","ftim","gprnt","hprnt","hmesh"],
           ["cflf","cflim","vis2","vis4","adis","qdis","bc","hmf"],
            ["cstp"],["cdis"],["smoopi","smoopj","ksmoop","vt"],
            ["iprec","epsf","epsc","diag"],
            ["cflc","fcoll","fadd","vis0","hmc","fbc","lcyc"]]
    flo_p=[["rm","al","fcl","clt","cd0"],["re","prn","prt","t0","xtran","kvis"]]
    geo_p=[["boundx","boundy","bunch"],["xte","ylim1","ylim2","ax","ay","sy"],
           ["aplus","ncut"],["isym","nu","nl"],["trail","slopt","xsing","ysing"]]
    in_v=["xn","yn"]

    
    # Constructor
    def __init__(self, filename):
        #Reading in file
        self.max_cols=0
        self.df=pd.DataFrame()
        self.max_no_cols(filename)
        self.read(filename,self.max_cols)
        
        #Param dictionaries
        self.dims={}
        self.solv_param={}
        self.flo_param={}
        self.geo_param={}
        self.in_var={}
        
        #Updating dictionaries

        #dims
        self.update_dict(self.df,self.dims,self.dim_p,2)
    
        #solv_param
        self.update_dict(self.df,self.solv_param,self.solv_p,4)

        cstp=self.solv_param["cstp"]
        cdis=self.solv_param["cdis"]

        self.solv_param["mstage"]=len(cstp)

        #making cdis and cstp on length 6 by adding zeros
        self.solv_param["cdis"]=np.append(cdis,np.zeros(6-len(cdis)))
        self.solv_param["cstp"]=np.append(cstp,np.zeros(6-len(cstp)))

        if self.solv_param["cflc"]==0.0:
            self.solv_param["cflc"]=self.solv_param["cflf"]
        
        #flow_param
        self.update_dict(self.df,self.flo_param,self.flo_p,18)
        flo=self.flo_param
        t0=flo["t0"]
        rm=flo["rm"]
        alpha = flo["alpha"] = flo["al"]*np.pi/180 #convert to radians
        ca = flo["ca"] = np.cos(flo["alpha"])
        sa = flo["sa"] = np.sin(flo["alpha"])
        
        kvis=flo["kvis"]
        if kvis > 1:
            flo["re"]=flo["re"]*1e-6 
        
        #set constants and far-field values
        gamma = flo["gamma"] = 1.4
        rho0 = flo["rho0"]  = 1.
        p0 = flo["p0"]    = 1.
        c0 = flo["c0"]    = np.sqrt(gamma*p0/rho0)
        ei0 = flo["ei0"]  = p0/((gamma  -1.)*rho0)
        u0 = flo["u0"]    = rm*c0*ca
        v0 = flo["v0"]    = rm*c0*sa
        h0 = flo["h0"]    = gamma*ei0  +.5*(u0*u0  +v0*v0)
        mu_air= 1.461e-06 #kinematic viscosity of air at sea-level at STP
        mu0 = flo["mu0"]  = mu_air*((t0)**(3/2))/(t0+110.3)

        #geoparam
        self.update_dict(self.df,self.geo_param,self.geo_p,22)
        geo=self.geo_param
        geo["nn"]=geo["nu"] + geo["nl"] -1

        #airfoil coordinated (in_var)
        self.update_geom(self.df,self.in_var,self.in_v,32)


    #Methods

    #Get max number of columns in a row
    def max_no_cols(self,file):
        #Loop the data lines
        with open(file, 'r') as temp_f:
            # get No of columns in each line
            col_count = [ len(l.split()) for l in temp_f.readlines() ]
        self.max_cols=max(col_count)
        return 
    
    # Read .data file
    def read(self,file,max_cols):
        dfs = pd.read_csv(file, header=None, delimiter="\s+|;|:", names=range(max_cols),engine="python")
        #convert string to float
        self.df=dfs.apply(pd.to_numeric, errors='coerce')
    
        return

    # Update param dictionaries
    def update_dict(self,df,dict,params,strt_row):
        no_nan=np.array(df.count(axis=1))
        for i in range(len(params)):
            row=strt_row +2*i
            row_len=len(params[i])#lenth of row
            
            if row_len > 1:
                dict.update(zip(params[i],df.iloc[row,0:no_nan[row]]))
            else:
                dict[params[i][0]]=np.array(df.iloc[row,0:no_nan[row]])
        return

    # Update aerfoil geometry
    def update_geom(self,df,dict,params,strt_row):
        no_nan_r=np.array(df.count(axis=1))
        no_nan_c=np.array(df.count())
        for i in range(len(params)):
            upper=np.array(df.iloc[strt_row:96,i])
            lower=np.array(df.iloc[98:162,i])
            dict[params[i]]=np.concatenate((upper,lower),axis=0)
        return

    def add_dicts(dict1,dict2):
        sum_dict= {**dict1, **dict2}
        return sum
     

    

    
##plot airfoil
# import matplotlib.pyplot as plt

# input=Input("rae9e-s3.data")
# x=input.in_var["xn"]
# y=input.in_var["yn"]
# plt.plot(x,y,"+")
# plt.show()
# #print(input.solv_param["cstp"])


