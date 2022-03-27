""" This module calculates turbulent viscosity at the cell faces.

    Libraries/Modules:
        numpy\n
        """
import sys
sys.path.append('../RANS/bin')

import numpy as np
from bin.Field import maximum, minimum, Field
# from Grid import Grid

# class BaldwinLomax():
# @profile
def turbulent_viscosity(model, ws, state,ncyc=0):
    """ Baldwin-lomax turbulence model:  modtur = 2.
    Calculates turbulent viscosity at the cell faces.
    Averages to obtain cell center values fully vectorized routine.                                         *
    Calculates eddy viscosity, vorticity, total velocity, normal distance.
    Also calculates outer and innner eddy viscosity.

    Attributes:
        rev: eddy viscocity
        ylen: normal distance
        vor: vorticity
        vol: control volume
        amuto: outer eddy viscosity
        amuti: inner eddy viscosity

    Notes: 
        Adapted from subroutine turb2.f
    """
    [nx, ny] = ws.field_size()
    PAD = model.padding
    # necessary fields
    def wsget(varName):
        return ws.get_field(varName, model.className)
    def mget(varName):
        return model.params[varName]
    gamma = mget('gamma')
    rm = mget('rm')
    re = mget('re')
    rev = wsget('ev')
    
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [nxp, nyp] = [nx+4, ny+4]
    dims = ws.get_dims()
    itl = dims['itl']
    itu = dims['itu']
    w = state
    x = ws.get_field('x')
    p = wsget('p')

    xtran = mget('xtran') 
    vol = wsget('vol') 

    tauw = np.zeros(nxp)
    yscal = 1000000*np.ones(nxp)
    vor = np.zeros((nxp,nyp))
    avor = np.zeros(nyp)
    avorm = np.zeros(nxp)
    ravg = np.ones(nyp)
    amu = np.ones((nxp,nyp))
    amut = np.ones((nxp,nyp))
    amuto = np.ones(nyp)
    amuti = np.ones(nyp)
    yvor = np.ones(nyp)
    yvorm = np.ones(nxp)
    utot = np.zeros((nxp,nyp))
    utot1 = np.zeros(nyp)
    utotm = np.zeros(nxp)
    utmin = np.zeros(nxp)
    utmax = np.zeros(nxp)
    fkleb = np.ones(nyp)
    jedge = np.ones(nxp)
    u = np.ones((nxp,nyp))
    v = np.ones((nxp,nyp))
    t = np.ones((nxp,nyp))
    fcros = np.ones(nyp)
    rinv = np.ones((nxp,nyp))
    vola = np.ones((nxp,nyp))
    ylen = np.ones((nxp,nyp))
    ylenm = np.ones(nxp)

    i2        = ie
    j2        = je
    jlm       = jl- 1

    jstop     = (int)(3* (j2- 2)/5)
    itlp      = itl+ 1    

    cwk1      = 1.0
    ckleb     = 0.3
    ccp       = 1.6
    modbl     = 3
    restarr   = 0

    rey       = re
    sgam      = np.sqrt(gamma)
    sgrm      = sgam*rm
    sgrmi     = 1.0/sgrm

    rinv[1:i2+1,1:j2+1] = 1.0/w[1:i2+1,1:j2+1,0]
    t[1:i2+1,1:j2+1]    = p[1:i2+1,1:j2+1]* rinv[1:i2+1,1:j2+1]
    u[1:i2+1,1:j2+1]    = w[1:i2+1,1:j2+1,1]* rinv[1:i2+1,1:j2+1]
    v[1:i2+1,1:j2+1]    = w[1:i2+1,1:j2+1,2]* rinv[1:i2+1,1:j2+1]
    amu[1:i2+1,1:j2+1]  = t[1:i2+1,1:j2+1]
    amut[1:i2+1,1:j2+1] = 0.0

    '''
    Determination of eddy viscosity
    Turbulence model:
    Wall boundary layer --- baldwin-lomax model
    Wake region         --- baldwin-lomax model (cwake= 1.0)
    Calculates vorticity and total velocity
    '''
    
    vola[1:i2+1,1:jl+1] = 0.5* (vol[1:i2+1,1:jl+1]+ vol[1:i2+1,2:jl+2])

    xxa = x[1:nx+1,1,0]-x[0:nx,1,0] 
    yxa = x[1:nx+1,1,1]-x[0:nx,1,1]
    uy        = u[2:il+1,2]- u[2:il+1,1]
    vy        = v[2:il+1,2]- v[2:il+1,1]
    uavg      = 0.5* (u[2:il+1,1] + u[2:il+1,2])
    vavg      = 0.5* (v[2:il+1,1] + v[2:il+1,2])


    vor1      = (xxa*uy + yxa*vy)/vola[2:il+1,1]
    vor2      = 0.0
    vor3      = 0.0

    vort      = vor1-vor2-vor3
    vor[2:il+1,1]  = abs(vort)
    utotal    = uavg*uavg + vavg*vavg

    utot[2:il+1,1] = np.sqrt(utotal)


    xxa       = x[1:nx+1,1:jlm,0]-x[0:nx,1:jlm,0] 
    yxa       = x[1:nx+1,1:jlm,1]-x[0:nx,1:jlm,1] 
    uy        = u[2:il+1,3:jlm+2]- u[2:il+1,2:jlm+1]
    vy        = v[2:il+1,3:jlm+2]- v[2:il+1,2:jlm+1]
    uavg      = 0.5* (u[2:il+1,2:jlm+1]+ u[2:il+1,3:jlm+2])
    vavg      = 0.5* (v[2:il+1,2:jlm+1]+ v[2:il+1,3:jlm+2])
    '''
    thin-layer navier-stokes contribution to vorticity
    '''
    vor1      = (xxa*uy + yxa*vy)/vola[2:il+1,2:jlm+1]
    '''
    additional contributions to vorticity
    '''
    xyw       = 0.5* (x[0:nx,2:jlm+1,0]- x[0:nx,0:jlm-1,0])
    xye       = 0.5* (x[1:nx+1,2:jlm+1,0]- x[1:nx+1,0:jlm-1,0])
    yyw       = 0.5* (x[0:nx,2:jlm+1,1]- x[0:nx,0:jlm-1,1])
    yye       = 0.5* (x[1:nx+1,2:jlm+1,1]- x[1:nx+1,0:jlm-1,1])
    volawi    = 2.0/(vola[2:il+1,2:jlm+1]+ vola[1:il,2:jlm+1])
    volaei    = 2.0/(vola[2:il+1,2:jlm+1]+ vola[3:il+2,2:jlm+1])
    uxe       = 0.5* (u[3:il+2,2:jlm+1]+u[3:il+2,3:jlm+2]) - uavg
    vxe       = 0.5* (v[3:il+2,2:jlm+1]+v[3:il+2,3:jlm+2]) - vavg
    uxw       = uavg - 0.5* (u[1:il,2:jlm+1]+u[1:il,3:jlm+2])
    vxw       = vavg - 0.5* (v[1:il,2:jlm+1]+v[1:il,3:jlm+2])
    vor2      = 0.5* (xye* volaei* uxe+ xyw* volawi* uxw)
    vor3      = 0.5* (yye* volaei* vxe+ yyw* volawi* vxw)
    vort      = vor1- vor2- vor3

    vor[2:il+1,2:jlm+1]  = abs(vort)
    utotal    = uavg* uavg+ vavg* vavg
    utot[2:il+1,2:jlm+1] = np.sqrt(utotal)

    '''
    Determine transition index
    '''
    itr1      = 0
    itr2      = 0
    j         = 1
    for i in range(il):
        if (x[i,j,0] <= xtran):
            itr1      = i - 1
            break 

    itr1p     = itr1 + 1
    for i in range(itr1p-1,il):
        if (x[i,j,0] >= xtran):
            itr2      = i
            break

    for i in range(PAD,nx+PAD):
        avor[1:jlm+1]   = vor[i,1:jlm+1]
        utot1[1:jlm+1]  = utot[i,1:jlm+1]

        jmaxv     = np.argmax(avor) 
        if (jmaxv == 0): 
            jmaxv = 1 

        jminut    = np.argmin(utot1)
        jmaxut    = np.argmax(utot1)

        avorm[i]  = avor[jmaxv]
        utmin[i]  = utot1[jminut]
        utmax[i]  = maximum(utot1[jmaxut],1.e-3)
        utotm[i]  = utmax[i]

    if (modbl == 1):
        tur1    = 1.0
        tur2    = 0.0
        tur3    = 0.0
    elif (modbl == 2):
        tur1    = 0.0
        tur2    = 1.0
        tur3    = 0.0
    else:
        tur1    = 0.0
        tur2    = 0.0
        tur3    = 1.0

    xxa       = x[itlp-1:itu-1,0,0]- x[itlp-2:itu-2,0,0]
    yxa       = x[itlp-1:itu-1,0,1]- x[itlp-2:itu-2,0,1]
    volai     = 1.0/vola[itlp:itu,1]
    uy        = 2.0* u[itlp:itu,1]
    amub      = 0.5* (amu[itlp:itu,1]+ amu[itlp:itu,2])
    tauw[itlp:itu]   = amub* (xxa* uy)* volai
    avor1     = vor[itlp:itu,1]
    avora     = avor1
    avorb     = 0.5* (avor1+ avorm[itlp:itu])
    avorc     = avorm[itlp:itu]
    avor2     = tur1*avora + tur2*avorb + tur3*avorc
    yscal[itlp:itu]  = np.sqrt(rey* sgrmi* amub* avor2* w[itlp:itu,2,0])/(26.0*amub)
    '''
    Compute normal distance ylen[i,j] and function 'yvor'
    (yvor = y* vorticity)
    '''
    ylen[2:nx+2,1] = 0.0
    xc2       = 0.5* (x[1:nx+1,1:ny+1,0]+ x[0:nx,1:ny+1,0]-x[1:nx+1,0:ny,0]- x[0:nx,0:ny,0])
    yc2       = 0.5* (x[1:nx+1,1:ny+1,1]+ x[0:nx,1:ny+1,1]-x[1:nx+1,0:ny,1]- x[0:nx,0:ny,1])

    scalf  = np.sqrt(np.square(xc2) + np.square(yc2))
    ylen[PAD:nx+PAD,PAD:ny+PAD] = ylen[PAD:nx+PAD,1:ny+1]+ scalf

    for i in range(PAD,nx+PAD):
        ylen1     = 0.5* ylen[i,2]
        for j in range(1,jstop+1):
            y1        = yscal[i]* ylen[i,j]
            damp      = 1.0- np.exp(-y1)
            yvor[j]   = ylen[i,j]* vor[i,j]* damp

        jmaxyv    = np.argmax(yvor)
        jmaxyv    = maximum(jmaxyv,2)
        jedge[i]  = jmaxyv

        yvorm[i]  = maximum(yvor[jmaxyv],1.0e-6)
        ylenm[i]  = maximum(ylen[i,jmaxyv],ylen1)

        if (jedge[i] < jstop):
            ylenm1  = ylenm[i]

        if (ncyc>=10 or restarr==1.0):
            jmyv    = int(jedge[i])
            dyvm    = yvor[jmyv]-yvor[jmyv-1]
            dyvp    = yvor[jmyv]-yvor[jmyv+1]

            if (yvor[jmyv-1] < yvor[jmyv+1]):
                ylenm[i] = ylen[i,jmyv]+ 0.5*(ylen[i,jmyv+1]- ylen[i,jmyv])*(1- dyvp/dyvm)
            else:
                if dyvp == 0.0:
                    dyvp = 0.00001
                ylenm[i] = ylen[i,jmyv]- 0.5*(ylen[i,jmyv]- ylen[i,jmyv-1])*(1- dyvm/dyvp)
        
        else:
            ylenm[i]  = ylenm1

    '''
    Compute outer eddy viscosity
    '''
    for i in range(PAD,nx+PAD): 
        udiff     = abs(utmax[i]- utmin[i])
        udiff1    = cwk1* udiff
        for j in range(2,jstop+1): 
            ravg[j]   = 0.5* (w[i,j,0]+ w[i,j+1,0])
            coeff     = 0.0168* ccp
            fwake1    = coeff* yvorm[i]* ylenm[i]
            coeff2    = coeff* cwk1* cwk1
            fwake2    = coeff2* ylenm[i]* udiff* udiff/yvorm[i]
            fwake     = minimum(fwake1,fwake2)
            fkleb0    = ckleb* ylen[i,j]/ylenm[i]
            fkleb1    = minimum(fkleb0,1.0e5)
            fkleb[j]  = 1.0/(1.0+ 5.5* fkleb1**6)
            amuto[j]  = rey* sgrmi* ravg[j]* fwake* fkleb[j]
            amuto[j]  = abs(amuto[j])
        amuto[1]  = amuto[2]

        '''
        Compute inner eddy viscosity
        '''
        y1        = yscal[i]* ylen[i,2:jstop+1]
        damp      = 1.0- np.exp(-y1)
        tscali    = 0.4* ylen[i,1:j]* damp
        amuti1    = tscali* tscali* vor[i,2:jstop+1]
        amuti[2:jstop+1]   = rey* sgrmi* ravg[2:jstop+1]* amuti1
        amuti[2:jstop+1]   = abs(amuti[2:jstop+1])
        amuti[1]  = 0.0
        if (i<itl or i>=itu):
            amuti[1] = amuti[2]

        '''
        Load viscosity coeffs. into array, use inner value until
        match point is reached
        '''
        ivect     = 1
        if (ivect == 0): 
            icross    = 0
            amut[i,1] = amuti[1]
            for j in range(2,jstop+1): 
                if (amuti[j]<=amuto[j] and icross==0): 
                    amut[i,j] = amuti[j]
                else:
                    icross    = 1
                    amut[i,j] = amuto[j]
        else: 
            amut[i,1] = amuti[1]
            ystop     = jstop
            for j in range(1,jstop+1): 
                amudif    = amuti[j]- amuto[j]
                if (amudif >= 0): 
                    fcros[j] = (float)(j)
                else: 
                    fcros[j] = 1000.0
                    
            jcros    = np.argmin(fcros)
            if (jcros == 1):
                jcros = 2
            jcrosm   = jcros- 1
            
            amut[i,1:jcrosm+1] = amuti[1:jcrosm+1]
            amut[i,jcros:jstop+1] = amuto[jcros:jstop+1]
        '''
        Compute turbulent viscosity at cell center
        '''
        amutc     = 0.5* (amut[i,2:jstop+1]+ amut[i,1:jstop])
        amu[i,2:jstop+1]  = amutc
        amut[i,2:jstop+1] = amu[i,2:jstop+1]

        if (i>itr1 and i<=itr2):
            amut[i,2:jstop+1] = 0.0

    '''
    Copy amut to rev
    '''
    scale     = 1.0
    rev[1:ie+1,1:je+1]  = scale*amut[1:ie+1,1:je+1]

    return