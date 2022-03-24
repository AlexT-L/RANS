""" This module calculates turbulent viscosity at the cell faces.

    Libraries/Modules:
        numpy\n
        """
import numpy as np
# from Grid import Grid

# class BaldwinLomax():
# @profile
def turbulent_viscosity(model, ws, state):
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
    # necessary fields
    def wsget(varName):
        return ws.get_field(varName, model.className)
    def mget(varName):
        return model.params[varName]
    gamma = mget('gamma')
    rm = mget('rm')
    re = mget('re')
    ncyc = 0
    rev = wsget('ev')

    il = nx+1
    jl = ny+1
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

    # tauw = np.ones(nx)
    # yscal = np.ones(nx)
    # vor = np.ones((nx,ny))
    # avorm = np.ones(nx)
    # ravg = np.ones(nx)
    amu = np.ones((nxp,nyp))
    amut = np.ones((nxp,nyp))
    # amuto = np.ones(nx)
    # amuti = np.ones(nx)
    # yvor = np.ones(nx)
    # yvorm = np.ones(nx)
    # utot = np.ones((nx,ny))
    # utotm = np.ones(nx)
    # utmin = np.ones(nx)
    # fkleb = np.ones(nx)
    # jedge = np.ones(nx)
    # utmax = np.ones(nx)
    u = np.ones((nxp,nyp))
    v = np.ones((nxp,nyp))
    t = np.ones((nxp,nyp))
    # fcros = np.ones(nx)
    rinv = np.ones((nxp,nyp))
    # vola = np.ones((nx,ny))
    # ylen = np.ones((nx,ny))
    # ylenm = np.ones(nx)

    j2        = je
    jlm       = jl- 1

    jstop     = 3* (j2- 2)/5
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

    rinv[:,:] = 1.0/w[:,:,0]
    t[:,:]    = p[:,:]* rinv[:,:]
    u[:,:]    = w[:,:,1]* rinv[:,:]
    v[:,:]    = w[:,:,2]* rinv[:,:]
    amu[:,:]  = t[:,:]
    amut[:,:] = 0.0

    '''
    Determination of eddy viscosity
    Turbulence model:
    Wall boundary layer --- baldwin-lomax model
    Wake region         --- baldwin-lomax model (cwake= 1.0)
    Calculates vorticity and total velocity
    '''

    vola[:,:] = 0.5* (vol[:,:]+ vol[:,:+1])

    xxa = x[1:,1,1]-x[:il,1,1] 
    yxa = x[1:,1,2]-x[:il,1,2]
    uy        = u[1:,2]- u[1:,1]
    vy        = v[1:,2]- v[1:,1]
    uavg      = .5* (u[1:,1]+u[1:,2])
    vavg      = .5* (v[1:,1]+v[1:,2])


    vor1      = (xxa*uy + yxa*vy)/vola[1:,1]
    vor2      = 0.0
    vor3      = 0.0

    vort      = vor1-vor2-vor3
    vor[1:,1]  = abs(vort)
    utotal    = uavg*uavg + vavg*vavg

    utot[1:,1] = np.sqrt(utotal)


    xxa       = x[1:il,1:jlm,0]-x[0:il-1,1:jlm,0] 
    yxa       = x[1:il,1:jlm,1]-x[0:il-1,1:jlm,1] 
    uy        = u[1:il,2:jlm+1]- u[1:il,1:jlm]
    vy        = v[1:il,2:jlm+1]- v[1:il,1:jlm]
    uavg      = 0.5* (u[1:il,1:jlm]+ u[1:il,2:jlm+1])
    vavg      = 0.5* (v[1:il,1:jlm]+ v[1:il,2:jlm+1])
    '''
    thin-layer navier-stokes contribution to vorticity
    '''
    vor1      = (xxa*uy + yxa*vy)/vola[1:il,1:jlm]
    '''
    additional contributions to vorticity
    '''
    xyw       = 0.5* (x[0:il-1,2:jlm+1,0]- x[0:il-1,0:jlm-1,0])
    xye       = 0.5* (x[1:il,2:jlm+1,0]- x[1:il,0:jlm-1,0])
    yyw       = 0.5* (x[0:il-1,2:jlm+1,1]- x[0:il-1,0:jlm-1,1])
    yye       = 0.5* (x[1:il,2:jlm+1,1]- x[1:il,0:jlm-1,1])
    volawi    = 2.0/(vola[1:il,1:jlm]+ vola[0:il-1,1:jlm])
    volaei    = 2.0/(vola[1:il,1:jlm]+ vola[2:il+1,1:jlm])
    uxe       = 0.5* (u[2:il+1,1:jlm]+u[2:il+1,2:jlm+1]) - uavg
    vxe       = 0.5* (v[2:il+1,1:jlm]+v[2:il+1,2:jlm+1]) - vavg
    uxw       = uavg - 0.5* (u[0:il-1,1:jlm]+u[0:il-1,2:jlm+1])
    vxw       = vavg - 0.5* (v[0:il-1,1:jlm]+v[0:il-1,2:jlm+1])
    vor2      = 0.5* (xye* volaei* uxe+ xyw* volawi* uxw)
    vor3      = 0.5* (yye* volaei* vxe+ yyw* volawi* vxw)
    vort      = vor1- vor2- vor3

    vor[1:il,1:jlm]  = abs(vort)
    utotal    = uavg* uavg+ vavg* vavg
    utot[1:il,1:jlm] = np.sqrt(utotal)

    '''
    Determine transition index
    '''
    itr1      = 0
    itr2      = 0
    j         = 1
    for i in range(0,il):
        if (x[i,j,1] <= xtran):
            itr1      = i - 1
            break 

    itr1p     = itr1 + 1
    for i in range(itr1p-1,il):
        if (x[i,j,1] >= xtran):
            itr2      = i
            break

    for i in range(1,il):
        avor   = vor[i,0:jlm]
        utot1  = utot[i,0:jlm]

        jmaxv     = np.argmax(avor) 
        if (jmaxv == 0): 
            jmaxv = 1 

        jminut    = np.argmin(utot1)
        jmaxut    = np.argmax(utot1)

        avorm[i]  = avor[jmaxv]
        utmin[i]  = utot1[jminut]
        utmax[i]  = max(utot1[jmaxut],1.e-3)
        utotm[i]  = utmax[i]
        yscal[i]  = 1000000.

    if (modbl == 1):
        tur1    = 1.0
        tur2    = 0.
        tur3    = 0.
    elif (modbl == 2):
        tur1    = 0.
        tur2    = 1.0
        tur3    = 0.
    else:
        tur1    = 0.
        tur2    = 0.
        tur3    = 1.0

    xxa       = x[itlp-1:itu,0,0]- x[itlp-1:itu-1,0,0]
    yxa       = x[itlp-1:itu,0,1]- x[itlp-1:itu-1,0,1]
    volai     = 1.0/vola[itlp-1:itu,0]
    uy        = 2.0* u[itlp-1:itu,1]
    amub      = .5* (amu[itlp-1:itu,0]+ amu[itlp-1:itu,1])
    tauw[itlp-1:itu]   = amub* (xxa* uy)* volai
    avor1     = vor[itlp-1:itu,0]
    avora     = avor1
    avorb     = 0.5* (avor1+ avorm[itlp-1:itu])
    avorc     = avorm[itlp-1:itu]
    avor2     = tur1*avora + tur2*avorb + tur3*avorc
    yscal[itlp-1:itu]  = np.sqrt(rey* sgrmi* amub* avor2* w[itlp-1:itu,1,0])/(26.*amub)
    '''
    Compute normal distance ylen[i,j] and function 'yvor'
    (yvor = y* vorticity)
    '''
    ylen[1:il,1] = 0.0
    xc2       = .50* (x[1:il,1:jlm,0]+ x[0:il-1,1:jlm,0]-x[1:il,0:jlm-1,0]- x[0:il-1,0:jlm-1,0])
    yc2       = .50* (x[1:il,1:jlm,1]+ x[0:il-1,1:jlm,1]-x[1:il,0:jlm-1,1]- x[0:il-1,0:jlm-1,1])

    scalf  = np.sqrt(np.square(xc2) + np.square(yc2))
    ylen[1:il,1:jlm] = ylen[1:il,0:jlm-1]+ scalf

    for i in range(1,il):
        ylen1     = 0.5* ylen[i,1]
        for j in range(0,int(np.floor(jstop))):
            y1        = yscal[i]* ylen[i,j]
            damp      = 1.0- np.exp(-y1)
            yvor[j]   = ylen[i,j]* vor[i,j]* damp

        jmaxyv    = np.argmax(yvor)
        jmaxyv    = max(jmaxyv,2)
        jedge[i]  = jmaxyv

        yvorm[i]  = max(yvor[jmaxyv],1.e-6)
        ylenm[i]  = max(ylen[i,jmaxyv],ylen1)

        if (jedge[i] < int(np.floor(jstop))):
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
    for i in range(1,il): 
        udiff     = abs(utmax[i]- utmin[i])
        udiff1    = cwk1* udiff
        for j in range(1,int(np.floor(jstop))): 
            ravg[j]   = 0.5* (w[i,j,0]+ w[i,j+1,0])
            coeff     = 0.0168* ccp
            fwake1    = coeff* yvorm[i]* ylenm[i]
            coeff2    = coeff* cwk1* cwk1
            fwake2    = coeff2* ylenm[i]* udiff* udiff/yvorm[i]
            fwake     = min(fwake1,fwake2)
            fkleb0    = ckleb* ylen[i,j]/ylenm[i]
            fkleb1    = min(fkleb0,1.e5)
            fkleb[j]  = 1.0/(1.0+ 5.5* fkleb1**6)
            amuto[j]  = rey* sgrmi* ravg[j]* fwake* fkleb[j]
            amuto[j]  = abs(amuto[j])
        amuto[1]  = amuto[2]

        '''
        Compute inner eddy viscosity
        '''
        j=int(np.floor(jstop))
        y1        = yscal[i]* ylen[i,1:j]
        damp      = 1.0- np.exp(-y1)
        tscali    = 0.4* ylen[i,1:j]* damp
        amuti1    = tscali* tscali* vor[i,1:j]
        amuti[1:j]   = rey* sgrmi* ravg[1:j]* amuti1
        amuti[1:j]   = abs(amuti[1:j])
        amuti[1]  = 0.0
        if (i<=itl or i>itu):
            amuti[1] = amuti[2]

        '''
        Load viscosity coeffs. into array, use inner value until
        match point is reached
        '''
        ivect     = 1
        if (ivect == 0): 
            icross    = 0
            amut[i,0] = amuti[0]
            for j in range(1,int(np.floor(jstop))): 
                if (amuti[j]<=amuto[j] and icross==0): 
                    amut[i,j] = amuti[j]
                else:
                    icross    = 1
                    amut[i,j] = amuto[j]
        else: 
            amut[i,0] = amuti[0]
            ystop     = int(np.floor(jstop))
            for j in range(0,int(np.floor(jstop))): 
                amudif    = amuti[j]- amuto[j]
                if (amudif >= 0): 
                    fcros[j] = j
                else: 
                    fcros[j] = 1000
            jcros    = np.argmin(fcros)
            if (jcros == 1):
                jcros = 2
            jcrosm   = jcros- 1
            amut[i,0:jcrosm] = amuti[0:jcrosm]
            j = int(np.floor(jstop))
            amut[i,jcros-1:j] = amuto[jcros-1:j]
        '''
        Compute turbulent viscosity at cell center
        '''
        j=int(np.floor(jstop))
        amutc     = 0.5* (amut[i,1:j]+ amut[i,0:j-1])
        amu[i,1:j]  = amutc
        amut[i,1:j] = amu[i,1:j]

        if (i>itr1 and i<=itr2):
            amut[i,1:j] = 0.

    '''
    Copy amut to rev
    '''
    scale     = 1.
    rev[0:ie,0:je]  = scale*amut[0:ie,0:je]

    return


