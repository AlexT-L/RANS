import numpy as np
from Grid import Grid

# class BaldwinLomax():
def turbulent_viscosity(params, dims):
    # from subroutine turb2.f

    #  **********************************************************************
    #  *   baldwin-lomax turbulence model:  modtur = 2                      *
    #  *                                                                    *
    #  *   calculates turbulent viscosity at the cell faces and then        *
    #  *   averages to obtain cell center values                            *
    #  *   fully vectorized routine                                         *
    #  **********************************************************************

    # "uses:" 
    # dims, flo_var, solv_var, mesh_var, psm_var, flo_param, solv_param

    # inputs
    ie = params['ie'] # Mesh dimension
    je = params['je'] # Mesh dimension
    kvis = params['kvis']
    gamma = params['gamma']
    rm = params['rm']
    re = params['re']
    ncyc = params['ncyc']
    rev = params['rev']
    cmesh = params['cmesh'] # doesn't seem to appear anywhere other than here?
    ncyci1 = params['ncyci1'] # also doesn't seem to appear anywhere other than here
    il = dims['il']
    jl = dims['jl']
    itl = params['itl']
    itu = params['itu']
    x = params['x']
    w = params['w']
    p = params['p']
    # w = grid.w
    # p = grid.p
    xtran = params['xtran'] # needs to be from flo_param

    vol = params['vol'] # new
    


    # initializing, defined later
    uedge = []
    dim_var = 100
    tauw = np.ones(dim_var)
    yscal = np.ones(dim_var)
    scalf = np.ones(dim_var)
    vor = np.ones((dim_var,dim_var))
    avor = np.ones(dim_var)
    avorm = np.ones(dim_var)
    ravg = np.ones(dim_var)
    amut = np.ones((dim_var,dim_var))
    amuto = np.ones(dim_var)
    amuti = np.ones(dim_var)
    yvor = np.ones(dim_var)
    yvorm = np.ones(dim_var)
    utot = np.ones((dim_var,dim_var))
    utotm = np.ones(dim_var)
    utot1 = np.ones(dim_var)
    ylenm1 = dim_var
    utmin = np.ones(dim_var)
    fkleb = np.ones(dim_var)
    jedge = np.ones(dim_var)
    utmax = np.ones(dim_var)
    amu = np.ones((dim_var,dim_var))
    u = np.ones((dim_var,dim_var))
    v = np.ones((dim_var,dim_var))
    t = np.ones((dim_var,dim_var))
    fcros = np.ones(dim_var)
    rinv = np.ones((dim_var,dim_var))
    vola = np.ones((dim_var,dim_var))
    ylen = np.ones((dim_var,dim_var))
    ylenm = np.ones(dim_var)

    i2        = ie
    j2        = je
    il        = i2- 1
    jl        = j2- 1
    jlm       = jl- 1

    jstop     = 3* (j2- 2)/5
    if (cmesh < 0.0):
        jstop = jl- 1
    itlp      = itl+ 1
    iwrit     = 6

    aplusi    = 1./26.
    cwk1      = 1.0
    ckleb     = 0.3
    ccp       = 1.6
    modbl     = 3
    restarr   = 0

    rey       = re
    rein      = 1.0/rey
    sgam      = np.sqrt(gamma)
    sgrm      = sgam*rm
    sgrmi     = 1.0/sgrm

# c     **********************************************************************
# c     *    next 2 lines of code is activated to "test" turbulence model    *
# c     *      effect on convergence                                         *
# c     *                                                                    *
# c     *    turbulent viscosity is frozen after ncyct cycles                *
# c     **********************************************************************

    # ncyct     = 10 # commented out in turb2
    # if (ncyc > ncyct) return

    for j in range(0,j2):
        for i in range(0,i2):
            rinv[i,j] = 1.0/w[i,j,0]
            t[i,j]    = p[i,j]* rinv[i,j]
            u[i,j]    = w[i,j,1]* rinv[i,j]
            v[i,j]    = w[i,j,2]* rinv[i,j]
            amu[i,j]  = t[i,j]
            amut[i,j] = 0.0

# c     **********************************************************************
# c     *   determination of eddy viscosity                                  *
# c     *                                                                    *
# c     *   turbulence model:                                                *
# c     *     wall boundary layer --- baldwin-lomax model                    *
# c     *     wake region         --- baldwin-lomax model (cwake= 1.0)       *
# c     *                                                                    *
# c     *   calculate vorticity and total velocity                           *
# c     **********************************************************************

    for j in range(0,jl):
        for i in range (0,i2):
            vola[i,j] = 0.5* (vol[i,j]+ vol[i,j+1])

    for i in range(1,il):
        xxa       = x[i,1,1]-x[i-1,1,1]
        yxa       = x[i,1,2]-x[i-1,1,2]
        uy        = u[i,2]- u[i,1]
        vy        = v[i,2]- v[i,1]
        uavg      = .5* (u[i,1]+u[i,2])
        vavg      = .5* (v[i,1]+v[i,2])
        vor1      = (xxa*uy + yxa*vy)/vola[i,1]
        vor2      = 0.0
        vor3      = 0.0

        vort      = vor1-vor2-vor3
        vor[i,1]  = abs(vort)
        utotal    = uavg*uavg + vavg*vavg
        utot[i,1] = np.sqrt(utotal)

    for j in range(1,jlm):
        for i in range( 1,il):
            xxa       = x[i,j,0]-x[i-1,j,0]
            yxa       = x[i,j,1]-x[i-1,j,1]
            uy        = u[i,j+1]- u[i,j]
            vy        = v[i,j+1]- v[i,j]
            uavg      = 0.5* (u[i,j]+ u[i,j+1])
            vavg      = 0.5* (v[i,j]+ v[i,j+1])

        #  thin-layer navier-stokes contribution to vorticity

            vor1      = (xxa*uy + yxa*vy)/vola[i,j]


        #  additional contributions to vorticity

            xyw       = 0.5* (x[i-1,j+1,0]- x[i-1,j-1,0])
            xye       = 0.5* (x[i,j+1,0]- x[i,j-1,0])
            yyw       = 0.5* (x[i-1,j+1,1]- x[i-1,j-1,1])
            yye       = 0.5* (x[i,j+1,1]- x[i,j-1,1])
            volawi    = 2.0/(vola[i,j]+ vola[i-1,j])
            volaei    = 2.0/(vola[i,j]+ vola[i+1,j])
            uxe       = 0.5* (u[i+1,j]+u[i+1,j+1]) - uavg
            vxe       = 0.5* (v[i+1,j]+v[i+1,j+1]) - vavg
            uxw       = uavg - 0.5* (u[i-1,j]+u[i-1,j+1])
            vxw       = vavg - 0.5* (v[i-1,j]+v[i-1,j+1])
            vor2      = 0.5* (xye* volaei* uxe+ xyw* volawi* uxw)
            vor3      = 0.5* (yye* volaei* vxe+ yyw* volawi* vxw)
            vort      = vor1- vor2- vor3

            vor[i,j]  = abs(vort)
            utotal    = uavg* uavg+ vavg* vavg
            utot[i,j] = np.sqrt(utotal)

    #  determine transition index

    itr1      = 0
    itr2      = 0
    j         = 1
    for i in range(0,il):
        if (x[i,j,1] <= xtran):
            itr1      = i - 1
            break # seems like it might continue, 
                # but if it continues then it changes nothing, so break?

    itr1p     = itr1 + 1
    for i in range(itr1p-1,il):
        if (x[i,j,1] >= xtran):
            itr2      = i
            break

    for i in range(1,il):
        for j in range(0,jlm):
            avor[j]   = vor[i,j]
            utot1[j]  = utot[i,j]
        # j loop ends here

        # effect of using jlm or jstop needs to be checked;
        # it does not seem to make a difference
        # jmaxv     = ismax(jstop,avor,1)
        jmaxv     = np.argmax(avor) # replacing ismax with np.argmax.
        if (jmaxv == 0): 
            jmaxv = 1 # Seems weird to me

        jminut    = np.argmin(utot1)
        jmaxut    = np.argmax(utot1)

        avorm[i]  = avor[jmaxv]
        utmin[i]  = utot1[jminut]
        utmax[i]  = max(utot1[jmaxut],1.e-3)
        utotm[i]  = utmax[i]
        yscal[i]  = 1000000.
    # i loop ends here

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

    for i in range(itlp-1,itu):
        xxa       = x[i,0,0]- x[i-1,0,0]
        yxa       = x[i,0,1]- x[i-1,0,1]
        volai     = 1.0/vola[i,0]
        uy        = 2.0* u[i,1]
        amub      = .5* (amu[i,0]+ amu[i,1])
        tauw[i]   = amub* (xxa* uy)* volai
        avor1     = vor[i,0]
        avora     = avor1
        avorb     = 0.5* (avor1+ avorm[i])
        avorc     = avorm[i]
        avor2     = tur1*avora + tur2*avorb + tur3*avorc
        yscal[i]  = np.sqrt(rey* sgrmi* amub* avor2* w[i,1,0])/(26.*amub)

        # **********************************************************************
        # *   compute normal distance ylen[i,j] and function  yvor             *
        # *   (yvor = y* vorticity)                                            *
        # **********************************************************************

    if (ncyc == ncyci1):
        for i in range(1,il):
            ylen[i,1] = 0.0
            for j in range(1,jl):
                xc2       = .50* (x[i,j,0]+ x[i-1,j,0]-x[i,j-1,0]- x[i-1,j-1,0])
                yc2       = .50* (x[i,j,1]+ x[i-1,j,1]-x[i,j-1,1]- x[i-1,j-1,1])
                xyc       = xc2
                yyc       = yc2
                scalf[j]  = np.sqrt(xyc*xyc + yyc*yyc)
                ylen[i,j] = ylen[i,j-1]+ scalf[j]

    for i in range(1,il):
        ylen1     = 0.5* ylen[i,1]
        for j in range(0,round(jstop)):
            y1        = yscal[i]* ylen[i,j]
            damp      = 1.0- np.exp(-y1)
            yvor[j]   = ylen[i,j]* vor[i,j]* damp
        # end j loop

        # i loop continues
        jmaxyv    = np.argmax(yvor)
        jmaxyv    = max(jmaxyv,2)
        jedge[i]  = jmaxyv

        # next line of code replaced because it caused convergence
        # stall when m = 0.001 - 12/10/05 (check this further !!!!)

        yvorm[i]  = max(yvor[jmaxyv],1.e-6)
        ylenm[i]  = max(ylen[i,jmaxyv],ylen1)

        if (jedge[i] < round(jstop)):
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

        # end i loop

        # **********************************************************************
        # *   compute outer eddy viscosity                                     *
        # *                                                                    *
        # *   outer do loop                                                    *
        # **********************************************************************

    for i in range(1,il): #start of outer i loop #####################
        udiff     = abs(utmax[i]- utmin[i])
        udiff1    = cwk1* udiff
        for j in range(1,round(jstop)): # loop 60
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
        # end loop 60
        amuto[1]  = amuto[2]

        # **********************************************************************
        # *   compute inner eddy viscosity                                     *
        # **********************************************************************

        for j in range(1,round(jstop)): # loop 70
            y1        = yscal[i]* ylen[i,j]
            damp      = 1.0- np.exp(-y1)
            tscali    = 0.4* ylen[i,j]* damp
            amuti1    = tscali* tscali* vor[i,j]
            amuti[j]   = rey* sgrmi* ravg[j]* amuti1
            amuti[j]   = abs(amuti[j])
        # end of loop 70
        amuti[1]  = 0.0
        if (i<=itl or i>itu):
            amuti[1] = amuti[2]

        # load viscosity coeffs. into array, use inner value until
        # match point is reached
        # scalar coding

        ivect     = 1
        if (ivect == 0): # start of big if statement
            icross    = 0
            amut[i,0] = amuti[0]
            for j in range(1,round(jstop)): # loop 75
                if (amuti[j]<=amuto[j] and icross==0): # nested if
                    amut[i,j] = amuti[j]
                else:
                    icross    = 1
                    amut[i,j] = amuto[j]
                # end nested if
            # end loop 75
        else: # else from the big if statement
            amut[i,0] = amuti[0]
            ystop     = round(jstop)
            for j in range(0,round(jstop)): # loop 80
                amudif    = amuti[j]- amuto[j]
                if (amudif >= 0): # if statement instead of cvmgp function
                    fcros[j] = j
                else: 
                    fcros[j] = 1000
            # end loop 80
            jcros    = np.argmin(fcros)
            if (jcros == 1):
                jcros = 2
            jcrosm   = jcros- 1
    
            for j in range(0,jcrosm): # loop 90
                amut[i,j] = amuti[j]
            # end loop 90
    
            for j in range(jcros-1,round(jstop)): # loop 100
                amut[i,j] = amuto[j]
            # end loop 100
        # end of big if statement

        # **********************************************************************
        # *   compute turbulent viscosity at cell center                       *
        # **********************************************************************

        for j in range(1,round(jstop)):
            amutc     = 0.5* (amut[i,j]+ amut[i,j-1])
            amu[i,j]  = amutc

        for j in range(1,round(jstop)):
            amut[i,j] = amu[i,j]

        if (i>itr1 and i<=itr2):
            for j in range(1,round(jstop)):
                amut[i,j] = 0.

# **********************************************************************
# *   debugging check (activate with jwrit = 1)                        *
# **********************************************************************
        # since jwrit is set to 0 anyways, just commenting this out for now
        #         jwrit     = 0
        #         if (jwrit == 1) then
        #         write(6,8000) i
        # 8000     format(5x,'i = ',i5)
        #         write(6,8050)
        # 8050     format(5x,'j,y1,fkleb,ravg,muti,muto,mut')
        # 8100     format(i5,6e15.5)
        #         do 195 j=2,jstop
        #             y1        = yscal[i]*ylen[i,j]
        #             cmuti     = amuti[j]
        #             cmuto     = amuto[j]
        #             write(6,8100) j,y1,fkleb[j],ravg[j],cmuti,cmuto,amut[i,j]
        # 195     continue
        #         end if

    # end of outer i loop #####################

    # copy amut to rev
    scale     = 1.
    for j in range(0,je):
        for i in range(0,ie):
            rev[i,j]  = scale*amut[i,j]

    print(np.mean(rev))
    return

    # 'return' statement here in turb2. So with current implementation, would end after this
    # therefore, commenting out rest for now

#     # **********************************************************************
#     # *   printout for boundary layer details                              *
#     # *   nturbw  -  input for frequency of printing b.l. data             *
#     # *              (currently set at large value)                        *
#     # **********************************************************************

#     nturbw = 5000

#     if (mod(ncyc,nturbw) == 0) then
#         do 250 i=2,il
#         if(i == 2) write(iwrit,840)
#         write (iwrit,850) i,jedge[i],ylenm[i],utotm[i],yvorm[i],
#     .                      yscal[i]
# 250   continue

#         do 300 i=itlp,itu
#         if(i==itlp) write(6,860)
#         je      = jedge[i]
#         delta1  = ylenm[i]
#         dstar1  = 1.0
#         uedge[i]= w(i,je,2)/w(i,je,1)
#         ue1     = uedge[i]* sgrmi
#         tauw1   = abs(tauw[i])
#         aylen   = abs(ylen[i,2])
#         yplus   = yscal[i]* aylen* 26.
#         yplus   = 0.5* yplus
#         cf1     = 2.0* rein* sgrmi* tauw1
#         write(6,870) i,je,ue1,delta1,dstar1,tauw1,cf1,yplus
# 300   continue

#         iloc    = itu
#         atauw   = abs(tauw(iloc))
#         atauwr  = atauw* rinv(iloc,2)
#         ustar   = sqrt(rein* rm* sgam* atauwr)
#         if(ustar==0.0) ustar = 1.0
#         ustari  = 1.0/ustar
# c
#         do 310 j=2,j2
#         yval    = 0.5* (ylen(iloc,j-1)+ ylen(iloc,j))
#         ayval   = abs(yval)
#         yplus   = yscal(iloc)* ayval* 26.
#         uval    = w(iloc,j,2)* rinv(iloc,j)
#         uplus   = uval* ustari
#         write(6,880) j,yplus,uplus,ylen(iloc,j),ylen(iloc,j-1),yval,
#     .                 yscal(iloc),yplus
# 310   continue
#     end if


#         # **********************************************************************
#         # *   additional checks of routine (activate with kwrite = 1)          *
#         # **********************************************************************

#     kwrite1   = 0
#     if (kwrite1==1 and ncyc==mcyc) then
#         do i=itl,itl-4,-4
#         do j=1,10
#         write(970,971) i,j,amut[i,j]
# 971     format (5x,'i,j,amut = ',2i5,e15.6)
#         end do
#         end do

#         itup      = itu + 1
#         do i=itup,itup+4,4
#         do j=1,10
#         write(970,971) i,j,amut[i,j]
#         end do
#         end do

# # ---->    finish of kwrite if
#         end if

#     return

# 800 format(1h ,'ylen,ravg,amuto1,amuto2,amuto,udiff,fkleb')
# 810 format(1h ,5x,i5,7e15.6)
# 820 format(1h ,'ylen,utot,vor,yvor,amuti,amuto,fkleb')
# 830 format(1h ,5x,i5,7e15.6)
# 840 format(1h ,'ylenm,utotm,yvorm,yscal')
# 850 format(1h ,5x,2i5,4e15.6)
# 860 format(1h ,1x,4h  i ,4h je ,1x,
#     .  10h    ue    ,10h   delta  ,10h   dstar  ,
#     .  10h   tauw   ,12h     cf     ,10h   yplus  )
# 870 format(1x,2i4,3f10.6,2e15.6,f10.6)
# 880 format(1h ,5x,i5,7f12.4)
#     end

dim_var = 100
params = {
  "ie": dim_var,
  "je": dim_var,
  "kvis": -1,
  "gamma": 1,
  "rm": 1,
  "re": 1,
  "ncyc": dim_var,
  "rev": np.ones((dim_var,dim_var)),
  "cmesh": 1,
  "ncyci1": -1,
  "itl": dim_var, 
  "itu": dim_var,
  "x": np.ones((dim_var,dim_var,3)),
  "w": np.ones((dim_var,dim_var,3)),
  "p": np.ones((dim_var,dim_var)),
  "vol": np.ones((dim_var,dim_var)),
  "xtran": 0,

  
}
dims = {
    "il": dim_var, 
    "jl": dim_var,
}

turbulent_viscosity(params, dims)
