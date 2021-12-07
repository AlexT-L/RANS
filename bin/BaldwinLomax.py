import numpy as np
from Grid import Grid

class BaldwinLomax():
    def turbulent_viscosity(grid: Grid, params, dims):
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


      dimension       uedge(il),tauw(il),yscal(il),scalf(il),avor(il),
     .                avorm(il),ravg(il),amuto(il),amuti(il),yvor(il),
     .                yvorm(il),utotm(il),ylenm(il),fkleb(il),jedge(il)
      dimension       utmin(il),utmax(il),utot1(il),fcros(il)
      dimension       amu(ib,jb),u(ib,jb),v(ib,jb),vor(ib,jb),
     .                rinv(ib,jb),utot(ib,jb),vola(ib,jb),
     .                t(ib,jb),ylen(ib,jb)
      dimension       amut(ie,je)
      dimension       test(10)

      i2        = ie
      j2        = je
    il        = i2- 1
    jl        = j2- 1
      jlm       = jl- 1

      jstop     = 3* (j2- 2)/5
    if (cmesh .lt. 0.0) jstop = jl- 1
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
      sgam      = sqrt(gamma)
      sgrm      = sgam*rm
      sgrmi     = 1.0/sgrm

# c     **********************************************************************
# c     *    next 2 lines of code is activated to "test" turbulence model    *
# c     *      effect on convergence                                         *
# c     *                                                                    *
# c     *    turbulent viscosity is frozen after ncyct cycles                *
# c     **********************************************************************

     ncyct     = 10
     if (ncyc .gt. ncyct) return

      do 5 j=1,j2
      do 5 i=1,i2
        rinv(i,j) = 1.0/w(i,j,1)
        t(i,j)    = p(i,j)* rinv(i,j)
        u(i,j)    = w(i,j,2)* rinv(i,j)
        v(i,j)    = w(i,j,3)* rinv(i,j)
        amu(i,j)  = t(i,j)
        amut(i,j) = 0.0
    5 continue

# c     **********************************************************************
# c     *   determination of eddy viscosity                                  *
# c     *                                                                    *
# c     *   turbulence model:                                                *
# c     *     wall boundary layer --- baldwin-lomax model                    *
# c     *     wake region         --- baldwin-lomax model (cwake= 1.0)       *
# c     *                                                                    *
# c     *   calculate vorticity and total velocity                           *
# c     **********************************************************************

      do 10 j=1,jl
      do 10 i=1,i2
        vola(i,j) = 0.5* (vol(i,j)+ vol(i,j+1))
   10 continue

      do 15 i=2,il
        xxa       = x(i,1,1)-x(i-1,1,1)
        yxa       = x(i,1,2)-x(i-1,1,2)
        uy        = u(i,2)- u(i,1)
        vy        = v(i,2)- v(i,1)
        uavg      = .5* (u(i,1)+u(i,2))
        vavg      = .5* (v(i,1)+v(i,2))
        vor1      = (xxa*uy + yxa*vy)/vola(i,1)
        vor2      = 0.0
        vor3      = 0.0
        vort      = vor1-vor2-vor3
        vor(i,1)  = abs(vort)
        utotal    = uavg*uavg + vavg*vavg
        utot(i,1) = sqrt(utotal)
   15 continue

      do 20 j=2,jlm
      do 20 i=2,il
        xxa       = x(i,j,1)-x(i-1,j,1)
        yxa       = x(i,j,2)-x(i-1,j,2)
        uy        = u(i,j+1)- u(i,j)
        vy        = v(i,j+1)- v(i,j)
        uavg      = 0.5* (u(i,j)+ u(i,j+1))
        vavg      = 0.5* (v(i,j)+ v(i,j+1))

    #  thin-layer navier-stokes contribution to vorticity

        vor1      = (xxa*uy + yxa*vy)/vola(i,j)

    #  additional contributions to vorticity

        xyw       = 0.5* (x(i-1,j+1,1)- x(i-1,j-1,1))
        xye       = 0.5* (x(i,j+1,1)- x(i,j-1,1))
        yyw       = 0.5* (x(i-1,j+1,2)- x(i-1,j-1,2))
        yye       = 0.5* (x(i,j+1,2)- x(i,j-1,2))
        volawi    = 2.0/(vola(i,j)+ vola(i-1,j))
        volaei    = 2.0/(vola(i,j)+ vola(i+1,j))
        uxe       = 0.5* (u(i+1,j)+u(i+1,j+1)) - uavg
        vxe       = 0.5* (v(i+1,j)+v(i+1,j+1)) - vavg
        uxw       = uavg - 0.5* (u(i-1,j)+u(i-1,j+1))
        vxw       = vavg - 0.5* (v(i-1,j)+v(i-1,j+1))
        vor2      = 0.5* (xye* volaei* uxe+ xyw* volawi* uxw)
        vor3      = 0.5* (yye* volaei* vxe+ yyw* volawi* vxw)
        vort      = vor1- vor2- vor3
        vor(i,j)  = abs(vort)
        utotal    = uavg* uavg+ vavg* vavg
        utot(i,j) = sqrt(utotal)
   20 continue

    #  determine transition index

      itr1      = 0
      itr2      = 0
      j         = 1
      do i=1,il
        if (x(i,j,1) .le. xtran) then
          itr1      = i - 1
          go to 21
        end if
      end do
   21 continue

      itr1p     = itr1 + 1
      do i=itr1p,il
        if (x(i,j,1) .ge. xtran) then
          itr2      = i
          go to 22
        end if
      end do
   22 continue

      do 30 i=2,il

        do 25 j=1,jlm
          avor(j)   = vor(i,j)
          utot1(j)  = utot(i,j)
   25   continue

# c     effect of using jlm or jstop needs to be checked;
# c     it does not seem to make a difference
c
c       jmaxv     = ismax(jstop,avor,1)
        jmaxv     = ismax(jlm,avor,1)
        if (jmaxv .eq. 0) jmaxv = 1

        jminut    = ismin(jlm,utot1,1)

        jmaxut    = ismax(jstop,utot1,1)

        avorm(i)  = avor(jmaxv)
        utmin(i)  = utot1(jminut)
        utmax(i)  = max(utot1(jmaxut),1.e-3)
        utotm(i)  = utmax(i)
        yscal(i)  = 1000000.
   30 continue

      if (modbl .eq. 1) then
        tur1    = 1.0
        tur2    = 0.
        tur3    = 0.
      else if (modbl .eq. 2) then
        tur1    = 0.
        tur2    = 1.0
        tur3    = 0.
      else
        tur1    = 0.
        tur2    = 0.
        tur3    = 1.0
      end if

      do 35 i=itlp,itu
        xxa       = x(i,1,1)- x(i-1,1,1)
        yxa       = x(i,1,2)- x(i-1,1,2)
        volai     = 1.0/vola(i,1)
        uy        = 2.0* u(i,2)
        amub      = .5* (amu(i,1)+ amu(i,2))
        tauw(i)   = amub* (xxa* uy)* volai
        avor1     = vor(i,1)
        avora     = avor1
        avorb     = 0.5* (avor1+ avorm(i))
        avorc     = avorm(i)
        avor2     = tur1*avora + tur2*avorb + tur3*avorc
        yscal(i)  = sqrt(rey* sgrmi* amub* avor2* w(i,2,1))/(26.*amub)
   35 continue

        # **********************************************************************
        # *   compute normal distance ylen(i,j) and function  yvor             *
        # *   (yvor = y* vorticity)                                            *
        # **********************************************************************

    if (ncyc .eq. ncyci1) then
        do 39 i=2,il
          ylen(i,1) = 0.0
          do 37 j=2,jl
            xc2       = .50* (x(i,j,1)+ x(i-1,j,1)
     1                       -x(i,j-1,1)- x(i-1,j-1,1))
            yc2       = .50* (x(i,j,2)+ x(i-1,j,2)
     1                       -x(i,j-1,2)- x(i-1,j-1,2))
            xyc       = xc2
            yyc       = yc2
            scalf(j)  = sqrt(xyc*xyc + yyc*yyc)
            ylen(i,j) = ylen(i,j-1)+ scalf(j)
   37     continue
   39   continue
    end if

      do 50 i=2,il
        ylen1     = 0.5* ylen(i,2)
        do 40 j=1,jstop
          y1        = yscal(i)* ylen(i,j)
          damp      = 1.0- exp(-y1)
          yvor(j)   = ylen(i,j)* vor(i,j)* damp
   40   continue

        jmaxyv    = ismax(jstop,yvor,1)
        jmaxyv    = max(jmaxyv,2)
        jedge(i)  = jmaxyv

        # next line of code replaced because it caused convergence
        # stall when m = 0.001 - 12/10/05 (check this further !!!!)

        yvorm(i)  = max(yvor(jmaxyv),1.e-6)
        ylenm(i)  = max(ylen(i,jmaxyv),ylen1)

        if (jedge(i) .lt. jstop) then
          ylenm1  = ylenm(i)

          if (ncyc.ge.10 .or. restarr.eq.1.0) then
            jmyv    = jedge(i)
            dyvm    = yvor(jmyv)-yvor(jmyv-1)
            dyvp    = yvor(jmyv)-yvor(jmyv+1)

            if (yvor(jmyv-1) .lt. yvor(jmyv+1)) then
              ylenm(i) = ylen(i,jmyv)+ .5*(ylen(i,jmyv+1)- ylen(i,jmyv))
     .                   *(1.- dyvp/dyvm)
            else
              ylenm(i) = ylen(i,jmyv)- .5*(ylen(i,jmyv)- ylen(i,jmyv-1))
     .                   *(1.- dyvm/dyvp)
            end if

          else
            ylenm(i)  = ylenm1
          end if

        end if
   50 continue

        # **********************************************************************
        # *   compute outer eddy viscosity                                     *
        # *                                                                    *
        # *   outer do loop                                                    *
        # **********************************************************************

      do 200 i=2,il

        udiff     = abs(utmax(i)- utmin(i))
        udiff1    = cwk1* udiff
        do 60 j=2,jstop
          ravg(j)   = 0.5* (w(i,j,1)+ w(i,j+1,1))
          coeff     = 0.0168* ccp
          fwake1    = coeff* yvorm(i)* ylenm(i)
          coeff2    = coeff* cwk1* cwk1
          fwake2    = coeff2* ylenm(i)* udiff* udiff/yvorm(i)
          fwake     = min(fwake1,fwake2)
          fkleb0    = ckleb* ylen(i,j)/ylenm(i)
          fkleb1    = min(fkleb0,1.e5)
          fkleb(j)  = 1.0/(1.0+ 5.5* fkleb1**6)
          amuto(j)  = rey* sgrmi* ravg(j)* fwake* fkleb(j)
          amuto(j)  = abs(amuto(j))
   60   continue
        amuto(1)  = amuto(2)

        # **********************************************************************
        # *   compute inner eddy viscosity                                     *
        # **********************************************************************

        do 70 j=2,jstop
          y1        = yscal(i)* ylen(i,j)
          damp      = 1.0- exp(-y1)
          tscali    = 0.4* ylen(i,j)* damp
          amuti1    = tscali* tscali* vor(i,j)
          amuti(j)   = rey* sgrmi* ravg(j)* amuti1
          amuti(j)   = abs(amuti(j))
   70   continue
        amuti(1)  = 0.0
        if (i.le.itl .or. i.gt.itu) amuti(1) = amuti(2)

        # load viscosity coeffs. into array, use inner value until
        # match point is reached
        # scalar coding

        ivect     = 1
        if (ivect .eq. 0) then
          icross    = 0
          amut(i,1) = amuti(1)
          do 75 j=2,jstop
            if (amuti(j).le.amuto(j) .and. icross.eq.0) then
              amut(i,j) = amuti(j)
            else
              icross    = 1
              amut(i,j) = amuto(j)
            end if
   75     continue
        else
c
          amut(i,1) = amuti(1)
          ystop     = jstop
          do 80 j=1,jstop
            amudif    = amuti(j)- amuto(j)
            fcros(j)  = cvmgp(float(j),1000.,amudif)
   80     continue
          jcros    = ismin(jstop,fcros,1)
          if (jcros .eq. 1) jcros = 2
          jcrosm   = jcros- 1
c
          do 90 j=1,jcrosm
            amut(i,j) = amuti(j)
   90     continue
c
          do 100 j=jcros,jstop
            amut(i,j) = amuto(j)
  100     continue
        end if

        # **********************************************************************
        # *   compute turbulent viscosity at cell center                       *
        # **********************************************************************

        do 110 j=2,jstop
          amutc     = 0.5* (amut(i,j)+ amut(i,j-1))
          amu(i,j)  = amutc
  110   continue

        do 120 j=2,jstop
          amut(i,j) = amu(i,j)
  120   continue

        if (i.gt.itr1 .and. i.le.itr2) then
          do 130 j=2,jstop
            amut(i,j) = 0.
  130     continue
        end if

# **********************************************************************
# *   debugging check (activate with jwrit = 1)                        *
# **********************************************************************

        jwrit     = 0
        if (jwrit .eq. 1) then
          write(6,8000) i
 8000     format(5x,'i = ',i5)
          write(6,8050)
 8050     format(5x,'j,y1,fkleb,ravg,muti,muto,mut')
 8100     format(i5,6e15.5)
          do 195 j=2,jstop
            y1        = yscal(i)*ylen(i,j)
            cmuti     = amuti(j)
            cmuto     = amuto(j)
            write(6,8100) j,y1,fkleb(j),ravg(j),cmuti,cmuto,amut(i,j)
  195     continue
        end if

       # finish of outer loop
  200 continue

       # copy amut to rev

      scale     = 1.
      do j=1,je
      do i=1,ie
         rev(i,j)  = scale*amut(i,j)
      end do
      end do
      return

    # **********************************************************************
    # *   printout for boundary layer details                              *
    # *   nturbw  -  input for frequency of printing b.l. data             *
    # *              (currently set at large value)                        *
    # **********************************************************************

      nturbw = 5000

      if (mod(ncyc,nturbw) .eq. 0) then
        do 250 i=2,il
          if(i .eq. 2) write(iwrit,840)
          write (iwrit,850) i,jedge(i),ylenm(i),utotm(i),yvorm(i),
     .                      yscal(i)
  250   continue

        do 300 i=itlp,itu
          if(i.eq.itlp) write(6,860)
          je      = jedge(i)
          delta1  = ylenm(i)
          dstar1  = 1.0
          uedge(i)= w(i,je,2)/w(i,je,1)
          ue1     = uedge(i)* sgrmi
          tauw1   = abs(tauw(i))
          aylen   = abs(ylen(i,2))
          yplus   = yscal(i)* aylen* 26.
          yplus   = 0.5* yplus
          cf1     = 2.0* rein* sgrmi* tauw1
          write(6,870) i,je,ue1,delta1,dstar1,tauw1,cf1,yplus
  300   continue

        iloc    = itu
        atauw   = abs(tauw(iloc))
        atauwr  = atauw* rinv(iloc,2)
        ustar   = sqrt(rein* rm* sgam* atauwr)
        if(ustar.eq.0.0) ustar = 1.0
        ustari  = 1.0/ustar
c
        do 310 j=2,j2
          yval    = 0.5* (ylen(iloc,j-1)+ ylen(iloc,j))
          ayval   = abs(yval)
          yplus   = yscal(iloc)* ayval* 26.
          uval    = w(iloc,j,2)* rinv(iloc,j)
          uplus   = uval* ustari
          write(6,880) j,yplus,uplus,ylen(iloc,j),ylen(iloc,j-1),yval,
     .                 yscal(iloc),yplus
  310   continue
      end if


        # **********************************************************************
        # *   additional checks of routine (activate with kwrite = 1)          *
        # **********************************************************************

      kwrite1   = 0
      if (kwrite1.eq.1 .and. ncyc.eq.mcyc) then
        do i=itl,itl-4,-4
        do j=1,10
          write(970,971) i,j,amut(i,j)
  971     format (5x,'i,j,amut = ',2i5,e15.6)
        end do
        end do

        itup      = itu + 1
        do i=itup,itup+4,4
        do j=1,10
          write(970,971) i,j,amut(i,j)
        end do
        end do

# ---->    finish of kwrite if
        end if

      return

  800 format(1h ,'ylen,ravg,amuto1,amuto2,amuto,udiff,fkleb')
  810 format(1h ,5x,i5,7e15.6)
  820 format(1h ,'ylen,utot,vor,yvor,amuti,amuto,fkleb')
  830 format(1h ,5x,i5,7e15.6)
  840 format(1h ,'ylenm,utotm,yvorm,yscal')
  850 format(1h ,5x,2i5,4e15.6)
  860 format(1h ,1x,4h  i ,4h je ,1x,
     .  10h    ue    ,10h   delta  ,10h   dstar  ,
     .  10h   tauw   ,12h     cf     ,10h   yplus  )
  870 format(1x,2i4,3f10.6,2e15.6,f10.6)
  880 format(1h ,5x,i5,7f12.4)
      end
