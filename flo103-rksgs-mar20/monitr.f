      subroutine monitr
c
c     ******************************************************************
c     *                                                                *
c     *   monitor the convergence                                      *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use mesh_var
      use out_var
c
c     ******************************************************************
c
      use flo_param
      use solv_param
      use mg_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      common/res/ res(1001),fsup(1001),count(1001)
c
c     ******************************************************************
c
      common/tit/ title
c
c     ******************************************************************
c
      character(80) :: title
c
c     ******************************************************************
c
      real     :: res,fsup,count
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: iwrit
c
c     ******************************************************************
c
      real     :: rad,al,dalpha
      real     :: rate1,rate2
      real     :: fvis,fsmoop
      real     :: prec
c
c     ******************************************************************
c
      iwrit     = 6

      rad       = 45./atan(1.)
      al        = alpha*rad

      if (imon.gt.1) go to 41
      if (imon.gt.0) go to 31
c
c     option to print the mesh
c
      if (iprnt.gt.1) call prntx (lprnt)
c
c     ******************************************************************
c     *                                                                *
c     *   data for the flow solution                                   *
c     *                                                                *
c     ******************************************************************
c
      write (iwrit,600)
      write (iwrit,12)
   12 format(1x,'Time dependent solution')
      write (iwrit,14)
   14 format(3x,' mach number',' ang attack ',' reynolds no',
     .          '    fvis    ','     fcl    ','     clt    ')
      fvis      = kvis
      write (iwrit,611) rm,al,re,fvis,fcl,clt
      write (iwrit,16)
   16 format(3x,'     nx     ','     ny     ','   n mesh   ',
     .          '    n cyc   ','   n stage  ')
      write (iwrit,640) nx,ny,nmesh,lcyc,mstage
      write (iwrit,18)
   18 format(3x,' cfl number ','   fsmoop   ','   smoopi   ',
     .          '   smoopj   ','  h factor  ','     vt     ')
      fsmoop    = ksmoop
      write (iwrit,610) cflf,fsmoop,smoopi,smoopj,hmf,vt
      write (iwrit,20)
   20 format(3x,'   cstp(1)  ','   cstp(2)  ','   cstp(3)  ',
     .          '   cstp(4)  ','   cstp(5)  ','   cstp(6)  ')
      write (iwrit,610) cstp(1),cstp(2),cstp(3),cstp(4),cstp(5),cstp(6)
      write (iwrit,22)
   22 format(3x,'   cdis(1)  ','   cdis(2)  ','   cdis(3)  ',
     .          '   cdis(4)  ','   cdis(5)  ','   cdis(6)  ')
      write (iwrit,610) cdis(1),cdis(2),cdis(3),cdis(4),cdis(5),cdis(6)
      write (iwrit,24)
   24 format(3x,'    vis 2   ','    vis 4   ','    vis 0   ',
     .          '    adis    ','     bc     ')
      write (iwrit,610) vis2,vis4,vis0,adis,bc
      write (iwrit,26)
   26 format(3x,'    cflc    ','     fbc    ','     hmc    ')
      write (iwrit,610) cflc,fbc,hmc
      prec = iprec
      write (iwrit,28)
   28 format(3x,'    fcoll   ','    fadd    ','    iprec    ',
     .          '    epsf    ','    epsc    ','    diag     ')
      write (iwrit,610) fcoll,fadd,prec,epsf,epsc,diag
      write (iwrit,30)
   30 format(3x,'Cycle',
     .          ' Max dr/dt','  i ','  j ',' Avg dr/dt',
     .          'Max enthpy','  i ','  j ','Avg enthpy',
     .          '   Alpha  ','    CL    ','    CD    ',
     .          '  N sup ')

      return
c
c     ******************************************************************
c     *                                                                *
c     *   monitor the flow solution                                    *
c     *                                                                *
c     ******************************************************************
c
   31 write (iwrit,650) ncyc,rtmax,irt,jrt,rtrms,
     .                  hmax,ih,jh,hrms,al,cl,cd,nsup
      call flush (iwrit)
      if (ncyc.eq.1) rtrms0 = rtrms
      nres      = nres  +1
c     count(nres) = tot  -1.
      count(nres) = ncyc  -1
      res(nres)   = rtrms
      fsup(nres)  = nsup
c
c     option to calculate the flow at a given lift coefficient
c
      if (fcl.le.0.) return

      if (kvis.eq.0) then
c        dalpha    = .0625*fcl*(clt  -cl)
c        dalpha    = .01*fcl*(clt  -cl)
c        dalpha    = .02*fcl*(clt  -cl)
         dalpha    = .025*fcl*(clt  -cl)
c        dalpha    = .03*fcl*(clt  -cl)
c        dalpha    = .04*fcl*(clt  -cl)
      else
c        dalpha    = .03125*fcl*(clt  -cl)
c        dalpha    = .002*fcl*(clt  -cl)
c        dalpha    = .003*fcl*(clt  -cl)
c        dalpha    = .004*fcl*(clt  -cl)
         dalpha    = .005*fcl*(clt  -cl)
      end if
      alpha     = alpha  +dalpha
c     alpha     = min(alpmax,alpha)
      ca        = cos(alpha)
      sa        = sin(alpha)
      u0        = rm*c0*ca
      v0        = rm*c0*sa

      return
c
c     ******************************************************************
c     *                                                                *
c     *   convergence of the flow solution                             *
c     *                                                                *
c     ******************************************************************
c
c     calculate the average contraction ratio of the iterations
c
   41 rate1     = (rtrms/rtrms0)**(1./(tot  -1.))
      rate2     = (rtrms/rtrms0)**(1./float(ncyc  -1))
      write (iwrit,42)
   42 format(3x,'   Rtrms 1  ','   Rtrms 2  ','    Work    ',
     .          'Reductn/work',' Reductn/cyc')
      write (iwrit,660) rtrms0,rtrms,tot,rate1,rate2
c
c     plot the convergence history
c
      call rplot (nres,res,fsup,count,title,rm,al,nx,ny)

      return
c
c     ******************************************************************
c
  600 format(1x)
  610 format(1x,6f12.4)
  611 format(1x,2f12.4,f12.2,3f12.4)
  640 format(1x,6i12)
  650 format(1x,i5,e10.3,2i4,2e10.3,2i4,e10.3,3f10.5,i7)
  660 format(1x,2e12.4,3f12.4)

      end
