      program flo103
c
c     ******************************************************************
c     *                                                                *
c     *   Analysis of profiles for transonic flow                      *
c     *                                                                *
c     ******************************************************************
c
c     ******************************************************************
c     *                                                                *
c     *   Copyright (c) Antony Jameson and Luigi Martinelli            *
c     *                                                                *
c     ******************************************************************
c
c     ******************************************************************
c     *                                                                *
c     *   References :                                                 *
c     *   (1) A. Jameson,"solution of the euler equations by a         *
c     *       multigrid method",Applied Mathematics and Computation,   *
c     *       vol. 13,1983,pp. 327-356                                 *
c     *   (2) A. Jameson,"transonic flow calculations for aircraft",   *
c     *       in numerical methods in fluid dynamics,edited by         *
c     *       F. Brezzi,Lecture Notes in Mathematics,vol. 1127,        *
c     *       Springer Verlag,1985,pp.156-242                          *
c     *   (3) A. Jameson,"multigrid algorithms for compressible flow   *
c     *       calculations",Proc. Second European Conference on        *
c     *       Multigrid Methods,Cologne,1985,edited by U. Trottenberg  *
c     *       and W. Hackbusch,Lecture Notes in Mathematics,vol. 1228, *
c     *       Springer Verlag,1986,pp. 166-201                         *
c     *   (4) A. Jameson,"computational transonics",                   *
c     *       Comm. Pure Appl. Math.,vol. 16,1988,pp. 507-549          *
c     *                                                                *
c     ******************************************************************
c
c     ******************************************************************
c     *                                                                *
c     *   The program has five main components                         *
c     *                                                                *
c     *               input and control module                         *
c     *               mesh generation module                           *
c     *               flow solver                                      *
c     *               multigrid module                                 *
c     *               output module                                    *
c     *                                                                *
c     ******************************************************************
c
c     Programmed by Antony Jameson and Luigi Martinelli
c     Flow solution method derived from flo82,
c     programmed by Antony Jameson,1986
c
c     ******************************************************************
c     *                                                                *
c     *   Comments in subroutine input describe the input parameters   *
c     *   and specify recommended values for general use               *
c     *                                                                *
c     ******************************************************************
c
c     ******************************************************************
c     *                                                                *
c     *   flow variables                                               *
c     *                                                                *
c     ******************************************************************
c
c     w(i,j,1)    = density
c     w(i,j,2)    = momentum in x direction
c     w(i,j,3)    = momentum in y direction
c     w(i,j,4)    = total energy
c
c     ******************************************************************
c     *                                                                *
c     *   multigrid control parameters                                 *
c     *                                                                *
c     ******************************************************************
c
c     nmesh         the number of meshes in the multigrid cycle
c     lcyc.eq.1     v cycle
c     lcyc.eq.2     w cycle
c     kxpand.eq.1   descend to a coarser mesh
c     kxpand.eq.-1  ascend to a finer mesh
c     mode.eq.0     the finest mesh in the cycle
c     mode.eq.1     any coarser mesh in the cycle
c     kode.eq.1     perform a normal update of the solution
c     kode.eq.2     first entry to a coarser mesh during the cycle
c     kode.eq.-1    calculate residuals without updating the solution
c     fcoll         relaxation factor for residual collection
c     fadd          smoothing factor for multigrid corrections
c     fbc.eq.1.     update the far field on the coarse meshes
c     fbc.eq.0.     freeze the far field on the coarse meshes
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use in_var
      use mesh_var
      use flo_var
      use psm_var
c
c     ******************************************************************
c
      use geo_param
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
c     local variables
c
c     ******************************************************************
c
      integer  :: iwrit
      integer  :: imesh
      integer  :: mlift,nlift,l
      integer  :: kxpand,kcyc(7)
c
c     ******************************************************************
c
c     real*4   :: tim,cutime
      real     :: tim,cutime
      real     :: work
c
c     ******************************************************************
c
      iwrit     = 6

      tim       = cutime(0)
      write (iwrit,700) tim

      write (iwrit,600)
      write (iwrit,2)
    2 format(1x,'Program flo103',50x,'Antony Jameson'/
     .       1x,'airfoil analysis in transonic flow',
     .          ' by solution of the unsteady RANS equations')
c
c     ******************************************************************
c     *                                                                *
c     *   license check                                                *
c     *                                                                *
c     ******************************************************************
c
c     call syslock
c
c     ******************************************************************
c     *                                                                *
c     *   set up the aerodynamic and geometric data                    *
c     *                                                                *
c     ******************************************************************
c
c     read the input data
c
      call input

      gcl       = fcl
      k1        = 1
      kx        = 1
      kw        = 1

      mode      = -1
c
c     increase the dimensions to those of the finest mesh
c     in the sequence
c
      if (mmesh.gt.1) then
         do lmesh=2,mmesh
            nx        = nx  +nx
            ny        = ny  +ny
         end do
      end if

      il        = nx  +1
      jl        = ny  +1
      ie        = nx  +2
      je        = ny  +2
      ib        = nx  +3
      jb        = ny  +3
c
c     allocate the arrays
c
      call alloc_flo
      call alloc_mesh
      call alloc_solv
      call alloc_psm
      call alloc_mg
      call alloc_geo
      call alloc_out
c
c     ******************************************************************
c     *                                                                *
c     *   generate the mesh for the final calculation                  *
c     *                                                                *
c     ******************************************************************
c
      call gmesh
c
c     expand the mesh to produce the first mesh in the sequence
c
      if (mmesh.gt.1) then
         do lmesh=2,mmesh
            call movco
         end do
      end if

      lmesh     = 1
c
c     set up the metric and initial data
c
   11 call setup
c
c     ******************************************************************
c     *                                                                *
c     *   calculate the flow                                           *
c     *                                                                *
c     ******************************************************************
c
      tim       = cutime(0)
      write (iwrit,700) tim
      mcyc      = fcyc(lmesh)
      mlift     = 0
      if (kvis.gt.0) mlift = 20
      nlift     = 1
      ntim      = ftim(lmesh)
      nout      = fout(lmesh)
      nprnt     = fprnt(lmesh)
      iprnt     = gprnt(lmesh)
      lprnt     = hprnt(lmesh)
      lprnt     = max(1,lprnt)
      nmesh     = hmesh(lmesh)

      do imesh=1,nmesh
         kcyc(imesh) = 0
      end do

      l         = (mcyc  -1)/500  +1
      nout      = max(nout,l)
      nres      = 0
      ncyc      = 0
      tot       = 0.
      imesh     = nmesh
      imon      = 0
c
c     print data to monitor convergence
c
      call monitr
      imon      = 1
c
c     ******************************************************************
c     *                                                                *
c     *   start of the multigrid cycle                                 *
c     *                                                                *
c     ******************************************************************
c
  101 ncyc      = ncyc  +1
      kode      = 1
      mode      = 0
      eps       = epsf
      cfl       = cflf
      cfl0      = cfl

      if (ncyc .gt. 5) fadd = 0.0
      if (iprec.ne.0 .and. ncyc.le.5) cfl = min(cfl0,10.)
c     if (iprec.ne.0 .and. ncyc.eq.1) cfl = min(cfl0,10.)

      hm        = hmf
      work      = 1.
      nsup      = 0

  111 kxpand    = 1
      kcyc(imesh) = kcyc(imesh)  +1

  121 tot       = tot  +work
c
c     calculate the viscosity coefficients
c
      if (kvis.gt.0.and.mode.eq.0) call viscf
c     if (kvis.gt.0.and.mode.eq.0.and.ncyc.le.25) call viscf
c     if (kvis.gt.0.and.mode.eq.0.and.ncyc.le.50) call viscf
c     if (kvis.gt.0.and.mode.eq.0.and.ncyc.le.100) call viscf
c
c     calculate the permissible time step
c
      if (mod(ncyc-1,ntim).eq.0) call step
c
c     update the flow
c
      fcl       = 0.
      if (ncyc.gt.mlift.and.mod(ncyc,nlift).eq.0.and.ncyc.lt.mcyc)
     .fcl       = gcl
      call euler
      kode      = 1
c
c     print data to monitor convergence
c
      if (imesh.eq.nmesh.and.mod(ncyc-1,nout).eq.0)
     .call monitr

      if (ncyc.ge.mcyc) go to 201

      if (mod(ncyc,nprnt).eq.0.and.imesh.eq.nmesh) call output

      if (kxpand.eq.-1.and.imesh.lt.nmesh) go to 141
      if (imesh.eq.1) go to 101

      kode      = -1
c
c     calculate the residuals
c
      call euler
c
c     transfer to the next coarser mesh
c
      call movco
      kode      = 2
      mode      = 1
      imesh     = imesh  -1
      eps       = epsc
      cfl       = cflc
      cfl0      = cfl

      if (iprec.ne.0 .and. ncyc.le.5) cfl = min(cfl0,10.)
c     if (iprec.ne.0 .and. ncyc.eq.1) cfl = min(cfl0,10.)

      hm        = hmc
      work      = .25*work

      if (imesh.gt.1) go to 111

      kxpand    = -1
      go to 121
c
c     transfer to the next finer mesh
c
  141 call movfin
      imesh     = imesh  +1
      work      = 4.*work

      if (imesh.eq.nmesh) go to 101

      if (kcyc(imesh).lt.lcyc) go to 111
      kcyc(imesh) = 0
      go to 141
c
c     ******************************************************************
c     *                                                                *
c     *   end of the multigrid cycle                                   *
c     *                                                                *
c     ******************************************************************
c
c     print data to monitor convergence
c
  201 fcl       = gcl
      imon      = 2
      call monitr
c
c     ******************************************************************
c     *                                                                *
c     *   postprocess the result                                       *
c     *                                                                *
c     ******************************************************************
c
      call output
      tim       = cutime(0)
      write (iwrit,700) tim
c
c     ******************************************************************
c     *                                                                *
c     *   transfer to a finer mesh if m.lt.mmesh                       *
c     *                                                                *
c     ******************************************************************
c
      if (lmesh.ge.mmesh) stop
      lmesh     = lmesh  +1
      call movfin
      go to 11
c
c     ******************************************************************
c
  600 format(1x)
  700 format(1x,'time ',f10.2,'seconds')

      end
