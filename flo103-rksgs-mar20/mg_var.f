      module mg_var
c
c     ******************************************************************
c     *                                                                *
c     *   dynamically allocated variables for the multigrid scheme     *
c     *                                                                *
c     ******************************************************************
c
      real, dimension(:)        , allocatable :: ww,ww1,wwr
      real, dimension(:)        , allocatable :: xx,xxc,volc
      real, dimension(:)        , allocatable :: dtlc,rflc
      real, dimension(:)        , allocatable :: rlvc,revc
      real, dimension(:)        , allocatable :: radic,radjc
      real, dimension(:)        , allocatable :: qq,qbndc

      end module mg_var
