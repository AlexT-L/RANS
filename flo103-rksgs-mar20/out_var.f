      module out_var
c
c     ******************************************************************
c     *                                                                *
c     *   dynamically allocated variables for the output               *
c     *                                                                *
c     ******************************************************************
c
      real, dimension(:)      , allocatable :: xp,yp,cp,cpt
      real, dimension(:)      , allocatable :: cf
      real, dimension(:)      , allocatable :: rcl,rcd,rld

      real, dimension(:,:)    , allocatable :: xp0,yp0,cp0

      real, dimension(:,:)    , allocatable :: xq,yq,qc

      integer, dimension(:,:) , allocatable :: icf

      end module out_var
