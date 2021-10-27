      module flo_var
c
c     ******************************************************************
c     *                                                                *
c     *   dynamically allocated variables for the flow solution        *
c     *                                                                *
c     ******************************************************************
c
      real, dimension(:,:,:)    , allocatable :: w
      real, dimension(:,:)      , allocatable :: p

      real, dimension(:,:)      , allocatable :: rlv,rev,aev,bev
      real, dimension(:)        , allocatable :: tw,cfric
      real, dimension(:)        , allocatable :: dsti,ynot,ssmax

      end module flo_var
