      module psm_var
c
c     ******************************************************************
c     *                                                                *
c     *   dynamically allocated variables for the solution             *
c     *                                                                *
c     ******************************************************************
c
      real, dimension(:,:,:)    , allocatable :: rs
      real, dimension(:,:)      , allocatable :: r00,p00,u00,v00

      end module psm_var
