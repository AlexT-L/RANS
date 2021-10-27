      module mesh_var
c
c     ******************************************************************
c     *                                                                *
c     *   dynamically allocated variables for the mesh                 *
c     *                                                                *
c     ******************************************************************
c
      real, dimension(:,:,:)    , allocatable :: x,xc
      real, dimension(:,:)      , allocatable :: vol
      real, dimension(:,:)      , allocatable :: pori,porj
      real, dimension(:,:)      , allocatable :: fint

      end module mesh_var
