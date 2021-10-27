      subroutine  bcresi
c
c*****   set boundary values to zero and implement cut bc's
c*****    for cell centered variables
c
      use dims
      use solv_var
c
      implicit none
      integer  :: i,j,ii
c
      do 10 i= 1,ie
        dw(i,1,1) = 0.0
        dw(i,1,2) = 0.0
        dw(i,1,3) = 0.0
        dw(i,1,4) = 0.0
c
        dw(i,je,1) = 0.0
        dw(i,je,2) = 0.0
        dw(i,je,3) = 0.0
        dw(i,je,4) = 0.0
   10 continue
c
      do 20 j= 1,je
        dw(1,j,1) = 0.0
        dw(1,j,2) = 0.0
        dw(1,j,3) = 0.0
        dw(1,j,4) = 0.0
c
        dw(ie,j,1) = 0.0
        dw(ie,j,2) = 0.0
        dw(ie,j,3) = 0.0
        dw(ie,j,4) = 0.0
   20 continue
c
      do 30 i= 1,ie
         if (i.le.itl .or. i.gt.itu) then
            ii    = ie-i+1
            dw(i,1,1) = dw(ii,2,1)
            dw(i,1,2) = dw(ii,2,2)
            dw(i,1,3) = dw(ii,2,3)
            dw(i,1,4) = dw(ii,2,4)
         end if
         if (i.gt.itl .or. i.lt.itl) then
            dw(i,1,1) = dw(i,2,1)
            dw(i,1,2) =-dw(i,2,2)
            dw(i,1,3) =-dw(i,2,3)
c           dw(i,1,2) = dw(i,2,2)
c           dw(i,1,3) = dw(i,2,3)
            dw(i,1,4) = dw(i,2,4)
         end if
   30 continue
c
      return
      end
