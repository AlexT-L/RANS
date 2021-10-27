      subroutine pscont(idim,jdim,is,ie,js,je,qt,qep)
c
c     this routine generates contour lines
c
c     arguments:
c
c     x(idim,*)       real  in ; x-coordinate of grid
c     y(idim,*)       real  in ; y-coordinate of grid
c     q(idim,*)       real  in ; physical scalar values on the grid
c     idim            int   in ; 1st dimension size of x, y and q
c     jdim            int   in ; 2nd dimension size of x, y and q
c     is, ie, js, je  int   in ; actual range
c     qt              real  in ; contour value
c     qep             real  in ; dimimum error allowance of q
c
c     local varibles:
c
c     qd(idim)        real  w  ; work array that shows qt-q
c     icf(idim,jdim)
c                     int   w  ; work array that shows contour crosses
c                                  that cell or not
c     xep(4,idim)     real  w  ; work array that keeps position that
c                                  edge and contour crosses -x.
c     yep(4,idim)     real  w  ; work array that keeps position that
c                                  edge and contour crosses -y.
c
c     ******************************************************************
c
      use out_var
c
c     implicit none
c
      integer is,ie,js,je
c     real  xq(idim,jdim),yq(idim,jdim),qc(idim,jdim)
      real  qt,qep
c
c     local variables
c
      real    qd(idim),xep(4,idim),yep(4,idim)
      real    qd1,qd2,qd3,qd4,frc,qave
c     integer icf(idim,jdim),imx1,jmx1,i,j,m1,m2
      integer imx1,jmx1,i,j,m1,m2
      integer inod(2,12)
      data    inod / 0, 0,  0, 0,  1, 2,  0, 0,
     &               1, 3,  2, 3,  0, 0,  0, 0,
     &               1, 4,  2, 4,  0, 0,  3, 4 /
c
c     ******************************************************************
c
      imx1 = ie - 1
      jmx1 = je - 1
c
c     set up
c
      do j=js, je
        do i=is, ie
          if (qc(i,j).eq.qt)  qc(i,j) = qc(i,j) + qep
        end do
      end do
c
      qd(is) = qt - qc(is,js)
      do 100 i=is,imx1
        qd4 = qt - qc(i+1,js)
        qd(i+1) = qd4
        if (qd(i)*qd4.lt.0.0)  then
c
c         crossed
c
          frc = qd4/(qc(i,js) - qc(i+1,js))
          xep(4,i) = (xq(i,js) - xq(i+1,js))*frc + xq(i+1,js)
          yep(4,i) = (yq(i,js) - yq(i+1,js))*frc + yq(i+1,js)
          icf(i,js) = 8
        else
          icf(i,js) = 0
        endif
  100 continue
c
c     main loop
c
      do 200 j=js,jmx1
        qd2 = qt - qc(is,j+1)
        if (qd(is)*qd2.lt.0.0)  then
c
c         crossed
c
          frc = qd2/(qc(is,j) - qc(is,j+1))
          xep(1,is) = (xq(is,j) - xq(is,j+1))*frc + xq(is,j+1)
          yep(1,is) = (yq(is,j) - yq(is,j+1))*frc + yq(is,j+1)
          icf(is,j) = icf(is,j) + 1
        endif
        do 210 i=is,imx1
          qd3 = qt - qc(i+1,j+1)
          if (qd3*qd2.lt.0.0)  then
c
c           crossed
c
            frc = qd3/(qc(i,j+1) - qc(i+1,j+1))
            xep(2,i) = (xq(i,j+1) - xq(i+1,j+1))*frc + xq(i+1,j+1)
            yep(2,i) = (yq(i,j+1) - yq(i+1,j+1))*frc + yq(i+1,j+1)
            icf(i,j) = icf(i,j) + 2
            icf(i,j+1) = 8
          else
            icf(i,j+1) = 0
          endif
c
c         for speed up
c
          if (icf(i,j).eq.0)  goto 230
c
          if (qd3*qd(i+1).lt.0.0)  then
c
c           crossed
c
            frc = qd3/(qc(i+1,j) - qc(i+1,j+1))
            xep(3,i) = (xq(i+1,j) - xq(i+1,j+1))*frc + xq(i+1,j+1)
            yep(3,i) = (yq(i+1,j) - yq(i+1,j+1))*frc + yq(i+1,j+1)
            icf(i,j) = icf(i,j) + 4
            icf(i+1,j) = icf(i+1,j) + 1
          endif
c
c         draw contour line
c
          if (icf(i,j).eq.15)  then
c
c           cross at 4 points
c
            qave = (qd(i) + qd(i+1) + qd2 + qd3)*0.25
            if (qave*qd3.gt.0.0)  then
c
c             draw 1-2 3-4
c
              call plot(xep(1,i), yep(1,i), 3)
              call plot(xep(2,i), yep(2,i), 2)
              call plot(xep(3,i), yep(3,i), 3)
              call plot(xep(4,i), yep(4,i), 2)
            else
c
c             draw 2-3 4-1
c
              call plot(xep(2,i), yep(2,i), 3)
              call plot(xep(3,i), yep(3,i), 2)
              call plot(xep(4,i), yep(4,i), 3)
              call plot(xep(1,i), yep(1,i), 2)
            endif
          elseif (icf(i,j).ne.0)  then
c
c           cross at 2 points
c
            m1 = inod(1,icf(i,j))
            m2 = inod(2,icf(i,j))
            call plot(xep(m1,i),yep(m1,i), 3)
            call plot(xep(m2,i),yep(m2,i), 2)
          endif
c
c         set-up for next loop
c
          xep(4,i) = xep(2,i)
          yep(4,i) = yep(2,i)
          xep(1,i+1) = xep(3,i)
          yep(1,i+1) = yep(3,i)
  230     continue
          qd(i) = qd2
          qd2 = qd3
  210     continue
c
c         end of i-loop
c
          qd(i) = qd3
  200   continue
c
        return
        end
