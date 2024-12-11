      subroutine momentqx(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, qx, qx2, idp)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer  , intent(in)  :: ncolx
      integer  , intent(in)  :: ncoly
      integer  , intent(in)  :: ncolz
      integer  , intent(in)  :: idp
      real(wp) , intent(in)  :: xclx(ncolx)
      real(wp) , intent(in)  :: xcly(ncoly)
      real(wp) , intent(in)  :: xclz(ncolz)
      real(wp) , intent(in)  :: wx(ncolx)
      real(wp) , intent(in)  :: wy(ncoly)
      real(wp) , intent(in)  :: wz(ncolz)
      real(wp) , intent(in)  :: rho(ncolx,ncoly,ncolz,2)
      real(wp) , intent(out) :: xcm(2)
      real(wp) , intent(out) :: ycm(2)
      real(wp) , intent(out) :: zcm(2)
      real(wp) , intent(out) :: qx(7)
      real(wp) , intent(out) :: qx2(7)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(wp), parameter :: pi = acos(-1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(wp) :: pnrq, xcmq, ycmq, zcmq, z, y, x, xcmn, ycmn, zcmn,&
                  pnrtot, pnrn, pnrp, zcmp, ycmp, xcmp, xcmtot, ycmtot, zcmtot,&
                  x2cmq, y2cmq, z2cmq, x2cmn, y2cmn, z2cmn, x2cmp, y2cmp, z2cmp
      real(wp) :: vol
      integer  :: iq, ix, iy, iz
!-----------------------------------------------------------------------
!        calculate mass, center of mass, and charge coordinates
!-----------------------------------------------------------------------
      qx = 0.0_wp
      qx2 = 0.0_wp
!
      do iq = 1, 2
         pnrq = 0.0_wp
         xcmq = 0.0_wp
         ycmq = 0.0_wp
         zcmq = 0.0_wp
         x2cmq = 0.0_wp
         y2cmq = 0.0_wp
         z2cmq = 0.0_wp
         do iz = 1, ncolz
            z = xclz(iz)
            do iy = 1, ncoly
               y = xcly(iy)
               do ix = 1, ncolx
                  x = xclx(ix)
                  vol = wx(ix)*wy(iy)*wz(iz)*rho(ix,iy,iz,iq)
                  pnrq = vol + pnrq
                  xcmq = vol*x + xcmq
                  ycmq = vol*y + ycmq
                  zcmq = vol*z + zcmq
                  x2cmq = vol*x*x + x2cmq
                  y2cmq = vol*y*y + y2cmq
                  z2cmq = vol*z*z + z2cmq
               end do
            end do
         end do
         if(iq == 1) then
            pnrn = pnrq
            zcmn = zcmq
            ycmn = ycmq
            xcmn = xcmq
            x2cmn = x2cmq
            y2cmn = y2cmq
            z2cmn = z2cmq
         else
            pnrp = pnrq
            zcmp = zcmq
            ycmp = ycmq
            xcmp = xcmq
            x2cmp = x2cmq
            y2cmp = y2cmq
            z2cmp = z2cmq
         endif
      end do
      pnrtot = pnrn + pnrp
      xcmtot = (xcmn + xcmp)/pnrtot
      ycmtot = (ycmn + ycmp)/pnrtot
      zcmtot = (zcmn + zcmp)/pnrtot
      xcm(1) = xcmn/pnrn
      xcm(2) = xcmp/pnrp
      ycm(1) = ycmn/pnrn
      ycm(2) = ycmp/pnrp
      zcm(1) = zcmn/pnrn
      zcm(2) = zcmp/pnrp
!
      if(idp == 1) then
! dipole constraints (q10, Re(q11), Im(q11), set all to zero)
         qx(1) = sqrt(3.0_wp/(4.0_wp*pi))*xcmtot
         qx(2) = -sqrt(3.0_wp/(8.0_wp*pi))*ycmtot
         qx(3) = sqrt(3.0_wp/(8.0_wp*pi))*zcmtot
         qx2(1) = 3.0_wp/(4.0_wp*pi)*(x2cmn + x2cmp)/pnrtot
         qx2(2) = 3.0_wp/(8.0_wp*pi)*(y2cmn + y2cmp)/pnrtot
         qx2(3) = 3.0_wp/(8.0_wp*pi)*(z2cmn + z2cmp)/pnrtot
!
         do iq = 1, 2
            do iz = 1, ncolz
               z = xclz(iz)
               do iy = 1, ncoly
                  y = xcly(iy)
                  do ix = 1, ncolx
                     x = xclx(ix)
                     vol = wx(ix)*wy(iy)*wz(iz)*rho(ix,iy,iz,iq)
! Q21 real and imaginary parts (set to zero)
                     qx(4) = qx(4) - sqrt(15.0_wp/(8.0_wp*pi))*vol*x*y
                     qx(5) = qx(5) + sqrt(15.0_wp/(8.0_wp*pi))*vol*x*z
                     qx(6) = qx(6) + sqrt(15.0_wp/(8.0_wp*pi))*vol*y*z
                     qx(7) = qx(7) + sqrt(15.0_wp/(32.0_wp*pi))*vol*(y**2-z**2)
                     qx2(4) = qx2(4) + 15.0_wp/(8.0_wp*pi)*vol*(x*y)**2
                     qx2(5) = qx2(5) + 15.0_wp/(8.0_wp*pi)*vol*(x*z)**2
                     qx2(6) = qx2(6) + 15.0_wp/(8.0_wp*pi)*vol*(y*z)**2
                     qx2(7) = qx2(7) + sqrt(15.0_wp/(32.0_wp*pi))*vol*(y**2-z**2)**2
                  end do
               end do
            end do
         end do
      elseif(idp == 2) then
! dipole constraints (q10, Re(q11), Im(q11), set all to zero)
         qx(1) = sqrt(3.0_wp/(4.0_wp*pi))*ycmtot
         qx(2) = -sqrt(3.0_wp/(8.0_wp*pi))*zcmtot
         qx(3) = sqrt(3.0_wp/(8.0_wp*pi))*xcmtot
         qx2(1) = 3.0_wp/(4.0_wp*pi)*(y2cmn + y2cmp)/pnrtot
         qx2(2) = 3.0_wp/(8.0_wp*pi)*(z2cmn + z2cmp)/pnrtot
         qx2(3) = 3.0_wp/(8.0_wp*pi)*(x2cmn + x2cmp)/pnrtot
!
         do iq = 1, 2
            do iz = 1, ncolz
               z = xclz(iz)
               do iy = 1, ncoly
                  y = xcly(iy)
                  do ix = 1, ncolx
                     x = xclx(ix)
                     vol = wx(ix)*wy(iy)*wz(iz)*rho(ix,iy,iz,iq)
! Q21 real and imaginary parts (set to zero)
                     qx(4) = qx(4) - sqrt(15.0_wp/(8.0_wp*pi))*vol*y*z
                     qx(5) = qx(5) + sqrt(15.0_wp/(8.0_wp*pi))*vol*x*y
                     qx(6) = qx(6) + sqrt(15.0_wp/(8.0_wp*pi))*vol*x*z
                     qx(7) = qx(7) + sqrt(15.0_wp/(32.0_wp*pi))*vol*(z**2-x**2)
                     qx2(4) = qx2(4) + 15.0_wp/(8.0_wp*pi)*vol*(y*z)**2
                     qx2(5) = qx2(5) + 15.0_wp/(8.0_wp*pi)*vol*(x*y)**2
                     qx2(6) = qx2(6) + 15.0_wp/(8.0_wp*pi)*vol*(x*z)**2
                     qx2(7) = qx2(7) + sqrt(15.0_wp/(32.0_wp*pi))*vol*(z**2-x**2)**2
                  end do
               end do
            end do
         end do
      elseif(idp == 3) then
! dipole constraints (q10, Re(q11), Im(q11), set all to zero)
         qx(1) = sqrt(3.0_wp/(4.0_wp*pi))*zcmtot
         qx(2) = -sqrt(3.0_wp/(8.0_wp*pi))*xcmtot
         qx(3) = sqrt(3.0_wp/(8.0_wp*pi))*ycmtot
         qx2(1) = 3.0_wp/(4.0_wp*pi)*(z2cmn + z2cmp)/pnrtot
         qx2(2) = 3.0_wp/(8.0_wp*pi)*(x2cmn + x2cmp)/pnrtot
         qx2(3) = 3.0_wp/(8.0_wp*pi)*(y2cmn + y2cmp)/pnrtot
!
         do iq = 1, 2
            do iz = 1, ncolz
               z = xclz(iz)
               do iy = 1, ncoly
                  y = xcly(iy)
                  do ix = 1, ncolx
                     x = xclx(ix)
                     vol = wx(ix)*wy(iy)*wz(iz)*rho(ix,iy,iz,iq)
! Q21 real and imaginary parts (set to zero)
                     qx(4) = qx(4) - sqrt(15.0_wp/(8.0_wp*pi))*vol*z*x
                     qx(5) = qx(5) + sqrt(15.0_wp/(8.0_wp*pi))*vol*z*y
                     qx(6) = qx(6) + sqrt(15.0_wp/(8.0_wp*pi))*vol*x*y
                     qx(7) = qx(7) + sqrt(15.0_wp/(32.0_wp*pi))*vol*(x**2-y**2)
                     qx2(4) = qx2(4) + 15.0_wp/(8.0_wp*pi)*vol*(z*x)**2
                     qx2(5) = qx2(5) + 15.0_wp/(8.0_wp*pi)*vol*(z*y)**2
                     qx2(6) = qx2(6) + 15.0_wp/(8.0_wp*pi)*vol*(x*y)**2
                     qx2(7) = qx2(7) + sqrt(15.0_wp/(32.0_wp*pi))*vol*(x**2-y**2)**2
                  end do
               end do
            end do
         end do
      end if
!
      end subroutine momentqx
