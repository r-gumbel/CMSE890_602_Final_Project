Subroutine inertia(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcmtot, ycmtot, zcmtot, rho, w, z)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: xcmtot
      Real(wp), Intent(In) :: ycmtot
      Real(wp), Intent(In) :: zcmtot
      Real(wp), Intent(Out) :: w(3)
      Real(wp), Intent(Out) :: z(3, 3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: xi(3, 3), xx, yy, zz, vol
      Integer  :: ix, iy, iz
!      Integer  :: i
!-------------------------------------------------------
!   calculate the traceless-symmetric quadrupole tensor
!-------------------------------------------------------
      xi = 0.0_wp
      w  = 0.0_wp
      z  = 0.0_wp
      Do iz = 1, ncolz
         zz = xclz(iz) - zcmtot
         Do iy = 1, ncoly
            yy = xcly(iy) - ycmtot
            Do ix = 1, ncolx
               xx = xclx(ix) - xcmtot
               vol = wx(ix) * wy(iy) * wz(iz) * (rho(ix, iy, iz, 1)+rho(ix, iy, iz, 2))
               xi(1, 1) = xi(1, 1) + (yy*yy+zz*zz) * vol
               xi(2, 2) = xi(2, 2) + (xx*xx+zz*zz) * vol
               xi(3, 3) = xi(3, 3) + (xx*xx+yy*yy) * vol
               xi(1, 2) = xi(1, 2) - xx * yy * vol
               xi(1, 3) = xi(1, 3) - xx * zz * vol
               xi(2, 3) = xi(2, 3) - yy * zz * vol
            End Do
         End Do
      End Do
      xi(2, 1) = xi(1, 2)
      xi(3, 1) = xi(1, 3)
      xi(3, 2) = xi(2, 3)
!
      Call jacobi(xi, 3, 3, w, z)
! we rotate about y-axis
!      WRITE(8,*) 'Cartesian inertia principal values and axes:'
!      WRITE(8,'(/,3(F10.2,''('',3F8.4,'')''))') (w(i),z(:,i),i=1,3)
!
End Subroutine inertia
