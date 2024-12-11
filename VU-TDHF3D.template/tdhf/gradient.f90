Subroutine gradient(ncolx, ncoly, ncolz, d1x, d1y, d1z, finn, fout)
!
!-----------------------------------------------------------------------
!       GRADIENT
!       gradient of function finn returned in vector-function fout
!       d1x, d1y, d1z are the derivative operators in cartesian coord.
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: d1x(ncolx, ncolx)
      Real(wp), Intent(In) :: d1y(ncoly, ncoly)
      Real(wp), Intent(In) :: d1z(ncolz, ncolz)
      Real(wp), Intent(In) :: finn(ncolx, ncoly, ncolz)
      Real(wp), Intent(Inout) :: fout(ncolx, ncoly, ncolz, 3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ix, iy, iz
!-----------------------------------------------
      fout = 0.0_wp
!
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               fout(ix, iy, iz, 1) = fout(ix, iy, iz, 1) + Sum(d1x(ix, :)*finn(:, iy, iz))
            End Do
         End Do
      End Do
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               fout(ix, iy, iz, 2) = fout(ix, iy, iz, 2) + Sum(d1y(iy, :)*finn(ix, :, iz))
            End Do
         End Do
      End Do
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               fout(ix, iy, iz, 3) = fout(ix, iy, iz, 3) + Sum(d1z(iz, :)*finn(ix, iy, :))
            End Do
         End Do
      End Do
!
      Return
End Subroutine gradient
