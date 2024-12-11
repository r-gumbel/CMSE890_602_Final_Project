Subroutine laplacian(ncolx, ncoly, ncolz, d2x, d2y, d2z, finn, fout)
!-----------------------------------------------------------------------
!       LAPLACIAN
!       laplacian of function finn returned in fout
!       d2x, d2y, d2z are second derivative operators in cartesian coord.
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: d2x(ncolx, ncolx)
      Real(wp), Intent(In) :: d2y(ncoly, ncoly)
      Real(wp), Intent(In) :: d2z(ncolz, ncolz)
      Real(wp), Intent(In) :: finn(ncolx, ncoly, ncolz)
      Real(wp), Intent(Inout) :: fout(ncolx, ncoly, ncolz)
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
               fout(ix, iy, iz) = fout(ix, iy, iz) + Sum(d2x(ix, :)*finn(:, iy, iz))
            End Do
         End Do
      End Do
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               fout(ix, iy, iz) = fout(ix, iy, iz) + Sum(d2y(iy, :)*finn(ix, :, iz))
            End Do
         End Do
      End Do
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               fout(ix, iy, iz) = fout(ix, iy, iz) + Sum(d2z(iz, :)*finn(ix, iy, :))
            End Do
         End Do
      End Do
!
      Return
End Subroutine laplacian
