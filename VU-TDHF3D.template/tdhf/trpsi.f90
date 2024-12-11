Subroutine trpsi(ncolx, ncoly, ncolz, psi, tpsi)
!-----------------------------------------------------------------------
!     find the time reversal of the states
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(Out) :: tpsi(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ix, iy, iz
!-----------------------------------------------------------------------
!     array for the wf and it's time reversed state
!-----------------------------------------------------------------------
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               tpsi(ix, iy, iz, 1) = - conjg(psi(ix, iy, iz, 2))
               tpsi(ix, iy, iz, 2) = conjg(psi(ix, iy, iz, 1))
            End Do
         End Do
      End Do
!
      Return
End Subroutine trpsi
