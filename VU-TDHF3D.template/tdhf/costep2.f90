Subroutine costep2(iq, ncolx, ncoly, ncolz, rlam, psin, pswk)
!-----------------------------------------------------------------------
!     density-constraint multiplication
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: iq
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: rlam(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(In) :: psin(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(Inout) :: pswk(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: is, iz, iy, ix
!-----------------------------------------------
      Do is = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pswk(ix, iy, iz, is) = rlam(ix, iy, iz, iq) * psin(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
!
      Return
End Subroutine costep2
