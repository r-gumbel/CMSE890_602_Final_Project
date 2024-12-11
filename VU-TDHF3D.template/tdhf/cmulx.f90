Subroutine cmulx(ncolx, ncoly, ncolz, xmat, pinn, pout, ifadd)
!-----------------------------------------------------------------------
!        multiply wavefunction with operator in x-space
!
!             ifadd = 0   -->  pout set to zero at the start
!             ifadd = 1   -->  result is added pout
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: ifadd
      Real(wp), Intent(In) :: xmat(ncolx, ncolx)
      Complex(wp), Intent(In) :: pinn(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(Inout) :: pout(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: is, ix, iy, iz, jx
!-----------------------------------------------
      If(ifadd == 0) Then
         pout = (0.0_wp, 0.0_wp)
      End If
!
      Do is = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  Do jx = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) + xmat(ix, jx) * pinn(jx, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End Do
!
      Return
End Subroutine cmulx
