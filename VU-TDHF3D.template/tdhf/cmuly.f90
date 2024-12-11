Subroutine cmuly(ncolx, ncoly, ncolz, ymat, pinn, pout, ifadd)
!-----------------------------------------------------------------------
!        multiply wavefunction with operator in y-space
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
      Real(wp), Intent(In) :: ymat(ncoly, ncoly)
      Complex(wp), Intent(In) :: pinn(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(Inout) :: pout(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: is, ix, iy, iz, jy
!-----------------------------------------------
      If(ifadd == 0) Then
         pout = (0.0_wp, 0.0_wp)
      End If
!
      Do is = 1, 2
         Do ix = 1, ncolx
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do jy = 1, ncoly
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) + ymat(iy, jy) * pinn(ix, jy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End Do
!
      Return
End Subroutine cmuly
