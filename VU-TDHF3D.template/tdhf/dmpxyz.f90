Subroutine dmpxyz(ncolx, ncoly, ncolz, psinn, psout, cdmpx, cdmpy, cdmpz)
!-----------------------------------------------------------------------
!        dampening with the inverse kinetic-energy
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: cdmpx(ncolx, ncolx)
      Real(wp), Intent(In) :: cdmpy(ncoly, ncoly)
      Real(wp), Intent(In) :: cdmpz(ncolz, ncolz)
      Complex(wp) :: psinn(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: psout(ncolx, ncoly, ncolz, 2)
!
!        x-product
!
      Call cmulx(ncolx, ncoly, ncolz, cdmpx, psinn, psout, 0)
!
!        y-product
!
      Call cmuly(ncolx, ncoly, ncolz, cdmpy, psout, psinn, 0)
!
!        z-product
!
      Call cmulz(ncolx, ncoly, ncolz, cdmpz, psinn, psout, 0)
!
      Return
End Subroutine dmpxyz
