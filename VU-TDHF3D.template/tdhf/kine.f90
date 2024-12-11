Subroutine kine(ncolx, ncoly, ncolz, lpsi, npsi, npmin, der2x, der2y, der2z, spkine, psi, pswk1, wxyz)
!-----------------------------------------------------------------------
!        kine:  calculates the expectation of laplacian for all
!               single particle states, with mass=1 and a '-' sign.
!               the result spkine is used to calculate the part of
!               the kinetic energy contribution to the total energy
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp), Intent(Out) :: spkine(lpsi, 2)
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, nst
      Complex(wp) :: cnormb
!-----------------------------------------------
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
!-----------------------------------------------------------------------
!        calculate and add second derivatives of psi
!-----------------------------------------------------------------------
            Call cmulx(ncolx, ncoly, ncolz, der2x, psi(1, 1, 1, 1, nst, iq), pswk1, 0)
            Call cmuly(ncolx, ncoly, ncolz, der2y, psi(1, 1, 1, 1, nst, iq), pswk1, 1)
            Call cmulz(ncolx, ncoly, ncolz, der2z, psi(1, 1, 1, 1, nst, iq), pswk1, 1)
!-----------------------------------------------------------------------
!        compute the expectation value
!-----------------------------------------------------------------------
            Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, iq), pswk1, cnormb, 1, wxyz)
            spkine(nst, iq) = - real(cnormb)
         End Do
      End Do
!
      Return
End Subroutine kine
