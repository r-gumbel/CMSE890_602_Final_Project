Subroutine tstep(ncolx, ncoly, ncolz, iq, mxp, terr, dt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, der1x, &
& der1y, der1z, der2x, der2y, der2z, psout, pswk1, pswk2, pswk3, wxyz, itimrev)
!-----------------------------------------------------------------------
!     tstep   =  one time step by exponential expansion of the
!                propagator.
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: iq
      Integer, Intent(In) :: mxp
      Real(wp), Intent(In) :: terr
      Real(wp), Intent(In) :: dt
      Real(wp), Intent(In) :: upot(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: bmass(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: xiq(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: hsigma(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: cq(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: dcq(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp), Intent(In) :: bmunu(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp), Intent(In) :: dbmu(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: dbmass(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: dxiq(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Complex(wp) :: psout(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk2(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk3(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: hbc = 197.3269885804381_wp
      Real(wp), Parameter :: esf = 0.0_wp
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: m
      Real(wp) :: error, cts1, denom
      Complex(wp) :: fmd, cts0
!
      cts1 = 2.0_wp
      error = 1.0_wp
      m = 0
!-----------------------------------------------------------------------
!        compute exp(-eye*dt*h) by power series expansion
!-----------------------------------------------------------------------
      pswk1 = psout
      Do while(error >= terr .And. m <= mxp)
         m = m + 1
         denom = m * hbc
         fmd = - eye * dt / denom
!
         Call hpsi(ncolx, ncoly, ncolz, upot(1, 1, 1, iq), bmass(1, 1, 1, iq), xiq(1, 1, 1, 1, iq), hsigma(1, 1, 1, 1, iq), cq(1, &
        & 1, 1, 1, iq), dcq(1, 1, 1, 1, 1, iq), bmunu(1, 1, 1, 1, 1, iq), dbmu(1, 1, 1, 1, iq), dbmass(1, 1, 1, 1, iq), dxiq(1, &
        & 1, 1, iq), der1x, der1y, der1z, der2x, der2y, der2z, esf, pswk1, pswk2, pswk3, itimrev)
!
         pswk1 = fmd * pswk2
         psout = psout + pswk1
!
         Call psnorm(ncolx, ncoly, ncolz, psout, psout, cts0, 1, wxyz)
!
         error = Abs(Abs(cts0)-cts1)
         cts1 = Abs(cts0)
      End Do
!      write (*, '(A,2I3,1P,2E12.4)') ' tstep: m,mxp,error,terr =', m, mxp, error, terr
!
End Subroutine tstep
