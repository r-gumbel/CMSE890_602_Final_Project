Subroutine grstep(ncolx, ncoly, ncolz, iq, mprint, iprint, itrs, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
& der1x, der1y, der1z, der2x, der2y, der2z, spe, efluct, denerg, psin, pswk1, pswk2, pswk3, cdmpx, cdmpy, cdmpz, x0dmp, itrbx, &
& irest, wxyz, xcm, ycm, zcm, xclx, xcly, xclz, xlam, xlamx, damp, iconstr, idp, itimrev, cnormb)
!-----------------------------------------------------------------------
!     grstep   =  one damped gradient iteration step for a given
!                 wave function of state nst with isospin iq.
!        psi   =  o[ psi - x0*cdmpxyz*[(h-spe)psi] ]
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: iq
      Integer, Intent(In) :: mprint
      Integer, Intent(In) :: iprint
      Integer, Intent(In) :: itrs
      Integer, Intent(In) :: itrbx
      Integer, Intent(In) :: irest
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: iconstr
      Integer, Intent(In) :: idp
      Real(wp), Intent(In) :: x0dmp
      Real(wp), Intent(In) :: xlam
      Real(wp), Intent(In) :: xlamx(5)
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
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
      Real(wp), Intent(In) :: cdmpx(ncolx, ncolx)
      Real(wp), Intent(In) :: cdmpy(ncoly, ncoly)
      Real(wp), Intent(In) :: cdmpz(ncolz, ncolz)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: damp(ncolx, ncoly, ncolz)
      Real(wp), Intent(Out) :: efluct
      Real(wp), Intent(Out) :: denerg
      Real(wp), Intent(Inout) :: spe
      Complex(wp) :: cnormb
      Complex(wp) :: psin(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk2(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk3(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: esf, ar, enrold, temarg
      Complex(wp) :: cexph2
      cnormb = (0.0_wp, 0.0_wp)
!-----------------------------------------------------------------------
!        action of (h-esf) on psi yields pswk1
!        action of (h-esf)**2 on psi yields pswk2.
!        esf is an energy shift which is later subtracted from spe
!        note: pswk2 is calculated only for writing output
!-----------------------------------------------------------------------
      esf = spe
      Call hpsi(ncolx, ncoly, ncolz, upot(1, 1, 1, iq), bmass(1, 1, 1, iq), xiq(1, 1, 1, 1, iq), hsigma(1, 1, 1, 1, iq), cq(1, 1, &
     & 1, 1, iq), dcq(1, 1, 1, 1, 1, iq), bmunu(1, 1, 1, 1, 1, iq), dbmu(1, 1, 1, 1, iq), dbmass(1, 1, 1, 1, iq), dxiq(1, 1, 1, &
     & iq), der1x, der1y, der1z, der2x, der2y, der2z, esf, psin, pswk1, pswk3, itimrev)
!-----------------------------------------------------------------------
!        calculate the expectation cnormb=<h-esf>
!-----------------------------------------------------------------------
      Call psnorm(ncolx, ncoly, ncolz, psin, pswk1, cnormb, 1, wxyz)
!-----------------------------------------------------------------------
!        for printing calculate fluctuation, i.e. <h*h>
!-----------------------------------------------------------------------
      If(mprint /= 0) Then
         If(iprint > 0 .And. (Mod(itrs, mprint) == 0 .Or. itrs == irest+itrbx)) Then
            Call hpsi(ncolx, ncoly, ncolz, upot(1, 1, 1, iq), bmass(1, 1, 1, iq), xiq(1, 1, 1, 1, iq), hsigma(1, 1, 1, 1, iq), &
           & cq(1, 1, 1, 1, iq), dcq(1, 1, 1, 1, 1, iq), bmunu(1, 1, 1, 1, 1, iq), dbmu(1, 1, 1, 1, iq), dbmass(1, 1, 1, 1, iq), &
           & dxiq(1, 1, 1, iq), der1x, der1y, der1z, der2x, der2y, der2z, esf, pswk1, pswk2, pswk3, itimrev)
!
            Call psnorm(ncolx, ncoly, ncolz, psin, pswk2, cexph2, 1, wxyz)
         End If
      End If
!-----------------------------------------------------------------------
!        repurpose psiwk3 for hmatrix generation
!-----------------------------------------------------------------------
      pswk3 = pswk1
!-----------------------------------------------------------------------
!        operate with the damping matrix and subtract from psi
!-----------------------------------------------------------------------
      ar = real(cnormb)
      pswk1 = pswk1 - ar * psin
! dampen the quadrupole constraint step as well
      If(iconstr == 1) Then ! just the quadrupole constraint
         Call costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm, xlam, damp, psin, pswk1, idp)
      ElseIf(iconstr == 3) Then ! quadrupole constraint + on <x>,<y>,<z>, and q21
         Call costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm, xlam, damp, psin, pswk1, idp)
         Call costep3(ncolx, ncoly, ncolz, xclx, xcly, xclz, xlamx, damp, psin, pswk1, idp)
      ElseIf(iconstr == 4) Then ! Everyhing i iconstr=3 + q22 axial constraint
         Call costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm, xlam, damp, psin, pswk1, idp)
         Call costep3(ncolx, ncoly, ncolz, xclx, xcly, xclz, xlamx, damp, psin, pswk1, idp)
         Call costep4(ncolx, ncoly, ncolz, xclx, xcly, xclz, xlamx, damp, psin, pswk1, idp)
      End If
!
      Call dmpxyz(ncolx, ncoly, ncolz, pswk1, pswk2, cdmpx, cdmpy, cdmpz)
!
      pswk1 = psin - x0dmp * pswk2
!-----------------------------------------------------------------------
!        copy function back into psi and calculate new spe and de/e
!-----------------------------------------------------------------------
      psin = pswk1
!-----------------------------------------------------------------------
!        calculate convergence and optionally the fluctuations
!-----------------------------------------------------------------------
      enrold = spe
      spe = cnormb + esf
      denerg = (enrold-spe) / Abs(spe+1.0e-25_wp)
!
      If(mprint /= 0) Then
         If(iprint > 0 .And. (Mod(itrs, mprint) == 0 .Or. itrs == irest+itrbx)) Then
            temarg = real(cexph2) - real(cnormb) ** 2
            efluct = Sqrt(Abs(temarg))
         End If
      End If
!
      Return
End Subroutine grstep
