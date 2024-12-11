Subroutine rdtpsi(lunit, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
& wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam, itimrev)
!-----------------------------------------------------------------------
!        read stored wave functions
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: lunit
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: itimrev
      Integer :: npsi(2)
      Integer :: npmin(2)
      Real(wp) :: vocc(lpsi, 2)
      Real(wp) :: xcm(2)
      Real(wp) :: ycm(2)
      Real(wp) :: zcm(2)
      Real(wp) :: spenrg(lpsi, 2)
      Real(wp) :: spnorm(lpsi, 2)
      Real(wp) :: spkine(lpsi, 2)
      Real(wp) :: ajz(lpsi, 2)
      Real(wp) :: spar(lpsi, 2)
      Real(wp) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp) :: xlam
      Real(wp) :: gapn, gapp
      Real(wp), Dimension(2) :: fermin, fermip, epairn, epairp
      Complex(wp) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, i, ix, iy, iz, nst, is, lpsis
!-----------------------------------------------
      If(itimrev == 1) Then
         lpsis = lpsi/2
      Else
         lpsis = lpsi
      End If
      Do iq = 1, 2
         Read(lunit) npmin(iq), npsi(iq), xcm(iq), ycm(iq), zcm(iq)
         Read(lunit) (vocc(i, iq), i=1, lpsis)
         Read(lunit) (ajz(i, iq), i=1, lpsis)
         Read(lunit) (spenrg(i, iq), i=1, lpsis)
         Read(lunit) (spar(i, iq), i=1, lpsis)
         Read(lunit) (spnorm(i, iq), i=1, lpsis)
         Read(lunit) (spkine(i, iq), i=1, lpsis)
      End Do
!
      Read(lunit) gapn, gapp, fermin, fermip, epairn, epairp, xlam
      Read(lunit) (((wcoul(ix, iy, iz), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
!
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            Do is = 1, 2
               Read(lunit) (((psi(ix, iy, iz, is, nst, iq), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
            End Do
         End Do
      End Do
!
      Return
End Subroutine rdtpsi
