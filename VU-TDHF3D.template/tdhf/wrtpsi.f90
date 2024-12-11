Subroutine wrtpsi(lunit, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
& wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam)
!-----------------------------------------------------------------------
!   write collocation wave functions
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: lunit
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Real(wp), Intent(In) :: vocc(lpsi, 2)
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(In) :: spenrg(lpsi, 2)
      Real(wp), Intent(In) :: spnorm(lpsi, 2)
      Real(wp), Intent(In) :: spkine(lpsi, 2)
      Real(wp), Intent(In) :: ajz(lpsi, 2)
      Real(wp), Intent(In) :: spar(lpsi, 2)
      Real(wp), Intent(In) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: xlam
      Real(wp), Intent(In) :: gapn, gapp
      Real(wp), Intent(In), Dimension(2) :: fermin, fermip, epairn, epairp
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, i, ix, iy, iz, nst, is
!-----------------------------------------------
      If(lunit == 12) Then
         Rewind(12)
      End If
!
      Do iq = 1, 2
         Write(lunit) npmin(iq), npsi(iq), xcm(iq), ycm(iq), zcm(iq)
         Write(lunit) (vocc(i, iq), i=1, lpsi)
         Write(lunit) (ajz(i, iq), i=1, lpsi)
         Write(lunit) (spenrg(i, iq), i=1, lpsi)
         Write(lunit) (spar(i, iq), i=1, lpsi)
         Write(lunit) (spnorm(i, iq), i=1, lpsi)
         Write(lunit) (spkine(i, iq), i=1, lpsi)
      End Do
      Write(lunit) gapn, gapp, fermin, fermip, epairn, epairp, xlam
      Write(lunit) (((wcoul(ix, iy, iz), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
!
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            Do is = 1, 2
               Write(lunit) (((psi(ix, iy, iz, is, nst, iq), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
            End Do
         End Do
      End Do
!
      Call flush(12)
      Call flush(15)
      Return
End Subroutine wrtpsi
