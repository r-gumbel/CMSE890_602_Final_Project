Subroutine putpsi(psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, wcoul, &
& mnof, tetrs, centf, npsif, npminf, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul)
!-----------------------------------------------------------------------
!        write collocation wave functions to unit 14
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: mnof
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Real(wp), Intent(In) :: tetrs
      Real(wp), Intent(In) :: rold
      Real(wp), Intent(In) :: tetold
      Real(wp), Intent(In) :: rscale
      Real(wp), Intent(In) :: q20old
      Real(wp), Intent(In) :: q20tot
      Real(wp), Intent(In) :: xdot
      Real(wp), Intent(In) :: zdot
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
      Real(wp), Intent(In) :: ecoul(mnof)
      Real(wp), Dimension(3, mnof) :: centf
      Integer, Dimension(2, mnof) :: npsif, npminf
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, i, ix, iy, iz, nst, is
!-----------------------------------------------
      Rewind(14)
!
      Write(14) tetrs
      Write(14) rold, tetold, rscale, q20old, q20tot, xdot, zdot
      Do iq = 1, 2
         Write(14) npmin(iq), npsi(iq), xcm(iq), ycm(iq), zcm(iq)
         Write(14) (vocc(i, iq), i=1, lpsi)
         Write(14) (ajz(i, iq), i=1, lpsi)
         Write(14) (spenrg(i, iq), i=1, lpsi)
         Write(14) (spar(i, iq), i=1, lpsi)
         Write(14) (spnorm(i, iq), i=1, lpsi)
         Write(14) (spkine(i, iq), i=1, lpsi)
      End Do
      Write(14) centf(1:3, 1:mnof)
      Write(14) npsif(1:2, 1:mnof)
      Write(14) npminf(1:2, 1:mnof)
      Write(14) ecoul(1:mnof)
      Write(14) (((wcoul(ix, iy, iz), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
!
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            Do is = 1, 2
               Write(14) (((psi(ix, iy, iz, is, nst, iq), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
            End Do
         End Do
      End Do
!
      Call flush(14)
      Return
End Subroutine putpsi
