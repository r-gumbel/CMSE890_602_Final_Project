Subroutine getpsi(psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, wcoul, &
& mnof, tetrs, centf, npsif, npminf, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul)
!-----------------------------------------------------------------------
!        read wave functions from unit 14
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
      Integer :: npsi(2)
      Integer :: npmin(2)
      Real(wp) :: tetrs
      Real(wp) :: rold
      Real(wp) :: tetold
      Real(wp) :: rscale
      Real(wp) :: q20old
      Real(wp) :: q20tot
      Real(wp) :: xdot
      Real(wp) :: zdot
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
      Real(wp) :: ecoul(mnof)
      Real(wp), Dimension(3, mnof) :: centf
      Integer, Dimension(2, mnof) :: npsif, npminf
      Complex(wp) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, i, ix, iy, iz, nst, is
!-----------------------------------------------
      Read(14) tetrs
      Read(14) rold, tetold, rscale, q20old, q20tot, xdot, zdot
      Do iq = 1, 2
         Read(14) npmin(iq), npsi(iq), xcm(iq), ycm(iq), zcm(iq)
         Read(14) (vocc(i, iq), i=1, lpsi)
         Read(14) (ajz(i, iq), i=1, lpsi)
         Read(14) (spenrg(i, iq), i=1, lpsi)
         Read(14) (spar(i, iq), i=1, lpsi)
         Read(14) (spnorm(i, iq), i=1, lpsi)
         Read(14) (spkine(i, iq), i=1, lpsi)
      End Do
      Read(14) centf(1:3, 1:mnof)
      Read(14) npsif(1:2, 1:mnof)
      Read(14) npminf(1:2, 1:mnof)
      Read(14) ecoul(1:mnof)
      Read(14) (((wcoul(ix, iy, iz), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
!
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            Do is = 1, 2
               Read(14) (((psi(ix, iy, iz, is, nst, iq), ix=1, ncolx), iy=1, ncoly), iz=1, ncolz)
            End Do
         End Do
      End Do
!
      Return
End Subroutine getpsi
