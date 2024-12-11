Subroutine bustf(lpsi, ncolx, ncoly, ncolz, npsif, npminf, mnof, xclx, xcly, xclz, boostf, centf, fmass, h2m, psi, iff)
!-----------------------------------------------------------------------
!        boosts the initial wavefunctions for dynamic calculation
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: mnof
      Integer, Intent(In) :: iff
      Integer, Intent(In) :: npsif(2, mnof)
      Integer, Intent(In) :: npminf(2, mnof)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: boostf(3, mnof), centf(3, mnof)
      Real(wp), Intent(In) :: fmass(mnof)
      Real(wp), Intent(In) :: h2m(2)
      Complex(wp), Intent(Inout) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, nst, is, iz, iy, ix
      Real(wp) :: sbx, sby, sbz, tbx, tby, tbz, akfx, akfy, akfz, z, y, x
!-----------------------------------------------------------------------
!        calculate the wavenumber for each single-particle state from
!        the array boostf andloop over all states and multiply with the
!        boost.
!-----------------------------------------------------------------------
      sbx = 0.0_wp
      sby = 0.0_wp
      sbz = 0.0_wp
      tbx = boostf(1, iff)
      tby = boostf(2, iff)
      tbz = boostf(3, iff)
      If(tbx /= 0.0_wp) sbx = tbx / Abs(tbx)
      If(tby /= 0.0_wp) sby = tby / Abs(tby)
      If(tbz /= 0.0_wp) sbz = tbz / Abs(tbz)
!
      Do iq = 1, 2
         akfx = sbx * Sqrt(1.0_wp/h2m(iq)*Abs(tbx)/fmass(iff))
         akfy = sby * Sqrt(1.0_wp/h2m(iq)*Abs(tby)/fmass(iff))
         akfz = sbz * Sqrt(1.0_wp/h2m(iq)*Abs(tbz)/fmass(iff))
!
         Do nst = npminf(iq, iff), npminf(iq, iff) + npsif(iq, iff) - 1
!
            Do is = 1, 2
               Do iz = 1, ncolz
                  z = xclz(iz) - centf(3, iff)
                  Do iy = 1, ncoly
                     y = xcly(iy) - centf(2, iff)
                     Do ix = 1, ncolx
                        x = xclx(ix) - centf(1, iff)
                        psi(ix, iy, iz, is, nst, iq) = psi(ix, iy, iz, is, nst, iq) * Exp(eye*akfx*x) * Exp(eye*akfy*y) * &
                       & Exp(eye*akfz*z)
                     End Do
                  End Do
               End Do
            End Do
!
         End Do
      End Do
!
      Return
End Subroutine bustf
