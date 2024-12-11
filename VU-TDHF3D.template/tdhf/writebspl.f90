Subroutine writebspl(psi, binvx, binvy, binvz, nx, ny, nz, ncolx, ncoly, ncolz, &
& lpsi, npsi, npmin, nordx, nordy, nordz, etheta, xcm, ycm, zcm, xclx, xcly, xclz, vocc)
!
!-------------------------------------------------------------------------------
!       subroutine to rotate the nucleus by rotating the coordinate system and
!       interpolating followed by a rotation in spin-space.
!-------------------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: nx
      Integer, Intent(In) :: ny
      Integer, Intent(In) :: nz
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: nordx
      Integer, Intent(In) :: nordy
      Integer, Intent(In) :: nordz
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: binvx(nx, nx)
      Real(wp), Intent(In) :: binvy(ny, ny)
      Real(wp), Intent(In) :: binvz(nz, nz)
      Real(wp), Intent(In) :: etheta(3)
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(In) :: vocc(lpsi, 2)
      Complex(wp) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
      Real(wp), Parameter :: small = 1.0e-25_wp
      Integer :: iq, nst, is, jz, jy, jx, iz, iy, ix
!     The bounds of bspl_coef should be the number of splines nx, ny, nz.
!     With periodic boundary conditions, these are equal to the number
!     of collocation points ncolx, ncoly, ncolz.
      Complex(wp) :: bspl_coef(nx, ny, nz, 2, lpsi, 2)
!-------------------------------------------------------------------------
!        Calculate the B-spline expansion coefficients
!        of each single-particle state. These will then be used to
!        interpolate the wavefunction onto a new mesh.
!        The size of the coefficient array is different for
!        periodic bc then fixed bc, the larger being the fixed bc.
!-------------------------------------------------------------------------
      bspl_coef = (0.0_wp, 0.0_wp)
!
      Do iq = 1, 2
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:bspl_coef)
         Do nst = npmin(iq), npsi(iq)
            Do is = 1, 2
!
               Do jz = 1, ncolz
                  Do jy = 1, ncoly
                     Do jx = 1, ncolx
!
                        Do iz = 1, nz
                           Do iy = 1, ny
                              Do ix = 1, nx
                                 bspl_coef(ix, iy, iz, is, nst, iq) = bspl_coef(ix, iy, iz, is, nst, iq) + binvz(iz, jz) * &
                                & binvy(iy, jy) * binvx(ix, jx) * psi(jx, jy, jz, is, nst, iq)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
!
            End Do
         End Do
!$OMP END PARALLEL DO
      End Do
!
      Rewind(19)
      Write(19) nordx, nordy, nordz
      Write(19) ncolx, ncoly, ncolz
      Write(19) xclx, xcly, xclz
      Write(19) etheta
      Write(19) xcm, ycm, zcm
      Write(19) vocc(npmin(1):npsi(1),1), vocc(npmin(2):npsi(2),2)
      Write(19) bspl_coef(:, :, :, :, npmin(1):npsi(1), 1)
      Write(19) bspl_coef(:, :, :, :, npmin(2):npsi(2), 2)
      Return
End Subroutine writebspl
