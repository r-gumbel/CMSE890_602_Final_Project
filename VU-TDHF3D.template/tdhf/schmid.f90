Subroutine schmid(lpsi, ncolx, ncoly, ncolz, npsi, npmin, spnorm, psi, tpsi, tpsir, trevpsi, itimrev, wxyz)
!-----------------------------------------------------------------------
!        gram-schmidt orthogonalization of single-particle states
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp), Intent(Out) :: spnorm(lpsi, 2)
      Complex(wp) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp) :: tpsi(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: tpsir(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: trevpsi(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: nst, npart, j
      Complex(wp) :: oij
!-----------------------------------------------------------------------
!        begin the gram-schmidt procedure - neutrons
!-----------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,oij)
      Do nst = npmin(1), npsi(1)
         Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, 1), psi(1, 1, 1, 1, nst, 1), oij, 0, wxyz)
      End Do
!$OMP END PARALLEL DO
      npart = npsi(1) - npmin(1) + 1
      If(npart /= 1) Then
         Do nst = npmin(1), npsi(1)
            tpsi = (0.0_wp,0.0_wp)
            tpsir = (0.0_wp,0.0_wp)
!$OMP PARALLEL DO PRIVATE(j,oij) REDUCTION(+:tpsi)
            Do j = npmin(1), nst - 1
               Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, j, 1), psi(1, 1, 1, 1, nst, 1), oij, 1, wxyz)
               tpsi(:,:,:,:) = tpsi(:,:,:,:) + oij*psi(:, :, :, :, j, 1)
            End Do
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------
!     find the time reversed states and also orthogonalize
!-----------------------------------------------------------------------
            If(itimrev /= 0) Then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,oij,trevpsi) REDUCTION(+:tpsir)
               Do j = npmin(1), nst - 1
                  Call trpsi(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, j, 1), trevpsi(1, 1, 1, 1))
                  Call psnorm(ncolx, ncoly, ncolz, trevpsi(1, 1, 1, 1), psi(1, 1, 1, 1, nst, 1), oij, 1, wxyz)
                  tpsir(:,:,:,:) = tpsir(:,:,:,:) + oij*trevpsi(:,:,:,:)
               End Do
!$OMP END PARALLEL DO
            End If
            psi(:, :, :, :, nst, 1) = psi(:, :, :, :, nst, 1) - tpsi(:,:,:,:) - tpsir(:,:,:,:)
            Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, 1), psi(1, 1, 1, 1, nst, 1), oij, 0, wxyz)
            spnorm(nst, 1) = real(oij)
         End Do
      End If
!-----------------------------------------------------------------------
!        begin the gram-schmidt procedure - protons
!-----------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,oij)
      Do nst = npmin(2), npsi(2)
         Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, 2), psi(1, 1, 1, 1, nst, 2), oij, 0, wxyz)
      End Do
!$OMP END PARALLEL DO
      npart = npsi(2) - npmin(2) + 1
      If(npart /= 1) Then
         Do nst = npmin(2), npsi(2)
            tpsi = (0.0_wp,0.0_wp)
            tpsir = (0.0_wp,0.0_wp)
!$OMP PARALLEL DO PRIVATE(j,oij) REDUCTION(+:tpsi)
            Do j = npmin(2), nst - 1
               Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, j, 2), psi(1, 1, 1, 1, nst, 2), oij, 1, wxyz)
               tpsi(:,:,:,:) = tpsi(:,:,:,:) + oij*psi(:, :, :, :, j, 2)
            End Do
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------
!     find the time reversed states and also orthogonalize
!-----------------------------------------------------------------------
            If(itimrev /= 0) Then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,oij,trevpsi) REDUCTION(+:tpsir)
               Do j = npmin(2), nst - 1
                  Call trpsi(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, j, 2), trevpsi(1, 1, 1, 1))
                  Call psnorm(ncolx, ncoly, ncolz, trevpsi(1, 1, 1, 1), psi(1, 1, 1, 1, nst, 2), oij, 1, wxyz)
                  tpsir(:,:,:,:) = tpsir(:,:,:,:) + oij*trevpsi(:,:,:,:)
               End Do
!$OMP END PARALLEL DO
            End If
            psi(:, :, :, :, nst, 2) = psi(:, :, :, :, nst, 2) - tpsi(:,:,:,:) - tpsir(:,:,:,:)
            Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, 2), psi(1, 1, 1, 1, nst, 2), oij, 0, wxyz)
            spnorm(nst, 2) = real(oij)
         End Do
      End If
      Return
End Subroutine schmid
