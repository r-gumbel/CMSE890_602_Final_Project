Subroutine pois(iperc, coulplan1, coulplan2, ncolx, ncoly, ncolz, nx2, ny2, nz2, dx, dy, dz, rho, rho2, q, wcoul)
!-----------------------------------------------------------------------
!       solve the poisson equation to obtain the coulomb potential on
!       the 3d lattice using Green's function FFT method.
!                  -d**2*wc=4*pi*e2*rho
!       The method is Hockney's domain doubling - FFT approach for
!       isolated charge distributions (ordinary HF and TDHF) and
!       standard method for periodic charge systems (neutron star
!       crust) assumed to be neutralized by the background electron
!       gas (remove monopole term).
!       output:
!           wcoul : coulomb potential at lattice points.
!       NOTE: Subroutines use the fftw3 library and the calling routines
!             set the dimension of rho2 and q depending on the type of
!             boundary condition.
!-----------------------------------------------------------------------
      Use, intrinsic :: iso_c_binding
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer(8), Intent(In) :: coulplan1, coulplan2 ! FFTW plans for Coulomb
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: nx2
      Integer, Intent(In) :: ny2
      Integer, Intent(In) :: nz2
      Integer, Intent(In) :: iperc
      Real(wp), Intent(In) :: dx
      Real(wp), Intent(In) :: dy
      Real(wp), Intent(In) :: dz
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz)
      Real(wp), Intent(Out) :: wcoul(ncolx, ncoly, ncolz)
      Complex(wp), Intent(In) :: q(nx2, ny2, nz2)
      Complex(wp) :: rho2(nx2, ny2, nz2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: e2 = 1.439978408596513_wp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ix, iy, iz
! put density into array of double size (periodic)
      Do iz = 1, nz2
         Do iy = 1, ny2
            Do ix = 1, nx2
               If(iz > ncolz .Or. iy > ncoly .Or. ix > ncolx) Then
                  rho2(ix, iy, iz) = 0.0_wp
               Else
                  rho2(ix, iy, iz) = rho(ix, iy, iz)
               End If
            End Do
         End Do
      End Do
! get the FFT of the new rho
      Call dfftw_execute_dft(coulplan1, rho2, rho2)
! fold and renormalize, also multiply by e2
      If(iperc == 1) Then
         rho2 = 4.0_wp * pi * e2 * REAL(q) * rho2
      Else
         rho2 = q * rho2 * (dx*dy*dz) * e2
      End If
! back transform to coordinate space
      Call dfftw_execute_dft(coulplan2, rho2, rho2)
! assign coulomb potential to real space
      wcoul = REAL(rho2(1:ncolx, 1:ncoly, 1:ncolz)) / (nx2*ny2*nz2)
!
End Subroutine pois
Subroutine coulinit(coulplan1, coulplan2, nx2, ny2, nz2, dx, dy, dz, q, iperc)
!-----------------------------------------------------------------------
!       initialize Fourier transform of the Greens function on mesh.
!       The method is Hockney's domain doubling - FFT approach for
!       isolated charge distributions (ordinary HF and TDHF) and
!       standard method for periodic charge systems (neutron star
!       crust) assumed to be neutralized by the background electron
!       gas (remove monopole term).
!-----------------------------------------------------------------------
      Use, intrinsic :: iso_c_binding
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
! include fftw3 library headers
      Include 'fftw3.f03'
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer(8), Intent(In) :: coulplan1, coulplan2 ! FFTW plans for Coulomb
      Integer, Intent(In) :: nx2
      Integer, Intent(In) :: ny2
      Integer, Intent(In) :: nz2
      Integer, Intent(In) :: iperc
      Real(wp), Intent(In) :: dx, dy, dz
      Complex(wp), Intent(Out) :: q(nx2, ny2, nz2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp), Allocatable :: iqx(:), iqy(:), iqz(:)
      Integer :: ix, iy, iz
!
      Allocate(iqx(nx2), iqy(ny2), iqz(nz2))
!
      Call dfftw_plan_dft_3d(coulplan1, nx2, ny2, nz2, q, q, FFTW_FORWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
      Call dfftw_plan_dft_3d(coulplan2, nx2, ny2, nz2, q, q, FFTW_BACKWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
!
      Call initiq(dx, nx2, iqx, iperc)
      Call initiq(dy, ny2, iqy, iperc)
      Call initiq(dz, nz2, iqz, iperc)
!
      Do iz = 1, nz2
         Do iy = 1, ny2
            Do ix = 1, nx2
               q(ix, iy, iz) = iqx(ix) + iqy(iy) + iqz(iz)
            End Do
         End Do
      End Do
      If(iperc == 1) Then
         q(1, 1, 1) = (1.0_wp, 0.0_wp)
         q = 1.0_wp / REAL(q)
! remove monopole term for charge neutral symmetric systems
         q(1, 1, 1) = (0.0_wp, 0.0_wp)
      Else
         q(1, 1, 1) = (1.0_wp, 0.0_wp)
         q = 1.0_wp / Sqrt(REAL(q))
         q(1, 1, 1) = 2.92_wp*3.0_wp/(dx+dy+dz)  !set q(0) to prescribed value
! get the FFT of the green's function
         Call dfftw_execute_dft(coulplan1, q, q)
      End If
!
      Deallocate(iqx, iqy, iqz)
!
End Subroutine coulinit
Subroutine initiq(d, n, iq, iperc)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: n
      Integer, Intent(In) :: iperc
      Real(wp), Intent(In) :: d
      Real(wp), Intent(Out) :: iq(n)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: i, ii
!
      Do i = 1, n
         If(i <= n/2) Then
            ii = i - 1
         Else
            ii = i - n - 1
         End If
         If(iperc == 1) Then
            iq(i) = (2.0_wp*pi*ii/(n*d)) ** 2
         Else
            iq(i) = (d*ii) ** 2
         End If
      End Do
!
End Subroutine initiq
