Program tdhf3d
!-----------------------------------------------------------------------
!    VU-TDHF3D
!
!    3-D Time-Dependent Hartree-Fock program@ in Fortran 95 written by A.S. Umar.
!
!    This code includes the entire Skyrme interaction, including all of the
!    time-odd terms as described in Engel et al.
!
!    Use of this program without the consent of the Vanderbilt Computational
!    Nuclear Theory Group would be a violation of scientific ethical rules.
!
!    The main papers describing the details of various parts of this code:
!
!    Engel, Brink, Goeke, Krieger, and Vautherin, Nucl. Phys. A249, 215 (1975).
!    Basis-Splines   : J. Comp. Phys. 93, 426 (1991).
!    3D HF           : Phys. Rev. C44, 2512 (1991).
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      Use spline_grid
      Use, intrinsic :: iso_c_binding
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Implicit None
      Integer, Parameter  :: wp = Kind(1.0D0)
      Include 'fftw3.f03'
      Integer, Parameter  :: mnof = 3 ! maximum number of nuclei
      Integer, Parameter  :: nfixb = 0 ! fix boost, use with caution for L violation
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: hbc = 197.3269885804381_wp
      Real(wp), Parameter :: e2 = 1.439978408596513_wp
!-----------------------------------------------------------------------
!        dimension statements for the spline routines
!-----------------------------------------------------------------------
! x-
      Real(wp), Dimension(:), Pointer :: xclx_tmp
      Real(wp), Dimension(:), Pointer :: xkvx_tmp
      Real(wp), Dimension(:, :), Pointer :: coefrx_tmp
      Real(wp), Dimension(:, :), Pointer :: coeflx_tmp
      Real(wp), Dimension(:), Allocatable :: xkvx
      Real(wp), Dimension(:), Allocatable :: xclx, wx
      Real(wp), Dimension(:, :), Allocatable :: coefrx, coeflx
      Real(wp), Dimension(:, :), Allocatable :: der1x, der2x, binvpx
! y-
      Real(wp), Dimension(:), Pointer :: xcly_tmp
      Real(wp), Dimension(:), Pointer :: xkvy_tmp
      Real(wp), Dimension(:, :), Pointer :: coefry_tmp
      Real(wp), Dimension(:, :), Pointer :: coefly_tmp
      Real(wp), Dimension(:), Allocatable :: xkvy
      Real(wp), Dimension(:), Allocatable :: xcly, wy
      Real(wp), Dimension(:, :), Allocatable :: coefry, coefly
      Real(wp), Dimension(:, :), Allocatable :: der1y, der2y, binvpy
! z-
      Real(wp), Dimension(:), Pointer :: xclz_tmp
      Real(wp), Dimension(:), Pointer :: xkvz_tmp
      Real(wp), Dimension(:, :), Pointer :: coefrz_tmp
      Real(wp), Dimension(:, :), Pointer :: coeflz_tmp
      Real(wp), Dimension(:), Allocatable :: xkvz
      Real(wp), Dimension(:), Allocatable :: xclz, wz
      Real(wp), Dimension(:, :), Allocatable :: coefrz, coeflz
      Real(wp), Dimension(:, :), Allocatable :: der1z, der2z, binvpz
!-----------------------------------------------------------------------
!        storage for wavefunctions and quantum numbers etc.
!-----------------------------------------------------------------------
      Complex(wp), Dimension(:, :, :, :, :, :), Allocatable :: psi, psitmp
      Complex(wp), Dimension(:, :, :, :), Allocatable :: pswk1, pswk2
      Complex(wp), Dimension(:, :, :, :), Allocatable :: pswk3, pswk4, tpsi
      Real(wp), Dimension(:, :), Allocatable   :: vocc, alz, asz, ajz, spar
      Real(wp), Dimension(:, :), Allocatable   :: vocct, alzt, aszt, ajzt, spart
      Real(wp), Dimension(:, :), Allocatable   :: spenrg, spkine, spnorm, speflu
      Real(wp), Dimension(:, :), Allocatable   :: spenrgt, spkinet, spnormt, speflut
      Real(wp), Dimension(2)                   :: fermin, fermip, epairn, epairp
      Integer, Dimension(:, :, :), Allocatable :: nshell
      Integer, Dimension(:, :), Allocatable    :: listis
      Integer, Dimension(2)                    :: npsi, npmin, npsitot
      Integer, Dimension(2, mnof)              :: npsif, npminf, ipairf
      Logical, Dimension(mnof)                 :: isrestart
!-----------------------------------------------------------------------
!     storage for fields,currents and densities
!
!     mean-field potentials:
!        upot    = local potential
!        bmass   = effective mass
!        dbmass  = gradient of the effective mass
!        xiq     = the current I_q, used in (del.I+I.del) term
!        dxiq    = divergence of xiq
!        cq      = time-odd C_q term
!        dcq     = gradient tensor of cq
!        bmunu   = B-tensor
!        dbmu    = divergence of B-tensor
!        hsigma  = time-odd upper case sigma_q
!        wcoul   = coulomb potential
!        wyuk    = Yukawa potential for BKN force
!        wctot   = used to save coulomb for initial guesses
!
!     densities and currents:
!       rho     = density (for protons and neutrons separately)
!       tau     = kinetic energy density
!       currnt  = current (is a space vector)
!       sodens  = spin-orbit density (is a space vector)
!       spinden = spin-density s
!       kinvden = kinetic energy density (vector part)
!       spincur = spin-current tensor
!-----------------------------------------------------------------------
      Real(wp), Dimension(:, :, :), Allocatable :: wcoul, wctot
      Real(wp), Dimension(:, :, :), Allocatable :: edenss, edenst
      Real(wp), Dimension(:, :, :), Allocatable :: wxyz, wyuk
      Real(wp), Dimension(:, :, :, :), Allocatable :: upot, bmass
      Real(wp), Dimension(:, :, :, :), Allocatable :: dxiq, rho, tau
      Real(wp), Dimension(:, :, :, :), Allocatable :: worka, rlam
      Real(wp), Dimension(:, :, :, :, :), Allocatable :: xiq, cq, hsigma, tmpvec
      Real(wp), Dimension(:, :, :, :, :), Allocatable :: dbmass
      Real(wp), Dimension(:, :, :, :, :), Allocatable :: currnt, sodens
      Real(wp), Dimension(:, :, :, :, :), Allocatable :: spinden, kinvden, dbmu
      Real(wp), Dimension(:, :, :, :, :, :), Allocatable :: spincur, delxj, bmunu, dcq
      Real(wp), Dimension(2) :: xcm, ycm, zcm
      Real(wp), Dimension(2) :: h2m
      Complex(wp), Dimension(:, :, :), Allocatable :: q, rho2
!-----------------------------------------------------------------------
!        local work arrays
!-----------------------------------------------------------------------
      Real(wp), Dimension(:, :), Allocatable        :: wsk1x, wsk2x
      Real(wp), Dimension(:, :), Allocatable        :: wsk1y, wsk2y
      Real(wp), Dimension(:, :), Allocatable        :: wsk1z, wsk2z
      Real(wp), Dimension(:,:,:,:,:), Allocatable   :: rhos, taus
      Real(wp), Dimension(:,:,:,:,:,:), Allocatable :: drhos, currnts
!-----------------------------------------------------------------------
!        storage for dynamic calculations
!-----------------------------------------------------------------------
      Real(wp), Dimension(3, mnof) :: centf, boostf, radinf
      Real(wp), Dimension(mnof) :: fmass, fcharg, euler_alpha, euler_beta, euler_gamma
!-----------------------------------------------------------------------
!        constants
!-----------------------------------------------------------------------
      Real(wp) :: xin, xout, yin, yout, zin, zout
      Real(wp) :: t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk
      Real(wp) :: v0neut, v0prot, rho0pr
      Real(wp) :: radinx, radiny, radinz
      Real(wp) :: serr, derr, x0dmp, e0dmp, c0, d0, q2in
      Real(wp) :: ecm, rsep, xli, terr, dt
      Real(wp) :: amass, zmass, cmfact
      Real(wp) :: amc2, xmu, ratio, akv, xb, rstr, xlg, epsi, ttt, tetr, cofl, cofr
      Real(wp) :: delx, dely, delz, roft0, roft, xcoul, xcent, tetc, tetf, tets
      Real(wp) :: tdot, ecm0, rdot, vx, vz, tcmx, tcmz, tketot, vlx, vlz, vrx, vrz
      Real(wp) :: tlx, tlz, trx, trz
      Real(wp) :: tlix, tliz, trix, triz, epslx, epslz, epsrx, epsrz
      Real(wp) :: atheta, btheta, gtheta, rmin, etheta(3)
      Real(wp) :: ecoultot, ecoul(mnof), ecx, xdt
!
      Integer :: nnotx, ncolx, nordx, nord1x, nsplx, igrdx, ianlyx, iperx
      Integer :: nnoty, ncoly, nordy, nord1y, nsply, igrdy, ianlyy, ipery
      Integer :: nnotz, ncolz, nordz, nord1z, nsplz, igrdz, ianlyz, iperz
      Integer :: iperc, iunit
      Integer :: nx2, ny2, nz2
      Integer :: ix, iy, iz
      Integer :: itrbx, mtrbx, iprint, mprint, mplots, mplott, icoul, irest, mrests, mrestt
      Integer :: itimrev, itimrevs, iodds, ihdiag
      Integer :: nneut, nprot, iconstr, iconstr0, mconstr, itheta, lpsi, imode
      Integer :: ipairn, ipairp, ipairn0, ipairp0, imode0
      Integer :: nt, nof, mxp, ifixb, irwgs, irwgs0
      Integer :: inuc, ifixcm, npsi1, npsi2, niter, nbloop
      Integer :: iq, nst, nst1, nst2
      Integer :: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
      Integer :: nprocs, maxthreads !, nfftw
      Integer :: nexadd, nextra_n(mnof), nextra_p(mnof), nexiter, n1, n2
      Integer(8) :: coulplan1, coulplan2 ! FFTW plans for Coulomb
      Character(len=80) :: directname
!-----------------------------------------------------------------------
!        open files for input
!-----------------------------------------------------------------------
      Open(Unit=3, File='bspl.inp', Status='old', Form='formatted', Position='asis')
      Open(Unit=5, File='tdhf3d.inp', Status='old', Form='formatted', Position='asis')
      Open(Unit=10, File='pair.inp', Status='old', Form='formatted', Position='asis')
      Open(Unit=11, File='skyrme.inp', Status='old', Form='formatted', Position='asis')
!-----------------------------------------------------------------------
!        call getin to read data initializing the calculation
!-----------------------------------------------------------------------
      Call getin(h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, v0neut, v0prot, rho0pr, &
     & radinf, itrbx, mtrbx, serr, derr, x0dmp, e0dmp, iprint, mprint, mplots, mplott, icoul, irest, mrests, &
     & mrestt, imode, ifixcm, itimrev, itimrevs, ipairf, c0, d0, iconstr, q2in, mconstr, itheta, nt, dt, mnof, nof, &
     & centf, boostf, euler_alpha, euler_beta, euler_gamma, fmass, fcharg, ecm, rsep, xb, mxp, terr, ifixb, irwgs, &
     & iodds, ihdiag, nexadd, nexiter, nextra_n, nextra_p, directname, isrestart)
!------------------------------------------------------------------------
!        setup knots and collocation points
!        NOTE: Here we make use of pointers to determine the yet unknown
!              dimension of the collocation/knot arrays. This is done via
!              the module spline_grid in the USE statement at the top.
!------------------------------------------------------------------------
! x-
      Call grid(coeflx_tmp, coefrx_tmp, xkvx_tmp, xclx_tmp, xin, xout, nnotx, ncolx, nordx, nord1x, nsplx, igrdx, ianlyx, iperx)
!
      Allocate(xkvx(nnotx), xclx(ncolx))
      xkvx = xkvx_tmp
      xclx = xclx_tmp
      Allocate(coeflx(nordx, nordx), coefrx(nordx, nordx))
      coeflx = coeflx_tmp
      coefrx = coefrx_tmp
      Nullify(xkvx_tmp, xclx_tmp, coeflx_tmp, coefrx_tmp)
! y-
      Call grid(coefly_tmp, coefry_tmp, xkvy_tmp, xcly_tmp, yin, yout, nnoty, ncoly, nordy, nord1y, nsply, igrdy, ianlyy, ipery)
!
      Allocate(xkvy(nnoty), xcly(ncoly))
      xkvy = xkvy_tmp
      xcly = xcly_tmp
      Allocate(coefly(nordy, nordy), coefry(nordy, nordy))
      coefly = coefly_tmp
      coefry = coefry_tmp
      Nullify(xkvy_tmp, xcly_tmp, coefly_tmp, coefry_tmp)
! z-
      Call grid(coeflz_tmp, coefrz_tmp, xkvz_tmp, xclz_tmp, zin, zout, nnotz, ncolz, nordz, nord1z, nsplz, igrdz, ianlyz, iperz)
!
      Allocate(xkvz(nnotz), xclz(ncolz))
      xkvz = xkvz_tmp
      xclz = xclz_tmp
      Allocate(coeflz(nordz, nordz), coefrz(nordz, nordz))
      coeflz = coeflz_tmp
      coefrz = coefrz_tmp
      Nullify(xkvz_tmp, xclz_tmp, coeflz_tmp, coefrz_tmp)
!-----------------------------------------------------------------------
!        allocate arrays for the spline derivaties and weights
!-----------------------------------------------------------------------
      Allocate(wx(ncolx))
      Allocate(der1x(ncolx, ncolx), der2x(ncolx, ncolx), binvpx(ncolx, ncolx))
      Allocate(wy(ncoly))
      Allocate(der1y(ncoly, ncoly), der2y(ncoly, ncoly), binvpy(ncoly, ncoly))
      Allocate(wz(ncolz))
      Allocate(der1z(ncolz, ncolz), der2z(ncolz, ncolz), binvpz(ncolz, ncolz))
!-----------------------------------------------------------------------
!        allocate local work arrays
!-----------------------------------------------------------------------
      Allocate(wsk1x(ncolx, ncolx), wsk2x(ncolx, ncolx))
      Allocate(wsk1y(ncoly, ncoly), wsk2y(ncoly, ncoly))
      Allocate(wsk1z(ncolz, ncolz), wsk2z(ncolz, ncolz))
      Allocate(rhos(ncolx,ncoly,ncolz,2,2), taus(ncolx,ncoly,ncolz,2,2) )
      Allocate(drhos(ncolx,ncoly,ncolz,3,2,2), currnts(ncolx,ncoly,ncolz,3,2,2) )
!-----------------------------------------------------------------------
!        setup splines etc. for the general problem
!-----------------------------------------------------------------------
      Call spline(der1x, der2x, binvpx, xclx, xkvx, coeflx, coefrx, wx, xin, xout, ncolx, nordx, nnotx, nsplx, igrdx,&
      &ianlyx, iperx, 'x')
      Call spline(der1y, der2y, binvpy, xcly, xkvy, coefly, coefry, wy, yin, yout, ncoly, nordy, nnoty, nsply, igrdy,&
      &ianlyy, ipery, 'y')
      Call spline(der1z, der2z, binvpz, xclz, xkvz, coeflz, coefrz, wz, zin, zout, ncolz, nordz, nnotz, nsplz, igrdz,&
      &ianlyz, iperz, 'z')
!-----------------------------------------------------------------------
!     storage for fields,currents and densities
!
!     mean-field potentials:
!        upot    = local potential
!        bmass   = effective mass
!        dbmass  = gradient of the effective mass
!        xiq     = the current I_q, used in (del.I+I.del) term
!        dxiq    = divergence of xiq
!        cq      = time-odd C_q term
!        dcq     = gradient tensor of cq
!        bmunu   = B-tensor
!        dbmu    = divergence of B-tensor
!        hsigma  = time-odd upper case sigma_q
!        wcoul   = coulomb potential
!        wyuk    = Yukawa potential for BKN
!        wctot   = used to save coulomb for initial guesses
!
!     densities and currents:
!       rho     = density (for protons and neutrons separately)
!       tau     = kinetic energy density
!       currnt  = current (is a space vector)
!       sodens  = spin-orbit density (is a space vector)
!       spinden = spin density (is a space vector)
!       kinvden = kinetic energy density (vector part)
!       spincur = spin-current tensor (has 9 components)
!
!     the 2 in the dimensions is the isospin (1=neutron, 2=proton),
!     the 3 in the dimensions denote the components of a space vector
!     the 3,3 in the dimensions indicate a 3x3 tensor
!-----------------------------------------------------------------------
      Allocate(wcoul(ncolx, ncoly, ncolz))
      Allocate(wyuk(ncolx, ncoly, ncolz))
      Allocate(wxyz(ncolx, ncoly, ncolz))
      Allocate(edenss(ncolx, ncoly, ncolz))
      Allocate(edenst(ncolx, ncoly, ncolz))
      Allocate(wctot(ncolx, ncoly, ncolz))
      Allocate(upot(ncolx, ncoly, ncolz, 2))
      Allocate(bmass(ncolx, ncoly, ncolz, 2))
      Allocate(dxiq(ncolx, ncoly, ncolz, 2))
      Allocate(rho(ncolx, ncoly, ncolz, 2))
      Allocate(tau(ncolx, ncoly, ncolz, 2))
      Allocate(worka(ncolx, ncoly, ncolz, 2))
      Allocate(rlam(ncolx, ncoly, ncolz, 2))
      Allocate(xiq(ncolx, ncoly, ncolz, 3, 2))
      Allocate(hsigma(ncolx, ncoly, ncolz, 3, 2))
      Allocate(cq(ncolx, ncoly, ncolz, 3, 2))
      Allocate(tmpvec(ncolx, ncoly, ncolz, 3, 2))
      Allocate(dbmass(ncolx, ncoly, ncolz, 3, 2))
      Allocate(currnt(ncolx, ncoly, ncolz, 3, 2))
      Allocate(sodens(ncolx, ncoly, ncolz, 3, 2))
      Allocate(spinden(ncolx, ncoly, ncolz, 3, 2))
      Allocate(kinvden(ncolx, ncoly, ncolz, 3, 2))
      Allocate(spincur(ncolx, ncoly, ncolz, 3, 3, 2))
      Allocate(dcq(ncolx, ncoly, ncolz, 3, 3, 2))
      Allocate(delxj(ncolx, ncoly, ncolz, 3, 3, 2))
      Allocate(bmunu(ncolx, ncoly, ncolz, 3, 3, 2))
      Allocate(dbmu(ncolx, ncoly, ncolz, 3, 2))
      npsif = 0
      npminf = 0
!-----------------------------------------------------------------------
!        determine/set number of processors on each node for openmp
!-----------------------------------------------------------------------
      Call OMP_SET_DYNAMIC(.False.)
      nprocs = OMP_GET_NUM_PROCS()
      Call OMP_SET_NUM_THREADS(nprocs)
      maxthreads = OMP_GET_MAX_THREADS()
      If(nprocs /= maxthreads) Then
         Write(8, '(//,a,//)') "WARNING: Number of processors different than number of threads"
      End If
!-----------------------------------------------------------------------
!        main static loop begins here.
!        several subsequent runs are possible.
!
!        allocate storage for wavefunctions and quantum numbers etc.
!        deallocate in case we are doing more than one nucleus at a time
!-----------------------------------------------------------------------
      If(imode == 0) Then
!         nfftw = fftw_init_threads()
!         call fftw_plan_with_nthreads(maxthreads)
!
         Do inuc = 1, nof
!
            If(allocated(psi)) deallocate(psi)
            If(allocated(nshell)) deallocate(nshell)
            If(allocated(listis)) deallocate(listis)
            If(allocated(pswk1)) deallocate(pswk1)
            If(allocated(pswk2)) deallocate(pswk2)
            If(allocated(pswk3)) deallocate(pswk3)
            If(allocated(pswk4)) deallocate(pswk4)
            If(allocated(tpsi)) deallocate(tpsi)
            If(allocated(vocc)) deallocate(vocc)
            If(allocated(alz)) deallocate(alz)
            If(allocated(asz)) deallocate(asz)
            If(allocated(ajz)) deallocate(ajz)
            If(allocated(spar)) deallocate(spar)
            If(allocated(spenrg)) deallocate(spenrg)
            If(allocated(spkine)) deallocate(spkine)
            If(allocated(spnorm)) deallocate(spnorm)
            If(allocated(speflu)) deallocate(speflu)
            If(allocated(q)) deallocate(q)
            If(allocated(rho2)) deallocate(rho2)
!
            nprot = Nint(fcharg(inuc))
            nneut = Nint(fmass(inuc)) - nprot
            npsi(1) = nneut
            npsi(2) = nprot
            If(nexadd == 1) Then
               npsi(1) = npsi(1) + nextra_n(inuc)
               npsi(2) = npsi(2) + nextra_p(inuc)
            End If
            ipairn = 0
            If(ipairf(1, inuc) /= 0) Then
               npsi(1) = NINT(nneut+1.65*FLOAT(nneut)**0.666667D0)!+nextra_n(inuc) !play with adding more states to pairing
               npsi(1) = npsi(1) - mod(npsi(1),2)
               ipairn = ipairf(1, inuc)
            End If
            ipairp = 0
            If(ipairf(2, inuc) /= 0) Then
               npsi(2) = NINT(nprot+1.65*FLOAT(nprot)**0.666667D0)!+nextra_p(inuc) !play with adding more states to pairing
               npsi(2) = npsi(2) - mod(npsi(2),2)
               ipairp = ipairf(2, inuc)
            End If
!-----------------------------------------------------------------------
!      only need 1/2 # of states if time reversal is present
!-----------------------------------------------------------------------
            If(itimrev /= 0) Then
               npsi(1) = npsi(1) / 2
               npsi(2) = npsi(2) / 2
            End If
!
            npsitot(1) = npsi(1)
            npsitot(2) = npsi(2)
            npmin(1) = 1
            npmin(2) = 1
!
            lpsi = Max(npsitot(1), npsitot(2))
            Allocate(psi(ncolx, ncoly, ncolz, 2, lpsi, 2))
            psi = (0.0_wp, 0.0_wp)
            Allocate(nshell(3, lpsi, 2), listis(lpsi, 2))
            Allocate(pswk1(ncolx, ncoly, ncolz, 2))
            Allocate(pswk2(ncolx, ncoly, ncolz, 2))
            Allocate(pswk3(ncolx, ncoly, ncolz, 2))
            Allocate(pswk4(ncolx, ncoly, ncolz, 2))
            Allocate(tpsi(ncolx, ncoly, ncolz, 2))
            Allocate(vocc(lpsi, 2), alz(lpsi, 2), asz(lpsi, 2), ajz(lpsi, 2), spar(lpsi, 2))
            Allocate(spenrg(lpsi, 2), spkine(lpsi, 2), spnorm(lpsi, 2), speflu(lpsi, 2))
            spenrg = 0.0_wp
            spkine = 0.0_wp
            spnorm = 0.0_wp
            speflu = 0.0_wp
            vocc = 0.0_wp
            alz = 0.0_wp
            asz = 0.0_wp
            ajz = 0.0_wp
            spar = 0.0_wp
!-----------------------------------------------------------------------
!        following lines include the cm fixing term for binding energy
!        calculations which are to be compared by others.
!        Dynamic case and the static runs for tdhf should not have this term.
!        calculations show that the addition of this term does not
!        change the wavefunctions in any appreciable way.
!-----------------------------------------------------------------------
            amass = nneut + nprot
            zmass = nprot
!
            If(ifixcm == 1) Then
               cmfact = (amass-1.0_wp) / amass
               h2m(1) = h2m(1) * cmfact
               h2m(2) = h2m(2) * cmfact
            End If
!
            xcm(1) = centf(1, inuc)
            ycm(1) = centf(2, inuc)
            zcm(1) = centf(3, inuc)
            xcm(2) = xcm(1)
            ycm(2) = ycm(1)
            zcm(2) = zcm(1)
!
            atheta = euler_alpha(inuc)
            btheta = euler_beta(inuc)
            gtheta = euler_gamma(inuc)
            radinx = radinf(1, inuc)
            radiny = radinf(2, inuc)
            radinz = radinf(3, inuc)
!
            iperc = 0 ! generally always non-periodic for static nuclei
            If(iperc == 0) Then
               Allocate(q(2*ncolx, 2*ncoly, 2*ncolz), rho2(2*ncolx, 2*ncoly, 2*ncolz))
               nx2 = 2 * ncolx
               ny2 = 2 * ncoly
               nz2 = 2 * ncolz
               Write(8, '(/,a,/)') " Poisson for isolated charge distribution"
            Else If(iperc == 1) Then
               Allocate(q(ncolx, ncoly, ncolz), rho2(ncolx, ncoly, ncolz))
               nx2 = ncolx
               ny2 = ncoly
               nz2 = ncolz
               Write(8, '(/,a,/)') " Poisson for periodic charge distribution"
            Else
               Write(8, '(//,a,//)') "WRONG VALUE OF IPERC ENTERED"
               Stop
            End If
!
            wcoul = 0.0_wp
            If(icoul /= 0) Then
               Call coulinit(coulplan1, coulplan2, nx2, ny2, nz2, wx(1), wy(1), wz(1), q, iperc)
            End If
!-----------------------------------------------------------------------
!        setup damping matrices in x,y,z directions separately.
!-----------------------------------------------------------------------
            Call setdmc(ncolx, ncoly, ncolz, wsk1x, wsk1y, wsk1z, der2x, der2y, der2z, wsk2x, wsk2y, wsk2z, h2m, e0dmp)
!
            Call setpsi(lpsi, ncolx, ncoly, ncolz, npsitot, npmin, nshell, radinx, radiny, radinz, wx, wy, wz, xclx, xcly, xclz, &
           &             psi, tpsi, xcm, ycm, zcm, atheta, btheta, gtheta, nneut, nprot, vocc, listis, itimrev, inuc, &
           &             nsplx, nsply, nsplz, nnotx, nnoty, nnotz, nordx, nordy, nordz, xkvx, xkvy, xkvz, 19, .FALSE.)
!---------------------------------------------------------------------------
!        start the damped static calculations.
!        set irwgs to 0 to since this is the static only calculation branch.
!        the output file t0gs.restart will still be created but not used.
!---------------------------------------------------------------------------
            irwgs0 = 0
            Call shf3d(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, itrbx, mprint, iprint, serr, h2m, t0, t1, t2, t3, &
           & t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
           & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, dbmass, dxiq, der1x, der1y, der1z, der2x, der2y, &
           & der2z, spenrg, speflu, spkine, vocc, spnorm, spar, ajz, psi, pswk1, pswk2, pswk3, pswk4, wsk2x, wsk2y, wsk2z, wx, &
           & wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, x0dmp, e0dmp, icoul, irest, mrests, mplots, tpsi, alz, asz, itimrev, nneut, &
           & nprot, ipairn, ipairp, fermin, fermip, epairn, epairp, worka, rlam, wxyz, c0, d0, iconstr, q2in, itheta, imode, irwgs0, &
           & iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, v0prot, rho0pr, nexadd, nexiter, etheta, ihdiag, .FALSE.)
!-----------------------------------------------------------------------
!        If added remove the extra non-pairing states
!-----------------------------------------------------------------------
            If(nexadd == 1) Then
               If(itimrev == 0) Then
                  Do iq = 1, 2
                     nst1 = 0
                     Do nst = 1, npsi(iq)
                        If(Nint(vocc(nst,iq)) == 1) Then
                           nst1 = nst1 + 1
                           psi(:,:,:,:,nst1,iq) = psi(:,:,:,:,nst,iq)
                        End If
                     End Do
                  End Do
                  npsi(1) = npsi(1) - nextra_n(inuc)
                  npsi(2) = npsi(2) - nextra_p(inuc)
                  npsif(1, inuc) = npsi1 - nextra_n(inuc)
                  npsif(2, inuc) = npsi2 - nextra_p(inuc)
                  npsi1 = npsi1 - nextra_n(inuc)
                  npsi2 = npsi2 - nextra_p(inuc)
                  psi(:,:,:,:,npsi(1)+1:npsi(1)+nextra_n(inuc),1) = (0.0_wp,0.0_wp)
                  psi(:,:,:,:,npsi(2)+1:npsi(2)+nextra_p(inuc),2) = (0.0_wp,0.0_wp)
               Else If(itimrev /= 0) Then
                  Do iq = 1, 2
                     nst1 = 0
                     Do nst = 1, npsi(iq)
                        If(Nint(vocc(nst,iq)) == 2) Then
                           nst1 = nst1 + 1
                           psi(:,:,:,:,nst1,iq) = psi(:,:,:,:,nst,iq)
                        End If
                     End Do
                  End Do
                  n1 = nextra_n(inuc)/2
                  n2 = nextra_p(inuc)/2
                  npsi(1) = npsi(1) - n1
                  npsi(2) = npsi(2) - n2
                  npsif(1, inuc) = npsi1 - n1
                  npsif(2, inuc) = npsi2 - n2
                  npsif(1, inuc) = 2*npsif(1, inuc)
                  npsif(2, inuc) = 2*npsif(2, inuc)
                  npsi1 = npsi1 - n1
                  npsi2 = npsi2 - n2
                  psi(:,:,:,:,npsi(1)+1:npsi(1)+n1,1) = (0.0_wp,0.0_wp)
                  psi(:,:,:,:,npsi(2)+1:npsi(2)+n2,2) = (0.0_wp,0.0_wp)
               End If
            End If
!
            npsitot(1) = npsi(1)
            npsitot(2) = npsi(2)
            Call writebspl(psi, binvpx, binvpy, binvpz, ncolx, ncoly, ncolz, ncolx, ncoly, ncolz, &
            & lpsi, npsi, npmin, nordx, nordy, nordz, etheta, xcm, ycm, zcm, xclx, xcly, xclz, vocc)
         End Do
!
         Write(8, '(2/,A,2/)') ' normal termination of shf3d'
!
         Close(Unit=3, Status='keep')
         Close(Unit=5, Status='keep')
         Close(Unit=8, Status='keep')
         Close(Unit=11, Status='keep')
         Close(Unit=12, Status='keep')
         Close(Unit=13, Status='keep')
         Close(Unit=14, Status='keep')
         Close(Unit=15, Status='keep')
!-----------------------------------------------------------------------
!        dynamic branch
!-----------------------------------------------------------------------
      Else If(imode == 1 .And. irest == 0) Then
!         nfftw = fftw_init_threads()
!         call fftw_plan_with_nthreads(maxthreads)
!-----------------------------------------------------------------------
!        First establish the total number of wavefunctions and allocate.
!        Time-dependent case does not have time-reversal invariance in
!        general but the static initial states can (doubled later by
!        generating the time-reversed states via the T-R operator)
!        We also consider possibility of pairing - double states
!-----------------------------------------------------------------------
         npsi(1) = 0
         npsi(2) = 0
!
         Do inuc = 1, nof
            nprot = Nint(fcharg(inuc))
            nneut = Nint(fmass(inuc)) - nprot
            npsi(1) = npsi(1) + nneut
            npsi(2) = npsi(2) + nprot
            If(nexadd == 1) Then
               npsi(1) = npsi(1) + nextra_n(inuc)
               npsi(2) = npsi(2) + nextra_p(inuc)
            End If
            If(ipairf(1, inuc) /= 0) Then
               npsi(1) = npsi(1) + NINT(1.65*FLOAT(nneut)**0.666667D0)
               npsi(1) = npsi(1) - mod(npsi(1),2)
            End If
            If(ipairf(2, inuc) /= 0) Then
               npsi(2) = npsi(2) + NINT(1.65*FLOAT(nprot)**0.666667D0)
               npsi(1) = npsi(1) - mod(npsi(1),2)
            End If
         End Do
!
         npsitot(1) = npsi(1)
         npsitot(2) = npsi(2)
!
         lpsi = Max(npsitot(1), npsitot(2))
!
         Allocate(psi(ncolx, ncoly, ncolz, 2, lpsi, 2))
         psi = (0.0_wp, 0.0_wp)
         Allocate(nshell(3, lpsi, 2), listis(lpsi, 2))
         Allocate(pswk1(ncolx, ncoly, ncolz, 2))
         Allocate(pswk2(ncolx, ncoly, ncolz, 2))
         Allocate(pswk3(ncolx, ncoly, ncolz, 2))
         Allocate(pswk4(ncolx, ncoly, ncolz, 2))
         Allocate(tpsi(ncolx, ncoly, ncolz, 2))
         Allocate(vocc(lpsi, 2), alz(lpsi, 2), asz(lpsi, 2), ajz(lpsi, 2), spar(lpsi, 2))
         Allocate(vocct(lpsi, 2), alzt(lpsi, 2), aszt(lpsi, 2), ajzt(lpsi, 2), spart(lpsi, 2))
         Allocate(spenrg(lpsi, 2), spkine(lpsi, 2), spnorm(lpsi, 2), speflu(lpsi, 2))
         Allocate(spenrgt(lpsi, 2), spkinet(lpsi, 2), spnormt(lpsi, 2), speflut(lpsi, 2))
         spenrg = 0.0_wp
         spkine = 0.0_wp
         spnorm = 0.0_wp
         speflu = 0.0_wp
         vocc = 0.0_wp
         alz = 0.0_wp
         asz = 0.0_wp
         ajz = 0.0_wp
         spar = 0.0_wp
!-----------------------------------------------------------------------
!        two fragments with boost adjustment
!        calculate the centers of the two nuclei
!-----------------------------------------------------------------------
         If(nof == 2) Then
            amc2 = hbc ** 2 / (h2m(1)+h2m(2))
            xmu = fmass(1) * fmass(2) / (fmass(1)+fmass(2)) * amc2
            ratio = xmu / hbc ** 2
            akv = Sqrt(2.0_wp*ratio*ecm)
            xli = xb * akv
            tdot = xli / (xmu/hbc*rsep**2)
            rstr = 1.4_wp * (fmass(1)**(1.0_wp/3.0_wp)+fmass(2)**(1.0_wp/3.0_wp))
            xlg = rstr * akv
            epsi = Sqrt(1.0_wp+2.0_wp*ecm*xli**2/(ratio*(fcharg(1)*fcharg(2)*e2)**2))
            rmin = fcharg(1) * fcharg(2) * e2 / (2.0_wp*ecm+1.0e-25_wp) * (1.0_wp+epsi)
            If(rsep <= rmin .And. ecm /= 0.0_wp) Then
               Write(8, '(2/,A)') "Relative nuclear distance too small for Coulomb initialization"
               Write(8, '(A,F10.4,A,2/)') "Choose a distance > ", rmin, " fm"
               Stop
            End If
            ttt = xli ** 2 / (ratio*fcharg(1)*fcharg(2)*e2*rsep)
            tetr = Acos(1.0_wp/epsi) - Acos(1.0_wp/epsi*(1.0_wp+ttt))
            If(xb < 1.0e-2_wp) tetr = 0.0_wp
         End If
         If(nof == 2 .And. ifixb == 1) Then
            centf = 0.0_wp
            cofl = fmass(2) / (fmass(1)+fmass(2))
            cofr = fmass(1) / (fmass(1)+fmass(2))
            centf(1, 1) = - cofl * rsep * Cos(tetr)
            centf(3, 1) = - cofl * rsep * Sin(tetr)
            centf(1, 2) = cofr * rsep * Cos(tetr)
            centf(3, 2) = cofr * rsep * Sin(tetr)
         End If
!
         npmin(1) = 1
         npmin(2) = 1
         npsi(1) = 0
         npsi(2) = 0
         wctot = 0.0_wp
!
         Do inuc = 1, nof
!-----------------------------------------------------------------------
!        first establish the number of wavefunctions for inuc and
!        define the proper indexing of wavefunctions
!-----------------------------------------------------------------------
            nprot = Nint(fcharg(inuc))
            nneut = Nint(fmass(inuc)) - nprot
            npsi1 = nneut
            npsi2 = nprot
            If(nexadd == 1) Then
               npsi1 = npsi1 + nextra_n(inuc)
               npsi2 = npsi2 + nextra_p(inuc)
            End If
            ipairn = 0
            If(ipairf(1, inuc) /= 0) Then
               npsi1 = npsi1 + NINT(1.65*FLOAT(nneut)**0.666667D0)
               npsi1 = npsi1 - mod(npsi1,2)
               ipairn = ipairf(1, inuc)
            End If
            ipairp = 0
            If(ipairf(2, inuc) /= 0) Then
               npsi2 = npsi2 + NINT(1.65*FLOAT(nprot)**0.666667D0)
               npsi2 = npsi2 - mod(npsi2,2)
               ipairp = ipairf(2, inuc)
            End If
            npsif(1, inuc) = npsi1
            npsif(2, inuc) = npsi2
!-----------------------------------------------------------------------
!        We first do a T-R invariant calculation and then create the
!        remaining states via T-R operator and divide the occupation
!        weights by two. This works if we have L-N as well.
!        (except for odd A nuclei).
!-----------------------------------------------------------------------
!      only need 1/2 # of states initially if time reversal is present
!-----------------------------------------------------------------------
            If(itimrevs /= 0) Then
               npsi1 = npsi1 / 2
               npsi2 = npsi2 / 2
            End If
            iunit = 18 + inuc
            If(npsif(1, inuc)==npsi(1) .AND. npsif(2, inuc)==npsi(2)) iunit = 19
            npsi(1) = npsi(1) + npsi1
            npsi(2) = npsi(2) + npsi2
!
            npsitot(1) = npsi(1)
            npsitot(2) = npsi(2)
!
            xcm(1) = centf(1, inuc)
            ycm(1) = centf(2, inuc)
            zcm(1) = centf(3, inuc)
            xcm(2) = xcm(1)
            ycm(2) = ycm(1)
            zcm(2) = zcm(1)
!
            atheta = euler_alpha(inuc)
            btheta = euler_beta(inuc)
            gtheta = euler_gamma(inuc)
            radinx = radinf(1, inuc)
            radiny = radinf(2, inuc)
            radinz = radinf(3, inuc)
!
            If(allocated(q)) deallocate(q)
            If(allocated(rho2)) deallocate(rho2)
            iperc = 0
            If(iperc == 0) Then
               Allocate(q(2*ncolx, 2*ncoly, 2*ncolz), rho2(2*ncolx, 2*ncoly, 2*ncolz))
               nx2 = 2 * ncolx
               ny2 = 2 * ncoly
               nz2 = 2 * ncolz
               Write(8, '(/,a,/)') " Poisson for isolated charge distribution"
            Else If(iperc == 1) Then
               Allocate(q(ncolx, ncoly, ncolz), rho2(ncolx, ncoly, ncolz))
               nx2 = ncolx
               ny2 = ncoly
               nz2 = ncolz
               Write(8, '(/,a,/)') " Poisson for periodic charge distribution"
            Else
               Write(8, '(//,a,//)') "WRONG VALUE OF IPERC ENTERED"
               Stop
            End If
!
            wcoul = 0.0_wp
            If(icoul /= 0) Then
               Call coulinit(coulplan1, coulplan2, nx2, ny2, nz2, wx(1), wy(1), wz(1), q, iperc)
            End If
!-----------------------------------------------------------------------
!        setup damping matrices in x,y,z directions separately.
!-----------------------------------------------------------------------
            Call setdmc(ncolx, ncoly, ncolz, wsk1x, wsk1y, wsk1z, der2x, der2y, der2z, wsk2x, wsk2y, wsk2z, h2m, e0dmp)
!
            Call setpsi(lpsi, ncolx, ncoly, ncolz, npsitot, npmin, nshell, radinx, radiny, radinz, wx, wy, wz, xclx, xcly, xclz, &
           &            psi, tpsi, xcm, ycm, zcm, atheta, btheta, gtheta, nneut, nprot, vocc, listis, itimrevs, inuc, &
           &            nsplx, nsply, nsplz, nnotx, nnoty, nnotz, nordx, nordy, nordz, xkvx, xkvy, xkvz, iunit, isrestart(inuc))
!---------------------------------------------------------------------------
!        start the static calculations for for tdhf setup.
!        do not do constraint and do not plot static densities in this case
!---------------------------------------------------------------------------
            iconstr0 = 0
            imode0 = 0
            Call shf3d(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, itrbx, mprint, iprint, serr, h2m, t0, t1, t2, t3, &
           & t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
           & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, dbmass, dxiq, der1x, der1y, der1z, der2x, der2y, &
           & der2z, spenrg, speflu, spkine, vocc, spnorm, spar, ajz, psi, pswk1, pswk2, pswk3, pswk4, wsk2x, wsk2y, wsk2z, wx, &
           & wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, x0dmp, e0dmp, icoul, irest, mrests, mplots, tpsi, alz, asz, itimrevs, nneut, &
           & nprot, ipairn, ipairp, fermin, fermip, epairn, epairp, worka, rlam, wxyz, c0, d0, iconstr0, q2in, itheta, imode0, irwgs, &
           &  iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, v0prot, rho0pr, nexadd, nexiter, etheta, ihdiag, isrestart(inuc))
! calculate the coulomb energy for each nucleus to be used in boost calculation below
            ecx = 0.0_wp
            ecoul(inuc) = 0.0_wp
            If(icoul /= 0) Then
               Do iz = 1, ncolz
                  Do iy = 1, ncoly
                     Do ix = 1, ncolx
                        ecx = ecx + 0.5_wp*wx(ix)*wy(iy)*wz(iz)*wcoul(ix, iy, iz)*rho(ix, iy, iz, 2)
                     End Do
                  End Do
               End Do
               ecoul(inuc) = ecx
            End If
! save the total coulomb potential
            wctot = wctot + wcoul
!-----------------------------------------------------------------------
!        If added remove the extra non-pairing states
!-----------------------------------------------------------------------
            If(nexadd == 1) Then
               If(itimrevs == 0) Then
                  Do iq = 1, 2
                     nst1 = 0
                     Do nst = 1, npsi(iq)
                        If(nint(vocc(nst,iq)) == 1) Then
                           nst1 = nst1 + 1
                           psi(:,:,:,:,nst1,iq) = psi(:,:,:,:,nst,iq)
                           spenrgt(nst1, iq) = spenrg(nst, iq)
                           spkinet(nst1, iq) = spkine(nst, iq)
                           speflut(nst1, iq) = speflu(nst, iq)
                           spnormt(nst1, iq) = spnorm(nst, iq)
                           alzt(nst1, iq) = alz(nst, iq)
                           aszt(nst1, iq) = asz(nst, iq)
                           ajzt(nst1, iq) = ajz(nst, iq)
                           spart(nst1, iq) = spar(nst, iq)
                           vocct(nst1, iq) = vocc(nst, iq)
                        End If
                     End Do
                  End Do
                  Do iq = 1, 2
                     Do nst = npmin(iq), npsi(iq)
                        spenrg(nst, iq) = spenrgt(nst, iq)
                        spkine(nst, iq) = spkinet(nst, iq)
                        speflu(nst, iq) = speflut(nst, iq)
                        spnorm(nst, iq) = spnormt(nst, iq)
                        alz(nst, iq) = alzt(nst, iq)
                        asz(nst, iq) = aszt(nst, iq)
                        ajz(nst, iq) = ajzt(nst, iq)
                        spar(nst, iq) = spart(nst, iq)
                        vocc(nst, iq) = vocct(nst, iq)
                     End Do
                  End Do
                  npsi(1) = npsi(1) - nextra_n(inuc)
                  npsi(2) = npsi(2) - nextra_p(inuc)
                  npsif(1, inuc) = npsi1 - nextra_n(inuc)
                  npsif(2, inuc) = npsi2 - nextra_p(inuc)
                  npsi1 = npsi1 - nextra_n(inuc)
                  npsi2 = npsi2 - nextra_p(inuc)
                  psi(:,:,:,:,npsi(1)+1:npsi(1)+nextra_n(inuc),1) = (0.0_wp,0.0_wp)
                  psi(:,:,:,:,npsi(2)+1:npsi(2)+nextra_p(inuc),2) = (0.0_wp,0.0_wp)
               Else If(itimrevs /= 0) Then
                  Do iq = 1, 2
                     nst1 = 0
                     Do nst = 1, npsi(iq)
                        If(nint(vocc(nst,iq)) == 2) Then
                           nst1 = nst1 + 1
                           psi(:,:,:,:,nst1,iq) = psi(:,:,:,:,nst,iq)
                           spenrgt(nst1, iq) = spenrg(nst, iq)
                           spkinet(nst1, iq) = spkine(nst, iq)
                           speflut(nst1, iq) = speflu(nst, iq)
                           spnormt(nst1, iq) = spnorm(nst, iq)
                           alzt(nst1, iq) = alz(nst, iq)
                           aszt(nst1, iq) = asz(nst, iq)
                           ajzt(nst1, iq) = ajz(nst, iq)
                           spart(nst1, iq) = spar(nst, iq)
                           vocct(nst1, iq) = vocc(nst, iq)
                        End If
                     End Do
                  End Do
                  Do iq = 1, 2
                     Do nst = npmin(iq), npsi(iq)
                        spenrg(nst, iq) = spenrgt(nst, iq)
                        spkine(nst, iq) = spkinet(nst, iq)
                        speflu(nst, iq) = speflut(nst, iq)
                        spnorm(nst, iq) = spnormt(nst, iq)
                        alz(nst, iq) = alzt(nst, iq)
                        asz(nst, iq) = aszt(nst, iq)
                        ajz(nst, iq) = ajzt(nst, iq)
                        spar(nst, iq) = spart(nst, iq)
                        vocc(nst, iq) = vocct(nst, iq)
                     End Do
                  End Do
                  n1 = nextra_n(inuc)/2
                  n2 = nextra_p(inuc)/2
                  npsi(1) = npsi(1) - n1
                  npsi(2) = npsi(2) - n2
                  npsif(1, inuc) = npsi1 - n1
                  npsif(2, inuc) = npsi2 - n2
                  npsif(1, inuc) = 2*npsif(1, inuc)
                  npsif(2, inuc) = 2*npsif(2, inuc)
                  npsi1 = npsi1 - n1
                  npsi2 = npsi2 - n2
                  psi(:,:,:,:,npsi(1)+1:npsi(1)+n1,1) = (0.0_wp,0.0_wp)
                  psi(:,:,:,:,npsi(2)+1:npsi(2)+n2,2) = (0.0_wp,0.0_wp)
               End If
            End If
!
            npsitot(1) = npsi(1)
            npsitot(2) = npsi(2)
!-----------------------------------------------------------------------
!        generate the other states via T-R operator
!-----------------------------------------------------------------------
            If(itimrevs /= 0) Then
               Allocate(psitmp(ncolx, ncoly, ncolz, 2, lpsi, 2))
               psitmp = (0.0_wp, 0.0_wp)
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     nst1 = 2 * (nst-1) + 1 - (npmin(iq)-1)
                     nst2 = 2 * nst - (npmin(iq)-1)
                     psitmp(:ncolx, :ncoly, :ncolz, :2, nst1, iq) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                     Call trpsi(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, iq), psitmp(1, 1, 1, 1, nst2, iq))
                     spenrgt(nst1, iq) = spenrg(nst, iq)
                     spenrgt(nst2, iq) = spenrgt(nst1, iq)
                     spkinet(nst1, iq) = spkine(nst, iq)
                     spkinet(nst2, iq) = spkinet(nst1, iq)
                     speflut(nst1, iq) = speflu(nst, iq)
                     speflut(nst2, iq) = speflut(nst1, iq)
                     spnormt(nst1, iq) = spnorm(nst, iq)
                     spnormt(nst2, iq) = spnormt(nst1, iq)
                     alzt(nst1, iq) = alz(nst, iq)
                     alzt(nst2, iq) = alzt(nst1, iq)
                     aszt(nst1, iq) = asz(nst, iq)
                     aszt(nst2, iq) = aszt(nst1, iq)
                     ajzt(nst1, iq) = ajz(nst, iq)
                     ajzt(nst2, iq) = ajzt(nst1, iq)
                     spart(nst1, iq) = spar(nst, iq)
                     spart(nst2, iq) = spart(nst1, iq)
                     vocct(nst1, iq) = vocc(nst, iq) / 2.0_wp
                     vocct(nst2, iq) = vocct(nst1, iq)
                  End Do
               End Do
!
               npsi(1) = npsi(1) + npsi1
               npsi(2) = npsi(2) + npsi2
               npsitot(1) = npsi(1)
               npsitot(2) = npsi(2)
!
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     psi(:ncolx, :ncoly, :ncolz, :2, nst, iq) = psitmp(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                     spenrg(nst, iq) = spenrgt(nst, iq)
                     spkine(nst, iq) = spkinet(nst, iq)
                     speflu(nst, iq) = speflut(nst, iq)
                     spnorm(nst, iq) = spnormt(nst, iq)
                     alz(nst, iq) = alzt(nst, iq)
                     asz(nst, iq) = aszt(nst, iq)
                     ajz(nst, iq) = ajzt(nst, iq)
                     spar(nst, iq) = spart(nst, iq)
                     vocc(nst, iq) = vocct(nst, iq)
                  End Do
               End Do
               Deallocate(psitmp)
            End If
!
            npminf(1, inuc) = npmin(1)
            npminf(2, inuc) = npmin(2)
!
            If(itimrevs /= 0) Then
               npmin(1) = npmin(1) + 2 * npsi1
               npmin(2) = npmin(2) + 2 * npsi2
            Else
               npmin(1) = npmin(1) + npsi1
               npmin(2) = npmin(2) + npsi2
            End If
         End Do
         Allocate(psitmp(ncolx, ncoly, ncolz, 2, lpsi, 2))
         psitmp = (0.0_wp, 0.0_wp)
         psitmp = psi
!-----------------------------------------------------------------------
!        reset the wavefunction beginning counter
!-----------------------------------------------------------------------
         npmin(1) = 1
         npmin(2) = 1
         nneut = 0
         nprot = 0
         Do inuc = 1, nof
            nprot = nprot + Nint(fcharg(inuc))
            nneut = nneut + Nint(fmass(inuc)) - Nint(fcharg(inuc))
         End Do
!-----------------------------------------------------------------------
!        calculate the initial quantities again. this time using
!        the exact separation between the two centers.
!        use tevolv for one time step.
!        NOTE: for this case we turn off printing, constraint etc.
!-----------------------------------------------------------------------
         iconstr0 = 0
         ipairn0 = 0
         ipairp0 = 0
         niter = 0
         xdt = -1.0_wp
         Call tevolv(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, npsif, npminf, 1, 1, 0, 0, mxp, terr, xdt, h2m, &
         & t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, &
         & delxj, bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, edenst, dbmass, dxiq, der1x, &
         & der1y, der1z, der2x, der2y, der2z, spenrg, spkine, vocc, spnorm, spar, ajz, psitmp, pswk1, pswk2, pswk3, pswk4, &
         & wx, wy, wz, wxyz, xclx, xcly, xclz, xcm, ycm, zcm, centf, roft, rdot, tdot, vx, vz, mnof, nof, itheta, itimrev, &
         & itimrevs, icoul, ecoul, irest, mrestt, worka, rlam, mtrbx, serr, speflu, wsk2x, wsk2y, wsk2z, x0dmp, e0dmp, tpsi, alz, &
         & asz, nneut, nprot, ipairn0, ipairp0, c0, d0, iconstr0, mconstr, imode, ifixb, niter, nfixb, iodds, ihdiag, rsep, xb, &
         & coulplan1, coulplan2, rho2, q, iperc, directname, v0neut, v0prot, rho0pr, fermin, fermip, epairn, epairp)
!
         If(nof == 2) Then
            delx = centf(1, 2) - centf(1, 1)
            dely = centf(2, 2) - centf(2, 1)
            delz = centf(3, 2) - centf(3, 1)
            roft0 = Sqrt(delx**2+dely**2+delz**2)
            roft = roft0
            psitmp = psi
            wcoul = wctot
!-----------------------------------------------------------------------
!        get the exact coulomb potential.
!-----------------------------------------------------------------------
            rho = 0.0_wp
            tau = 0.0_wp
            currnt = 0.0_wp
            sodens = 0.0_wp
            spinden = 0.0_wp
            kinvden = 0.0_wp
            spincur = 0.0_wp
            rhos = 0.0_wp
            taus = 0.0_wp
            drhos = 0.0_wp
            currnts = 0.0_wp
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                 Call densit(lpsi, ncolx, ncoly, ncolz, nst, iq, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, psi(1, 1, 1, 1, nst, iq), pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
               End Do
            End Do
            Deallocate( rhos, taus, drhos, currnts )
            If(icoul /= 0) Then
               Call pois(iperc, coulplan1, coulplan2, ncolx, ncoly, ncolz, nx2, ny2, nz2, wx(1), wy(1), wz(1),&
                         rho(1, 1, 1, 2), rho2, q, wcoul)
! calculate the coulomb energy for the combined system to be used in boost calculation below
               ecoultot = 0.0_wp
                Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        ecoultot = ecoultot + 0.5_wp * wx(ix) * wy(iy) * wz(iz) * wcoul(ix, iy, iz) * rho(ix, iy, iz, 2)
                     End Do
                  End Do
               End Do
            End If
!
            If(roft <= rmin .And. ecm /= 0.0_wp) Then
               Write(8, '(2/,A)') "Dynamic nuclear separation too small for Coulomb initialization"
               Write(8, '(A,F10.4,A,2/)') "Choose a distance > ", rmin, " fm"
               Stop
            End If
            xcoul = fcharg(1) * fcharg(2) * e2 / roft
            If(icoul /= 0) Then
               xcoul = ecoultot - ecoul(1) - ecoul(2)
            End If
            xcent = xli ** 2 / (2.0*ratio*roft**2)
            ttt = xli ** 2 / (ratio*fcharg(1)*fcharg(2)*e2*roft)
            epsi = Sqrt(1.0_wp+2.0_wp*ecm*xli**2/(ratio*(fcharg(1)*fcharg(2)*e2)**2))
            tetr = Acos(1.0_wp/epsi) - Acos(1.0_wp/epsi*(1.0_wp+ttt))
            tetc = pi - tetr
            tetf = 2.0_wp * Acos(1.0_wp/epsi) - tetr
            tets = tetc - tetf
            tdot = xli / (xmu/hbc*roft**2)
            ecm0 = ecm - xcoul - xcent
            If(ecm0 > 0.0_wp) Then
               rdot = - Sqrt(2.0_wp*ecm0/xmu)
            Else
               rdot = 0.0_wp
            End If
            vx = rdot * Cos(tetr) - roft * Sin(tetr) * tdot
            vz = rdot * Sin(tetr) + roft * Cos(tetr) * tdot
            tcmx = 0.5_wp * xmu * vx ** 2
            tcmz = 0.5_wp * xmu * vz ** 2
            tketot = tcmx + tcmz
            Write(8, '(2/,A,/,6(A,F12.4,A,F12.4,A,F12.4,A,F12.4,/),2/)') ' initial coulomb kinematics:', '  red. mass=', xmu / &
           & amc2, '     l/hbar=', xli, ' lgraz/hbar=', xlg, '      b(fm)=', xb, '   ecm(mev)=', ecm, '   roft(fm)=', roft, '    &
           & rdot/c=', rdot, '  td/c(deg)=', tdot * 180.0_wp / pi, '  trke(mev)=', tketot, '  tetr(deg)=', tetr * 180.0_wp / pi, &
           & '  tetc(deg)=', tetc * 180.0_wp / pi, '  tets(deg)=', tets * 180.0_wp / pi, '    mass 1 =', fmass(1), '    mass 2 ='&
           & , fmass(2), '  charge 1 =', fcharg(1), '  charge 2 =', fcharg(2), '  center 1x=', centf(1, 1), '  center 1z=', &
           & centf(3, 1), '  center 2x=', centf(1, 2), '  center 2z=', centf(3, 2), ' vcoul(mev)=', xcoul, ' vcent(mev)=', xcent
         End If
!
         If(nof == 2 .And. ifixb == 1) Then
            cofl = fmass(2) / (fmass(1)+fmass(2))
            cofr = fmass(1) / (fmass(1)+fmass(2))
            vlx = - cofl * vx
            vlz = - cofl * vz
            vrx = cofr * vx
            vrz = cofr * vz
            tlx = 0.5_wp * fmass(1) * amc2 * vlx ** 2
            tlz = 0.5_wp * fmass(1) * amc2 * vlz ** 2
            trx = 0.5_wp * fmass(2) * amc2 * vrx ** 2
            trz = 0.5_wp * fmass(2) * amc2 * vrz ** 2
!
            If(vx /= 0.0_wp) Then
               boostf(1, 1) = tlx * (vlx/Abs(vlx))
               boostf(1, 2) = trx * (vrx/Abs(vrx))
            Else
               boostf(1, 1) = 0.0_wp
               boostf(1, 2) = 0.0_wp
            End If
            If(vz /= 0.0_wp .And. xb > 0.1_wp) Then
               boostf(3, 1) = tlz * (vlz/Abs(vlz))
               boostf(3, 2) = trz * (vrz/Abs(vrz))
            Else
               boostf(3, 1) = 0.0_wp
               boostf(3, 2) = 0.0_wp
            End If
!-----------------------------------------------------------------------
!        dynamic fitting of the c.m. energy begins here. this is only
!        neccesary if we want precise ecm value. since we correctly
!        calculate the initial coulomb energy the actual ecm should be
!        close to the input value. Otherwise make nfixb above 10.
!        if the parameter ifixb=0 then we use the boosts read from input.
!-----------------------------------------------------------------------
            nbloop = nfixb
            If(ecm == 0.0_wp) nbloop = 0
            Do niter = 1, nbloop
               roft = roft0
               psitmp = psi
               wcoul = wctot
!
               Do inuc = 1, nof
                  Call bustf(lpsi, ncolx, ncoly, ncolz, npsif, npminf, mnof, xclx, xcly, xclz, boostf, centf, fmass, h2m, psitmp, &
                 & inuc)
               End Do
!-----------------------------------------------------------------------
!        use tevolv for one time step.
!        NOTE: for this case we turn off printing, constraint etc.
!-----------------------------------------------------------------------
               iconstr0 = 0
               ipairn0 = 0
               ipairp0 = 0
               Call tevolv(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, npsif, npminf, 1, 1, 0, 0, mxp, terr, dt, h2m, &
              & t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, &
              & delxj, bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, edenst, dbmass, dxiq, der1x, &
              & der1y, der1z, der2x, der2y, der2z, spenrg, spkine, vocc, spnorm, spar, ajz, psitmp, pswk1, pswk2, pswk3, pswk4, &
              & wx, wy, wz, wxyz, xclx, xcly, xclz, xcm, ycm, zcm, centf, roft, rdot, tdot, vx, vz, mnof, nof, itheta, itimrev, &
              & itimrevs, icoul, ecoul, irest, mrestt, worka, rlam, mtrbx, serr, speflu, wsk2x, wsk2y, wsk2z, x0dmp, e0dmp, tpsi, &
              & alz, asz, nneut, nprot, ipairn0, ipairp0, c0, d0, iconstr0, mconstr, imode, ifixb, niter, nfixb, iodds, ihdiag, &
              & rsep, xb, coulplan1, coulplan2, rho2, q, iperc, directname, v0neut, v0prot, rho0pr, fermin, fermip, epairn, epairp)
! keep vz same
               vlx = - cofl * vx
               vlz = - cofl * vz
               vrx = cofr * vx
               vrz = cofr * vz
               tlix = 0.5_wp * fmass(1) * amc2 * vlx ** 2
               tliz = 0.5_wp * fmass(1) * amc2 * vlz ** 2
               trix = 0.5_wp * fmass(2) * amc2 * vrx ** 2
               triz = 0.5_wp * fmass(2) * amc2 * vrz ** 2
!
               epslx = tlx - tlix
               epslz = tlz - tliz
               epsrx = trx - trix
               epsrz = trz - triz
!
               boostf(1, 1) = boostf(1, 1) + epslx / 2.0_wp
               boostf(1, 2) = boostf(1, 2) - epsrx / 2.0_wp
               If(xb > 0.1_wp) Then
                  boostf(3, 1) = boostf(3, 1) + epslz / 2.0_wp
                  boostf(3, 2) = boostf(3, 2) - epsrz / 2.0_wp
               End If
!
            End Do
         End If
         Deallocate(psitmp)
!
         wcoul = wctot
         If(nof >= 1) Then
            Do inuc = 1, nof
               Call bustf(lpsi, ncolx, ncoly, ncolz, npsif, npminf, mnof, xclx, xcly, xclz, boostf, centf, fmass, h2m, psi, inuc)
            End Do
         End If
!
         ipairn0 = 0
         ipairp0 = 0
         roft = roft0
         niter = nfixb + 1
         Call tevolv(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, npsif, npminf, nt, mprint, mplott, iprint, mxp, terr, &
        & dt, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, &
        & spincur, delxj, bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, edenst, dbmass, dxiq, &
        & der1x, der1y, der1z, der2x, der2y, der2z, spenrg, spkine, vocc, spnorm, spar, ajz, psi, pswk1, pswk2, pswk3, pswk4, wx, &
        & wy, wz, wxyz, xclx, xcly, xclz, xcm, ycm, zcm, centf, roft, rdot, tdot, vx, vz, mnof, nof, itheta, itimrev, itimrevs, &
        & icoul, ecoul, irest, mrestt, worka, rlam, mtrbx, derr, speflu, wsk2x, wsk2y, wsk2z, x0dmp, e0dmp, tpsi, alz, asz, nneut,&
        & nprot, ipairn0, ipairp0, c0, d0, iconstr, mconstr, imode, ifixb, niter, nfixb, iodds, ihdiag, rsep, xb, coulplan1, &
        & coulplan2, rho2, q, iperc, directname, v0neut, v0prot, rho0pr, fermin, fermip, epairn, epairp)
!
         Write(8, '(2/,A,2/)') ' normal termination of tdhf3d'
         Deallocate(q, rho2)
!
         Close(Unit=3, Status='keep')
         Close(Unit=5, Status='keep')
         Close(Unit=8, Status='keep')
         Close(Unit=11, Status='keep')
         Close(Unit=12, Status='keep')
         Close(Unit=13, Status='keep')
         Close(Unit=14, Status='keep')
         Close(Unit=15, Status='keep')
!
      Else If(imode == 1 .And. irest /= 0) Then
!-----------------------------------------------------------------------
!        first establish the total number of wavefunctions and allocate
!-----------------------------------------------------------------------
         npsi(1) = 0
         npsi(2) = 0
!
         Do inuc = 1, nof
            nprot = Nint(fcharg(inuc))
            nneut = Nint(fmass(inuc)) - nprot
            npsi(1) = npsi(1) + nneut
            npsi(2) = npsi(2) + nprot
            If(ipairf(1, inuc) /= 0) Then
               npsi(1) = npsi(1) + NINT(1.65*FLOAT(nneut)**0.666667D0)
               npsi(1) = npsi(1) - mod(npsi(1),2)
            End If
            If(ipairf(2, inuc) /= 0) Then
               npsi(2) = npsi(2) + NINT(1.65*FLOAT(nprot)**0.666667D0)
               npsi(2) = npsi(2) - mod(npsi(2),2)
            End If
         End Do
!-----------------------------------------------------------------------
!      only need 1/2 # of states if time reversal is present
!-----------------------------------------------------------------------
         If(itimrev /= 0) Then
            npsi(1) = npsi(1) / 2
            npsi(2) = npsi(2) / 2
         End If
!
         npsitot(1) = npsi(1)
         npsitot(2) = npsi(2)
         nneut = 0
         nprot = 0
         Do inuc = 1, nof
            nprot = nprot + Nint(fcharg(inuc))
            nneut = nneut + Nint(fmass(inuc)) - Nint(fcharg(inuc))
         End Do
!
         lpsi = Max(npsitot(1), npsitot(2))
!
         Allocate(psi(ncolx, ncoly, ncolz, 2, lpsi, 2))
         Allocate(nshell(3, lpsi, 2), listis(lpsi, 2))
         Allocate(pswk1(ncolx, ncoly, ncolz, 2))
         Allocate(pswk2(ncolx, ncoly, ncolz, 2))
         Allocate(pswk3(ncolx, ncoly, ncolz, 2))
         Allocate(pswk4(ncolx, ncoly, ncolz, 2))
         Allocate(tpsi(ncolx, ncoly, ncolz, 2))
         Allocate(vocc(lpsi, 2), alz(lpsi, 2), asz(lpsi, 2), ajz(lpsi, 2), spar(lpsi, 2))
         Allocate(spenrg(lpsi, 2), spkine(lpsi, 2), spnorm(lpsi, 2), speflu(lpsi, 2))
!
         iperc = 0
         If(iperc == 0) Then
            Allocate(q(2*ncolx, 2*ncoly, 2*ncolz), rho2(2*ncolx, 2*ncoly, 2*ncolz))
            nx2 = 2 * ncolx
            ny2 = 2 * ncoly
            nz2 = 2 * ncolz
            Write(8, '(/,a,/)') " Poisson for isolated charge distribution"
         Else If(iperc == 1) Then
            Allocate(q(ncolx, ncoly, ncolz), rho2(ncolx, ncoly, ncolz))
            nx2 = ncolx
            ny2 = ncoly
            nz2 = ncolz
            Write(8, '(/,a,/)') " Poisson for periodic charge distribution"
         Else
            Write(8, '(//,a,//)') "WRONG VALUE OF IPERC ENTERED"
            Stop
         End If
!
         wcoul = 0.0_wp
         If(icoul /= 0) Then
            Call coulinit(coulplan1, coulplan2, nx2, ny2, nz2, wx(1), wy(1), wz(1), q, iperc)
         End If
!-----------------------------------------------------------------------
!        calculate the damping matrices in case we do density constraint
!-----------------------------------------------------------------------
         Call setdmc(ncolx, ncoly, ncolz, wsk1x, wsk1y, wsk1z, der2x, der2y, der2z, wsk2x, wsk2y, wsk2z, h2m, e0dmp)
!
         ipairn0 = 0
         ipairp0 = 0
         niter = nfixb + 1
         Call tevolv(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, npsif, npminf, nt, mprint, mplott, iprint, mxp, terr, &
        & dt, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, &
        & spincur, delxj, bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, edenst, dbmass, dxiq, &
        & der1x, der1y, der1z, der2x, der2y, der2z, spenrg, spkine, vocc, spnorm, spar, ajz, psi, pswk1, pswk2, pswk3, pswk4, wx, &
        & wy, wz, wxyz, xclx, xcly, xclz, xcm, ycm, zcm, centf, roft, rdot, tdot, vx, vz, mnof, nof, itheta, itimrev, itimrevs, &
        & icoul, ecoul, irest, mrestt, worka, rlam, mtrbx, derr, speflu, wsk2x, wsk2y, wsk2z, x0dmp, e0dmp, tpsi, alz, asz, nneut,&
        & nprot, ipairn0, ipairp0, c0, d0, iconstr, mconstr, imode, ifixb, niter, nfixb, iodds, ihdiag, rsep, xb, coulplan1, &
        & coulplan2, rho2, q, iperc, directname, v0neut, v0prot, rho0pr, fermin, fermip, epairn, epairp)
!
         Write(8, '(2/,A,2/)') ' normal termination of tdhf3d'
!
         Close(Unit=3, Status='keep')
         Close(Unit=5, Status='keep')
         Close(Unit=8, Status='keep')
         Close(Unit=11, Status='keep')
         Close(Unit=12, Status='keep')
         Close(Unit=13, Status='keep')
         Close(Unit=14, Status='keep')
         Close(Unit=15, Status='keep')
         Close(Unit=19, Status='keep')
         Close(Unit=20, Status='keep')
      Else
         Write(8, '(/A/)') ' Wrong value of imode entered'
         Stop
      End If
!
End Program tdhf3d
