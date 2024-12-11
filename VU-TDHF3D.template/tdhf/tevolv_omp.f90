Subroutine tevolv(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, npsif, npminf, nt, mprint, mplot, iprint, mxp, terr, dt, &
& h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, &
& bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, edenst, dbmass, dxiq, der1x, der1y, der1z, der2x, &
& der2y, der2z, spenrg, spkine, vocc, spnorm, spar, ajz, psi, pswk1, pswk2, pswk3, pswk4, wx, wy, wz, wxyz, xclx, xcly, xclz, &
& xcm, ycm, zcm, centf, roft, rdot, tdot, vx, vz, mnof, nof, itheta, itimrev, itimrevs, icoul, ecoul, irest, mrest, worka, rlam, &
& itrbx, serr, speflu, cdmpx, cdmpy, cdmpz, x0dmp, e0dmp, tpsi, alz, asz, nneut, nprot, ipairn, ipairp, c0, d0, iconstr, mconstr, &
& imode, ifixb, niter, nfixb, iodds, ihdiag, rsep, xb, coulplan1, coulplan2, rho2, q, iperc, directname, v0neut, v0prot, rho0pr, &
& fermin, fermip, epairn, epairp)
!
!-----------------------------------------------------------------------
!     tevolv  =  routine for dynamic tdhf evolution
!                iterates all wave functions one time-step by calling
!                tstep.
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
      Integer, Intent(In) :: nx2
      Integer, Intent(In) :: ny2
      Integer, Intent(In) :: nz2
      Integer, Intent(In) :: nt
      Integer, Intent(In) :: mprint
      Integer, Intent(In) :: mplot
      Integer, Intent(In) :: iprint
      Integer, Intent(In) :: mxp
      Integer, Intent(In) :: mnof
      Integer, Intent(In) :: nof
      Integer, Intent(In) :: itheta
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: itimrevs
      Integer, Intent(In) :: iodds
      Integer, Intent(In) :: ihdiag
      Integer, Intent(In) :: nneut
      Integer, Intent(In) :: nprot
      Integer, Intent(In) :: ipairn
      Integer, Intent(In) :: ipairp
      Integer, Intent(In) :: iconstr
      Integer, Intent(In) :: mconstr
      Integer, Intent(In) :: imode
      Integer, Intent(In) :: icoul
      Integer, Intent(In) :: irest
      Integer, Intent(In) :: mrest
      Integer, Intent(In) :: ifixb
      Integer, Intent(In) :: niter
      Integer, Intent(In) :: nfixb
      Integer, Intent(In) :: iperc
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Integer, Intent(In), Dimension(2, mnof) :: npsif, npminf
      Integer, Intent(In) :: itrbx
      Integer(8), Intent(In) :: coulplan1, coulplan2
      Real(wp), Intent(In) :: serr
      Real(wp), Intent(In) :: terr
      Real(wp), Intent(In) :: dt
      Real(wp), Intent(In) :: t0
      Real(wp), Intent(In) :: t1
      Real(wp), Intent(In) :: t2
      Real(wp), Intent(In) :: t3
      Real(wp), Intent(In) :: t4
      Real(wp), Intent(In) :: t4p
      Real(wp), Intent(In) :: x0
      Real(wp), Intent(In) :: x1
      Real(wp), Intent(In) :: x2
      Real(wp), Intent(In) :: x3
      Real(wp), Intent(In) :: alpha
      Real(wp), Intent(In) :: ayuk
      Real(wp), Intent(In) :: vyuk
      Real(wp), Intent(In) :: v0neut
      Real(wp), Intent(In) :: v0prot
      Real(wp), Intent(In) :: rho0pr
      Real(wp), Intent(In) :: x0dmp
      Real(wp), Intent(In) :: e0dmp
      Real(wp), Intent(In) :: c0
      Real(wp), Intent(In) :: d0
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: h2m(2)
      Real(wp), Intent(In) :: cdmpx(ncolx, ncolx)
      Real(wp), Intent(In) :: cdmpy(ncoly, ncoly)
      Real(wp), Intent(In) :: cdmpz(ncolz, ncolz)
      Real(wp), Intent(In) :: ecoul(mnof)
      Real(wp), Intent(In) :: rsep
      Real(wp), Intent(In) :: xb
      Real(wp) :: roft
      Real(wp) :: rdot
      Real(wp) :: tdot
      Real(wp) :: vx
      Real(wp) :: vz
      Real(wp) :: vocc(lpsi, 2)
      Real(wp) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp) :: tau(ncolx, ncoly, ncolz, 2)
      Real(wp) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: sodens(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: spinden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: kinvden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: spincur(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: delxj(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: bmunu(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: dbmu(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: upot(ncolx, ncoly, ncolz, 2)
      Real(wp) :: bmass(ncolx, ncoly, ncolz, 2)
      Real(wp) :: xiq(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: hsigma(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: cq(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: dcq(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: tmpvec(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp) :: wyuk(ncolx, ncoly, ncolz)
      Real(wp) :: edenss(ncolx, ncoly, ncolz)
      Real(wp) :: edenst(ncolx, ncoly, ncolz)
      Real(wp) :: dbmass(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: worka(ncolx, ncoly, ncolz, 2)
      Real(wp) :: rlam(ncolx, ncoly, ncolz, 2)
      Real(wp) :: dxiq(ncolx, ncoly, ncolz, 2)
      Real(wp) :: spenrg(lpsi, 2)
      Real(wp) :: speflu(lpsi, 2)
      Real(wp) :: spkine(lpsi, 2)
      Real(wp) :: spnorm(lpsi, 2)
      Real(wp) :: spar(lpsi, 2)
      Real(wp) :: ajz(lpsi, 2)
      Real(wp) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp) :: xcm(2)
      Real(wp) :: ycm(2)
      Real(wp) :: zcm(2)
      Real(wp) :: alz(lpsi, 2)
      Real(wp) :: asz(lpsi, 2)
      Real(wp) :: centf(3, mnof)
      Real(wp) :: fermin(2), fermip(2), epairn(2), epairp(2)
      Complex(wp) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk2(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk3(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk4(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: tpsi(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: q(nx2, ny2, nz2), rho2(nx2, ny2, nz2)
      Character(len=80) :: directname
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Integer, Parameter :: ifixcm = 0
      Real(wp), Parameter :: esf = 0.0_wp
      Real(wp), Parameter :: hbc = 197.3269885804381_wp
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, nst, lpsihalf, it, ix, iy, iz, is, iyy, iquant
      Integer :: npsit(2), nst1, irest0, irwgs0, mplot0, nexadd0, nexiter0
      Real(wp) :: time, xdt, xvcmc, amc2, akfx, x
      Real(wp) :: ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehfc, ecorc, ehfy, tetrs
      Real(wp) :: ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2
      Real(wp) :: rold, roftq2, tetold, rscale, q20old, q20tot, xdot, zdot, q2in, cursq2, grsq2, tautf!, cur2q
!      Real(wp) :: cur2n, cur2p, cur2t, rhon, rhop, rhot
      Complex(wp) :: cnormb
      Real(wp) :: isoedenst(ncolx, ncoly, ncolz, 2), etheta(3)
      Real(wp) :: dqmu(2,3), dqmu0(2), dqmu0_tf(2)
      Real(wp), Dimension(:, :), Allocatable                :: vocct, explt, spenrgt, spkinet
      Real(wp), Dimension(:, :, :), Allocatable             :: estar, taus_int, tausp_int, tautf_int
      Complex(wp), Dimension(:, :, :, :, :, :), Allocatable :: psihalf
      Real(wp), Dimension(:,:,:,:,:,:), Allocatable         :: rhos, taus, clocal
      Real(wp), Dimension(:,:,:,:,:,:,:), Allocatable       :: drhos, currnts
      Real(wp), Dimension(:, :), Allocatable                :: rhoplt
!-------------------------------------------------------------------------
! allocate temp psi for density contraint (we do this with T-R invariance)
!-------------------------------------------------------------------------
      If(iconstr /= 0 .And. itimrevs /= 0) Then
         lpsihalf = lpsi / 2
         Allocate(psihalf(ncolx, ncoly, ncolz, 2, lpsihalf, 2))
         Allocate(vocct(lpsihalf, 2), spenrgt(lpsihalf, 2), spkinet(lpsihalf, 2))
      Elseif(iconstr /= 0 .And. itimrevs == 0) Then
         Allocate(psihalf(ncolx, ncoly, ncolz, 2, lpsi, 2))
         Allocate(vocct(lpsi, 2), spenrgt(lpsi, 2), spkinet(lpsi, 2))
      End If
      Allocate(estar(ncolx, ncoly, ncolz))
      Allocate(explt(ncolx, ncolz))
      estar = 0.0_wp
      Allocate( rhos(ncolx,ncoly,ncolz,2,2,3), taus(ncolx,ncoly,ncolz,2,2,3), clocal(ncolx,ncoly,ncolz,2,2,3) )
      Allocate( taus_int(2,2,3), tausp_int(2,2,3), tautf_int(2,2,3)  )
      Allocate( drhos(ncolx,ncoly,ncolz,3,2,2,3), currnts(ncolx,ncoly,ncolz,3,2,2,3) )
      Allocate( rhoplt(ncolx, ncolz) )
!-----------------------------------------------------------------------
!       product of weight's
!-----------------------------------------------------------------------
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               wxyz(ix, iy, iz) = wx(ix) * wy(iy) * wz(iz)
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        initialize densities, currents, potential etc.
!-----------------------------------------------------------------------
      If(irest /= 0) Then
         Call getpsi(psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, wcoul, &
        & mnof, tetrs, centf, npsif, npminf, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul)
      End If
!
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
!
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, iq, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, psi(1, 1, 1, 1, nst, iq), pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
      End Do
      xdt = - 1.0_wp
      Call tinfo(ncolx, ncoly, ncolz, nx2, ny2, nz2, lpsi, 1, 1, 0, 0, wx, wy, wz, xclx, xcly, xclz, h2m, t3, x3, alpha, time, npsi, npmin, &
     & vocc, spkine, spenrg, spnorm, spar, ajz, xcm, ycm, zcm, centf, roft, roftq2, rdot, tdot, rho, tau, currnt, ehft, ehf0, ehf12, &
     & ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, ehfy, mnof, nof, &
     & ifixb, vx, vz, xb, tetrs, xdt, irest, psi, niter, nfixb, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul, wcoul, &
     & estar, xvcmc, isoedenst, der1x, der1y, der1z, coulplan1, coulplan2, rho2, q, iperc, fermin, fermip, epairn, epairp)
!
      Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x, &
     & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
     & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1, &
     & coulplan2, rho2, q, iperc)
!-----------------------------------------------------------------------
!        dynamic steps begin here
!        first do the half-time step
!-----------------------------------------------------------------------
      Do it = (irest+1), (irest+nt)
         time = it * dt
         xdt = dt / 2.0_wp
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
!-----------------------------------------------------------------------
!        perform the time stepping and compute densities, currents, etc.
!-----------------------------------------------------------------------
!
! ** neutron states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(1), npsi(1)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
            Call tstep(ncolx, ncoly, ncolz, 1, mxp, terr, xdt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 1, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
!$OMP END PARALLEL DO
!
! ** proton states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(2), npsi(2)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
            Call tstep(ncolx, ncoly, ncolz, 2, mxp, terr, xdt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 2, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
!$OMP END PARALLEL DO
         xdt = - 1.0_wp
         Call tinfo(ncolx, ncoly, ncolz, nx2, ny2, nz2, lpsi, 1, 1, 0, 0, wx, wy, wz, xclx, xcly, xclz, h2m, t3, x3, alpha, time, npsi, npmin, &
        & vocc, spkine, spenrg, spnorm, spar, ajz, xcm, ycm, zcm, centf, roft, roftq2, rdot, tdot, rho, tau, currnt, ehft, ehf0, ehf12, &
        & ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, ehfy, mnof, nof, &
        & ifixb, vx, vz, xb, tetrs, xdt, irest, psi, niter, nfixb, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul, &
        & wcoul, estar, xvcmc, isoedenst, der1x, der1y, der1z, coulplan1, coulplan2, rho2, q, iperc, fermin, fermip, epairn, epairp)
!
         Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x, &
        & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
        & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1, &
        & coulplan2, rho2, q, iperc)
!-----------------------------------------------------------------------
!        full-time step
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
!-----------------------------------------------------------------------
!        perform the time stepping and compute densities, currents, etc.
!-----------------------------------------------------------------------
!
! ** neutron states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(1), npsi(1)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
            Call tstep(ncolx, ncoly, ncolz, 1, mxp, terr, dt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 1, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
            psi(:ncolx, :ncoly, :ncolz, :2, nst, 1) = pswk4(:ncolx, :ncoly, :ncolz, :2)
         End Do
!$OMP END PARALLEL DO
!
! ** proton states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(2), npsi(2)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
            Call tstep(ncolx, ncoly, ncolz, 2, mxp, terr, dt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 2, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
            psi(:ncolx, :ncoly, :ncolz, :2, nst, 2) = pswk4(:ncolx, :ncoly, :ncolz, :2)
         End Do
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------
!        calculate the expectation cnormb=<h> and norms  (printing)
!-----------------------------------------------------------------------
         If(mprint /= 0) Then
            If(iprint > 0 .And. (Mod(it, mprint) == 0 .Or. it == 1)) Then
!
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     Call hpsi(ncolx, ncoly, ncolz, upot(1, 1, 1, iq), bmass(1, 1, 1, iq), xiq(1, 1, 1, 1, iq), hsigma(1, 1, 1, &
                    & 1, iq), cq(1, 1, 1, 1, iq), dcq(1, 1, 1, 1, 1, iq), bmunu(1, 1, 1, 1, 1, iq), dbmu(1, 1, 1, 1, iq), &
                    & dbmass(1, 1, 1, 1, iq), dxiq(1, 1, 1, iq), der1x, der1y, der1z, der2x, der2y, der2z, esf, psi(1, 1, 1, 1, &
                    & nst, iq), pswk1, pswk3, itimrev)
                     Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, iq), pswk1, cnormb, 1, wxyz)
                     spenrg(nst, iq) = cnormb
                     Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, iq), psi(1, 1, 1, 1, nst, iq), cnormb, 1, wxyz)
                     spnorm(nst, iq) = cnormb
                  End Do
               End Do
!
               Call kine(ncolx, ncoly, ncolz, lpsi, npsi, npmin, der2x, der2y, der2z, spkine, psi, pswk1, wxyz)
!
               Call energy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, &
              & der1y, der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, &
              & wyuk, edenst, worka, icoul, ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, &
              & ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, ehfy, itheta, itimrev, iodds)
               Call isoenergy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, &
              & der1y, der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, &
              & wyuk, isoedenst, worka, icoul, itheta, itimrev, iodds)

               If(iconstr /= 2) Then
!                  write(42,*) ((((rho(ix,iy,iz,1)+rho(ix,iy,iz,2)),ix=1,ncolx),iy=1,ncoly),iz=1,ncolz)
               End If
!-----------------------------------------------------------------------
!        calculate localization function
!-----------------------------------------------------------------------
               rhos = 0.0_wp
               taus = 0.0_wp
               drhos = 0.0_wp
               currnts = 0.0_wp
               Call densit_loc(lpsi, npsi, npmin, ncolx, ncoly, ncolz, vocc, der1x, der1y, der1z, psi, pswk4, pswk1, rhos, taus,&
            &                  drhos, currnts, itimrev)
!
               taus_int = 0.0_wp
               tausp_int = 0.0_wp
               tautf_int = 0.0_wp
               Do iq = 1, 2
                  Do is = 1, 2
                     Do iquant = 1, 3
                        Do ix = 1, ncolx
                           Do iy = 1, ncoly
                              Do iz = 1, ncolz
                                grsq2 = drhos(ix,iy,iz,1,is,iq,iquant)**2 + drhos(ix,iy,iz,2,is,iq,iquant)**2 + &
                              &         drhos(ix,iy,iz,3,is,iq,iquant)**2
                                cursq2 = currnts(ix,iy,iz,1,is,iq,iquant)**2 + currnts(ix,iy,iz,2,is,iq,iquant)**2 + &
                              &          currnts(ix,iy,iz,3,is,iq,iquant)**2
                                tautf = 3.0_wp/5.0_wp*(6.0_wp*pi**2)**(2.0_wp/3.0_wp)*rhos(ix,iy,iz,is,iq,iquant)**(5.0_wp/3.0_wp)
                                tautf_int(is,iq,iquant) = tautf_int(is,iq,iquant) + wxyz(ix,iy,iz)*(tautf+0.5_wp*grsq2)
                                taus_int(is,iq,iquant) = taus_int(is,iq,iquant) + wxyz(ix,iy,iz)*taus(ix,iy,iz,is,iq,iquant)
                                tausp_int(is,iq,iquant) = tausp_int(is,iq,iquant) + &
                              & wxyz(ix,iy,iz)*(0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant)+cursq2/rhos(ix,iy,iz,is,iq,iquant))
                                clocal(ix,iy,iz,is,iq,iquant) = taus(ix,iy,iz,is,iq,iquant) - &
                              &                                 0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) - &
                              &                                 cursq2/rhos(ix,iy,iz,is,iq,iquant)
                                clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(1.0_wp + (clocal(ix,iy,iz,is,iq,iquant)/tautf)**2)
!                                clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(clocal(ix,iy,iz,is,iq,iquant)/tautf)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
               write(8,*) "======================="
               dqmu0(1) = 22.158369931
               dqmu0(2) = 21.642406618
               dqmu0_tf(1) = 21.32314_wp
               dqmu0_tf(2) = 20.86163_wp
               do iq = 1, 2
                  do iquant = 1, 3
                     dqmu(iq,iquant) = 0.0_wp
                     do is = 1, 2
                        dqmu(iq,iquant) = dqmu(iq,iquant) + (taus_int(is,iq,iquant)-tausp_int(is,iq,iquant))
!                        dqmu(iq,iquant) = dqmu(iq,iquant) + (tautf_int(is,iq,iquant)-tausp_int(is,iq,iquant))
                     end do
                     write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-t: ", rsep, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0(iq))
!                     write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-t: ", rsep, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0_tf(iq))
                  end do
               end do
               write(8,*) "======================="
            End If
         End If
! branch here if we want to do a density constraint calculation
         If(iconstr /= 0 .And. (Mod(it, mconstr) == 0 .Or. it == 1)) Then
! save the original time-dependent wavefunction, use psihalf instead
! since static calculations preserve T-R invariance we can use 1/2 states only
            If(itimrevs /= 0) Then
               Do iq = 1, 2
                  npsit(iq) = npsi(iq) / 2
                  Do nst = npmin(iq), npsi(iq), 2
                     nst1 = (nst+1) / 2
                     vocct(nst1, iq) = 2 * vocc(nst, iq)
                     spenrgt(nst1, iq) = spenrg(nst, iq)
                     spkinet(nst1, iq) = spkine(nst, iq)
                     psihalf(:ncolx, :ncoly, :ncolz, :2, nst1, iq) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                  End Do
               End Do
            Else
               psihalf = psi
               vocct = vocc
               spenrgt = spenrg
               spkinet = spkine
            End If
! call static calculation
            irest0 = 0
            irwgs0 = 0
            mplot0 = 0
            nexadd0 = 0
            nexiter0 = 0
            q2in = 0.0_wp ! not used here
            If(itimrevs /= 0) Then
               Call shf3d(lpsihalf, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsit, npmin, itrbx, mprint, iprint, serr, h2m, t0, t1,&
              & t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj,&
              & bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, dbmass, dxiq, der1x, der1y, der1z,&
              & der2x, der2y, der2z, spenrgt, speflu, spkinet, vocct, spnorm, spar, ajz, psihalf, pswk1, pswk2, pswk3, pswk4,&
              & cdmpx, cdmpy, cdmpz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, x0dmp, e0dmp, icoul, irest0, mrest, mplot0, tpsi,&
              & alz, asz, itimrevs, nneut, nprot, ipairn, ipairp, fermin, fermip, epairn, epairp, worka, rlam, wxyz, c0, d0, &
              & iconstr, q2in, itheta, imode, irwgs0, iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, v0prot, rho0pr, &
              & nexadd0, nexiter0, etheta, ihdiag, .FALSE.)
            Else
               Call shf3d(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, itrbx, mprint, iprint, serr, h2m, t0, t1, t2,&
              & t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj,&
              & bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, dbmass, dxiq, der1x, der1y, der1z,&
              & der2x, der2y, der2z, spenrgt, speflu, spkinet, vocct, spnorm, spar, ajz, psihalf, pswk1, pswk2, pswk3, pswk4,&
              & cdmpx, cdmpy, cdmpz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, x0dmp, e0dmp, icoul, irest0, mrest, mplot0, tpsi,&
              & alz, asz, itimrevs, nneut, nprot, ipairn, ipairp, fermin, fermip, epairn, epairp, worka, rlam, wxyz, c0, d0, &
              & iconstr, q2in, itheta, imode, irwgs0, iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, v0prot, rho0pr, &
              & nexadd0, nexiter0, etheta, ihdiag, .FALSE.)
            End If
!-----------------------------------------------------------------------
!        calculate localization function
!-----------------------------------------------------------------------
            rhos = 0.0_wp
            taus = 0.0_wp
            drhos = 0.0_wp
            currnts = 0.0_wp
            Call densit_loc(lpsi, npsi, npmin, ncolx, ncoly, ncolz, vocc, der1x, der1y, der1z, psihalf, pswk4, pswk1, rhos, taus,&
         &                  drhos, currnts, itimrev)
!
            taus_int = 0.0_wp
            tausp_int = 0.0_wp
            tautf_int = 0.0_wp
            Do iq = 1, 2
               Do is = 1, 2
                  Do iquant = 1, 3
                     Do ix = 1, ncolx
                        Do iy = 1, ncoly
                           Do iz = 1, ncolz
                             grsq2 = drhos(ix,iy,iz,1,is,iq,iquant)**2 + drhos(ix,iy,iz,2,is,iq,iquant)**2 + &
                           &         drhos(ix,iy,iz,3,is,iq,iquant)**2
                             cursq2 = currnts(ix,iy,iz,1,is,iq,iquant)**2 + currnts(ix,iy,iz,2,is,iq,iquant)**2 + &
                           &          currnts(ix,iy,iz,3,is,iq,iquant)**2
                             tautf = 3.0_wp/5.0_wp*(6.0_wp*pi**2)**(2.0_wp/3.0_wp)*rhos(ix,iy,iz,is,iq,iquant)**(5.0_wp/3.0_wp)
                             tautf_int(is,iq,iquant) = tautf_int(is,iq,iquant) + wxyz(ix,iy,iz)*(tautf+0.5_wp*grsq2)
                             taus_int(is,iq,iquant) = taus_int(is,iq,iquant) + wxyz(ix,iy,iz)*taus(ix,iy,iz,is,iq,iquant)
                             tausp_int(is,iq,iquant) = tausp_int(is,iq,iquant) + &
                           & wxyz(ix,iy,iz)*(0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) + cursq2/rhos(ix,iy,iz,is,iq,iquant))
                             clocal(ix,iy,iz,is,iq,iquant) = taus(ix,iy,iz,is,iq,iquant) - &
                           &                                 0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) - &
                           &                                 cursq2/rhos(ix,iy,iz,is,iq,iquant)
                             clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(1.0_wp + (clocal(ix,iy,iz,is,iq,iquant)/tautf)**2)
!                                clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(clocal(ix,iy,iz,is,iq,iquant)/tautf)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
            write(8,*) "======================="
            dqmu0(1) = 22.158369931
            dqmu0(2) = 21.642406618
            dqmu0_tf(1) = 21.32314_wp
            dqmu0_tf(2) = 20.86163_wp
            do iq = 1, 2
               do iquant = 1, 3
                  dqmu(iq,iquant) = 0.0_wp
                  do is = 1, 2
                     dqmu(iq,iquant) = dqmu(iq,iquant) + (taus_int(is,iq,iquant)-tausp_int(is,iq,iquant))
!                        dqmu(iq,iquant) = dqmu(iq,iquant) + (tautf_int(is,iq,iquant)-tausp_int(is,iq,iquant))
                  end do
                  write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-dcfhf: ", rsep, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0(iq))
!                     write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-dcfhf: ", rsep, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0_tf(iq))
               end do
            end do
            write(8,*) "======================="

! restore back densities/currents/potentials etc from the original time-dependent wavefunction
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
!
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                  pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                  Call densit(lpsi, ncolx, ncoly, ncolz, nst, iq, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, &
                 & der1x, der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
               End Do
            End Do
!-----------------------------------------------------------------------
!        calculate the excitation energy (we choose to do this in tinfo with tke instead)
!-----------------------------------------------------------------------
            estar = 0.0_wp
!            Do iz = 1, ncolz
!               Do iy = 1, ncoly
!                  Do ix = 1, ncolx
!                     cur2t = (currnt(ix, iy, iz, 1, 1)+currnt(ix, iy, iz, 1, 2)) ** 2 + (currnt(ix, iy, iz, 2, 1)+currnt(ix, iy,&
!                    & iz, 2, 2)) ** 2 + (currnt(ix, iy, iz, 3, 1)+currnt(ix, iy, iz, 3, 2)) ** 2
!                     cur2n = currnt(ix, iy, iz, 1, 1) ** 2 + currnt(ix, iy, iz, 2, 1) ** 2 + currnt(ix, iy, iz, 3, 1) ** 2
!                     cur2p = currnt(ix, iy, iz, 1, 2) ** 2 + currnt(ix, iy, iz, 2, 2) ** 2 + currnt(ix, iy, iz, 3, 2) ** 2
!                     rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2) + 1.0e-25_wp
!                     rhon = rho(ix, iy, iz, 1) + 1.0e-25_wp
!                     rhop = rho(ix, iy, iz, 2) + 1.0e-25_wp
!                     estar(ix,iy,iz) = estar(ix,iy,iz) + (h2m(1)+h2m(2))/2.0_wp*cur2t/(rhot + 1.0e-25_wp)
!                     estar(ix, iy, iz) = estar(ix, iy, iz) + h2m(1) * cur2n / rhon + h2m(2) * cur2p / rhop
!                  End Do
!               End Do
!            End Do
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     estar(ix, iy, iz) = edenst(ix, iy, iz) - edenss(ix, iy, iz)! - estar(ix, iy, iz) !don't subtract tke here
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        plot the excitation energy
!-----------------------------------------------------------------------
            iyy = ncoly / 2
            explt(:ncolx, :ncolz) = estar(:ncolx, iyy, :ncolz)
!            call grid_plot(explt,"Excitation Density",xclx,xclz,ncolx,ncolz,16,2)
!            write(44,*) (((estar(ix,iy,iz),ix=1,ncolx),iy=1,ncoly),iz=1,ncolz)
            Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk,&
           & der1x, der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj,&
           & bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds,&
           & wxyz, coulplan1, coulplan2, rho2, q, iperc)
!
! end constraint branch
!
         End If
!
         Call tinfo(ncolx, ncoly, ncolz, nx2, ny2, nz2, lpsi, mprint, mplot, iprint, it, wx, wy, wz, xclx, xcly, xclz, h2m, t3, x3, alpha, time,&
          & npsi, npmin, vocc, spkine, spenrg, spnorm, spar, ajz, xcm, ycm, zcm, centf, roft, roftq2, rdot, tdot, rho, tau, currnt, ehft,&
          & ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc,&
          & ehfy, mnof, nof, ifixb, vx, vz, xb, tetrs, dt, irest, psi, niter, nfixb, rold, tetold, rscale, q20old, q20tot, xdot,&
          & zdot, ecoul, wcoul, estar, xvcmc, isoedenst, der1x, der1y, der1z, coulplan1, coulplan2, rho2, q, iperc, fermin, fermip, epairn, epairp)
!-------------------------------------------------------------------------------------
! check if fragments are separated by a distance larger than the initial one and stop
!-------------------------------------------------------------------------------------
         If((roftq2 > 1.2*rsep) .And. (rdot > 0.0_wp) .And. (it > 2)) Then
            Write(8, '(//,a,//)') ' Final separation distance reached'
            Write(8, '(//,a,//)') ' Printing final time step'
            Exit
         End If
          ! correct if the c.m. is moving with an opposite boost (if turned on in parameter statement)
         If(ifixcm == 1) then
            If(Mod(it, mprint) == 0) Then
               amc2 = hbc * hbc / (h2m(1)+h2m(2))
               akfx = - amc2 * xvcmc / hbc
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     Do is = 1, 2
                        Do iz = 1, ncolz
                           Do iy = 1, ncoly
                              Do ix = 1, ncolx
                                 x = xclx(ix)
                                 psi(ix, iy, iz, is, nst, iq) = psi(ix, iy, iz, is, nst, iq) * Exp(eye*akfx*x)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End If
         End If
!-----------------------------------------------------------------------
!        compute densities, currents, potentials etc.
!-----------------------------------------------------------------------
         Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x,&
        & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu,&
        & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1,&
        & coulplan2, rho2, q, iperc)
         If(Mod(it, mrest) == 0) Then
            Call putpsi(psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
           & wcoul, mnof, tetrs, centf, npsif, npminf, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul)
         End If
!
         If(nt > 1 .and. Mod(it, mprint) == 0) then
            Call write_density(it, time, ncolx, ncoly, ncolz, xclx, xcly, xclz, rho, currnt, clocal, directname)
         End if
! end time loop
      End Do
! Always do one more time step if sep is reached
      If((roftq2 > 1.2*rsep) .And. (rdot > 0.0_wp) .And. (it > 2)) Then
         it = it + 1
         time = it * dt
         xdt = dt / 2.0_wp
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
!-----------------------------------------------------------------------
!        perform the time stepping and compute densities, currents, etc.
!-----------------------------------------------------------------------
!
! ** neutron states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(1), npsi(1)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
            Call tstep(ncolx, ncoly, ncolz, 1, mxp, terr, xdt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 1, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
!$OMP END PARALLEL DO
!
! ** proton states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(2), npsi(2)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
            Call tstep(ncolx, ncoly, ncolz, 2, mxp, terr, xdt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 2, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
!$OMP END PARALLEL DO
         xdt = - 1.0_wp
         Call tinfo(ncolx, ncoly, ncolz, nx2, ny2, nz2, lpsi, 1, 1, 0, 0, wx, wy, wz, xclx, xcly, xclz, h2m, t3, x3, alpha, time, npsi, npmin, &
        & vocc, spkine, spenrg, spnorm, spar, ajz, xcm, ycm, zcm, centf, roft, roftq2, rdot, tdot, rho, tau, currnt, ehft, ehf0, ehf12, &
        & ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, ehfy, mnof, nof, &
        & ifixb, vx, vz, xb, tetrs, xdt, irest, psi, niter, nfixb, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul, &
        & wcoul, estar, xvcmc, isoedenst, der1x, der1y, der1z, coulplan1, coulplan2, rho2, q, iperc, fermin, fermip, epairn, epairp)
!
         Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x, &
        & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
        & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1, &
        & coulplan2, rho2, q, iperc)
!-----------------------------------------------------------------------
!        full-time step
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
!-----------------------------------------------------------------------
!        perform the time stepping and compute densities, currents, etc.
!-----------------------------------------------------------------------
!
! ** neutron states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(1), npsi(1)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
            Call tstep(ncolx, ncoly, ncolz, 1, mxp, terr, dt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 1, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
            psi(:ncolx, :ncoly, :ncolz, :2, nst, 1) = pswk4(:ncolx, :ncoly, :ncolz, :2)
         End Do
!$OMP END PARALLEL DO
!
! ** proton states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk2,pswk3,pswk4)&
!$OMP SCHEDULE(STATIC) REDUCTION(+:rho, tau, currnt, sodens, spinden, kinvden, spincur, rhos, taus, drhos, currnts)
         Do nst = npmin(2), npsi(2)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
            Call tstep(ncolx, ncoly, ncolz, 2, mxp, terr, dt, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, dxiq, &
           & der1x, der1y, der1z, der2x, der2y, der2z, pswk4, pswk1, pswk2, pswk3, wxyz, itimrev)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 2, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
            psi(:ncolx, :ncoly, :ncolz, :2, nst, 2) = pswk4(:ncolx, :ncoly, :ncolz, :2)
         End Do
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------
!        calculate the expectation cnormb=<h> and norms  (printing)
!-----------------------------------------------------------------------
         If(mprint /= 0) Then
            If(iprint > 0) Then
!
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     Call hpsi(ncolx, ncoly, ncolz, upot(1, 1, 1, iq), bmass(1, 1, 1, iq), xiq(1, 1, 1, 1, iq), hsigma(1, 1, 1, &
                    & 1, iq), cq(1, 1, 1, 1, iq), dcq(1, 1, 1, 1, 1, iq), bmunu(1, 1, 1, 1, 1, iq), dbmu(1, 1, 1, 1, iq), &
                    & dbmass(1, 1, 1, 1, iq), dxiq(1, 1, 1, iq), der1x, der1y, der1z, der2x, der2y, der2z, esf, psi(1, 1, 1, 1, &
                    & nst, iq), pswk1, pswk3, itimrev)
                     Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, iq), pswk1, cnormb, 1, wxyz)
                     spenrg(nst, iq) = cnormb
                     Call psnorm(ncolx, ncoly, ncolz, psi(1, 1, 1, 1, nst, iq), psi(1, 1, 1, 1, nst, iq), cnormb, 1, wxyz)
                     spnorm(nst, iq) = cnormb
                  End Do
               End Do
!
               Call kine(ncolx, ncoly, ncolz, lpsi, npsi, npmin, der2x, der2y, der2z, spkine, psi, pswk1, wxyz)
!
               Call energy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, &
              & der1y, der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, &
              & wyuk, edenst, worka, icoul, ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, &
              & ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, ehfy, itheta, itimrev, iodds)
               Call isoenergy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, &
              & der1y, der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, &
              & wyuk, isoedenst, worka, icoul, itheta, itimrev, iodds)

               If(iconstr /= 2) Then
!                  write(42,*) ((((rho(ix,iy,iz,1)+rho(ix,iy,iz,2)),ix=1,ncolx),iy=1,ncoly),iz=1,ncolz)
               End If
!-----------------------------------------------------------------------
!        calculate localization function
!-----------------------------------------------------------------------
               rhos = 0.0_wp
               taus = 0.0_wp
               drhos = 0.0_wp
               currnts = 0.0_wp
               Call densit_loc(lpsi, npsi, npmin, ncolx, ncoly, ncolz, vocc, der1x, der1y, der1z, psi, pswk4, pswk1, rhos, taus,&
            &                  drhos, currnts, itimrev)
!
               taus_int =0.0_wp
               tausp_int = 0.0_wp
               tautf_int = 0.0_wp
               Do iq = 1, 2
                  Do is = 1, 2
                     Do iquant = 1, 3
                        Do ix = 1, ncolx
                           Do iy = 1, ncoly
                              Do iz = 1, ncolz
                                grsq2 = drhos(ix,iy,iz,1,is,iq,iquant)**2 + drhos(ix,iy,iz,2,is,iq,iquant)**2 + &
                              &         drhos(ix,iy,iz,3,is,iq,iquant)**2
                                cursq2 = currnts(ix,iy,iz,1,is,iq,iquant)**2 + currnts(ix,iy,iz,2,is,iq,iquant)**2 + &
                              &          currnts(ix,iy,iz,3,is,iq,iquant)**2
                                tautf = 3.0_wp/5.0_wp*(6.0_wp*pi**2)**(2.0_wp/3.0_wp)*rhos(ix,iy,iz,is,iq,iquant)**(5.0_wp/3.0_wp)
                                tautf_int(is,iq,iquant) = tautf_int(is,iq,iquant) + wxyz(ix,iy,iz)*(tautf+0.5_wp*grsq2)
                                taus_int(is,iq,iquant) = taus_int(is,iq,iquant) + wxyz(ix,iy,iz)*taus(ix,iy,iz,is,iq,iquant)
                                tausp_int(is,iq,iquant) = tausp_int(is,iq,iquant) + &
                              & wxyz(ix,iy,iz)*(0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant)+cursq2/rhos(ix,iy,iz,is,iq,iquant))
                                clocal(ix,iy,iz,is,iq,iquant) = taus(ix,iy,iz,is,iq,iquant) - &
                              &                                 0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) - &
                              &                                 cursq2/rhos(ix,iy,iz,is,iq,iquant)
                                clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(1.0_wp + (clocal(ix,iy,iz,is,iq,iquant)/tautf)**2)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
               write(8,*) "======================="
               dqmu0(1) = 22.158369931
               dqmu0(2) = 21.642406618
               dqmu0_tf(1) = 21.32314_wp
               dqmu0_tf(2) = 20.86163_wp
               do iq = 1, 2
                  do iquant = 1, 3
                     dqmu(iq,iquant) = 0.0_wp
                     do is = 1, 2
                        dqmu(iq,iquant) = dqmu(iq,iquant) + (taus_int(is,iq,iquant)-tausp_int(is,iq,iquant))
!                        dqmu(iq,iquant) = dqmu(iq,iquant) + (tautf_int(is,iq,iquant)-tausp_int(is,iq,iquant))
                     end do
                     write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-t: ", roft, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0(iq))
!                     write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-t: ", roft, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0_tf(iq))
                  end do
               end do
               write(8,*) "======================="
            End If
         End If
! branch here if we want to do a density constraint calculation
         If(iconstr /= 0 .And. (Mod(it, mconstr) == 0 .Or. it == 1)) Then
! save the original time-dependent wavefunction, use psihalf instead
! since static calculations preserve T-R invariance we can use 1/2 states only
            If(itimrevs /= 0) Then
               Do iq = 1, 2
                  npsit(iq) = npsi(iq) / 2
                  Do nst = npmin(iq), npsi(iq), 2
                     nst1 = (nst+1) / 2
                     vocct(nst1, iq) = 2 * vocc(nst, iq)
                     spenrgt(nst1, iq) = spenrg(nst, iq)
                     spkinet(nst1, iq) = spkine(nst, iq)
                     psihalf(:ncolx, :ncoly, :ncolz, :2, nst1, iq) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                  End Do
               End Do
            Else
               psihalf = psi
               vocct = vocc
               spenrgt = spenrg
               spkinet = spkine
            End If
! call static calculation
            irest0 = 0
            irwgs0 = 0
            mplot0 = 0
            nexadd0 = 0
            nexiter0 = 0
            q2in = 0.0_wp ! not used here
            If(itimrevs /= 0) Then
               Call shf3d(lpsihalf, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsit, npmin, itrbx, mprint, iprint, serr, h2m, t0, t1,&
              & t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj,&
              & bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, dbmass, dxiq, der1x, der1y, der1z,&
              & der2x, der2y, der2z, spenrgt, speflu, spkinet, vocct, spnorm, spar, ajz, psihalf, pswk1, pswk2, pswk3, pswk4,&
              & cdmpx, cdmpy, cdmpz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, x0dmp, e0dmp, icoul, irest0, mrest, mplot0, tpsi,&
              & alz, asz, itimrevs, nneut, nprot, ipairn, ipairp, fermin, fermip, epairn, epairp, worka, rlam, wxyz, c0, d0, &
              & iconstr, q2in, itheta, imode, irwgs0, iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, v0prot, rho0pr, &
              & nexadd0, nexiter0, etheta, ihdiag, .FALSE.)
            Else
               Call shf3d(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, itrbx, mprint, iprint, serr, h2m, t0, t1, t2,&
              & t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj,&
              & bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edenss, dbmass, dxiq, der1x, der1y, der1z,&
              & der2x, der2y, der2z, spenrgt, speflu, spkinet, vocct, spnorm, spar, ajz, psihalf, pswk1, pswk2, pswk3, pswk4,&
              & cdmpx, cdmpy, cdmpz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, x0dmp, e0dmp, icoul, irest0, mrest, mplot0, tpsi,&
              & alz, asz, itimrevs, nneut, nprot, ipairn, ipairp, fermin, fermip, epairn, epairp, worka, rlam, wxyz, c0, d0, &
              & iconstr, q2in, itheta, imode, irwgs0, iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, v0prot, rho0pr, &
              & nexadd0, nexiter0, etheta, ihdiag, .FALSE.)
            End If
!-----------------------------------------------------------------------
!        calculate localization function
!-----------------------------------------------------------------------
            rhos = 0.0_wp
            taus = 0.0_wp
            drhos = 0.0_wp
            currnts = 0.0_wp
            Call densit_loc(lpsi, npsi, npmin, ncolx, ncoly, ncolz, vocc, der1x, der1y, der1z, psihalf, pswk4, pswk1, rhos, taus,&
         &                  drhos, currnts, itimrev)
!
            taus_int = 0.0_wp
            tausp_int = 0.0_wp
            tautf_int = 0.0_wp
            Do iq = 1, 2
               Do is = 1, 2
                  Do iquant = 1, 3
                     Do ix = 1, ncolx
                        Do iy = 1, ncoly
                           Do iz = 1, ncolz
                             grsq2 = drhos(ix,iy,iz,1,is,iq,iquant)**2 + drhos(ix,iy,iz,2,is,iq,iquant)**2 + &
                           &         drhos(ix,iy,iz,3,is,iq,iquant)**2
                             cursq2 = currnts(ix,iy,iz,1,is,iq,iquant)**2 + currnts(ix,iy,iz,2,is,iq,iquant)**2 + &
                           &          currnts(ix,iy,iz,3,is,iq,iquant)**2
                             tautf = 3.0_wp/5.0_wp*(6.0_wp*pi**2)**(2.0_wp/3.0_wp)*rhos(ix,iy,iz,is,iq,iquant)**(5.0_wp/3.0_wp)
                             tautf_int(is,iq,iquant) = tautf_int(is,iq,iquant) + wxyz(ix,iy,iz)*(tautf+0.5_wp*grsq2)
                             taus_int(is,iq,iquant) = taus_int(is,iq,iquant) + wxyz(ix,iy,iz)*taus(ix,iy,iz,is,iq,iquant)
                             tausp_int(is,iq,iquant) = tausp_int(is,iq,iquant) + &
                           & wxyz(ix,iy,iz)*(0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) + cursq2/rhos(ix,iy,iz,is,iq,iquant))
                             clocal(ix,iy,iz,is,iq,iquant) = taus(ix,iy,iz,is,iq,iquant) - &
                           &                                 0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) - &
                           &                                 cursq2/rhos(ix,iy,iz,is,iq,iquant)
                             clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(1.0_wp + (clocal(ix,iy,iz,is,iq,iquant)/tautf)**2)
!                                clocal(ix,iy,iz,is,iq,iquant) = 1.0_wp/(clocal(ix,iy,iz,is,iq,iquant)/tautf)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
            write(8,*) "======================="
            dqmu0(1) = 22.158369931
            dqmu0(2) = 21.642406618
            dqmu0_tf(1) = 21.32314_wp
            dqmu0_tf(2) = 20.86163_wp
            do iq = 1, 2
               do iquant = 1, 3
                  dqmu(iq,iquant) = 0.0_wp
                  do is = 1, 2
                     dqmu(iq,iquant) = dqmu(iq,iquant) + (taus_int(is,iq,iquant)-tausp_int(is,iq,iquant))
!                        dqmu(iq,iquant) = dqmu(iq,iquant) + (tautf_int(is,iq,iquant)-tausp_int(is,iq,iquant))
                  end do
                  write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-dcfhf: ", rsep, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0(iq))
!                     write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-dcfhf: ", rsep, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0_tf(iq))
               end do
            end do
            write(8,*) "======================="
! restore back densities/currents/potentials etc from the original time-dependent wavefunction
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
!
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                  pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                  Call densit(lpsi, ncolx, ncoly, ncolz, nst, iq, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, &
                 & der1x, der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
               End Do
            End Do
!-----------------------------------------------------------------------
            estar = 0.0_wp
!           Do iz = 1, ncolz
!                       Do iy = 1, ncoly
!                              Do ix = 1, ncolx
!                                 cur2t = (currnt(ix, iy, iz, 1, 1)+currnt(ix, iy, iz, 1, 2)) ** 2 + (currnt(ix, iy, iz, 2, 1)+currnt(ix, iy,&
!                     & iz, 2, 2)) ** 2 + (currnt(ix, iy, iz, 3, 1)+currnt(ix, iy, iz, 3, 2)) ** 2
!                     cur2n = currnt(ix, iy, iz, 1, 1) ** 2 + currnt(ix, iy, iz, 2, 1) ** 2 + currnt(ix, iy, iz, 3, 1) ** 2
!                     cur2p = currnt(ix, iy, iz, 1, 2) ** 2 + currnt(ix, iy, iz, 2, 2) ** 2 + currnt(ix, iy, iz, 3, 2) ** 2
!                     rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2) + 1.0e-25_wp
!                     rhon = rho(ix, iy, iz, 1) + 1.0e-25_wp
!                     rhop = rho(ix, iy, iz, 2) + 1.0e-25_wp
!                     estar(ix,iy,iz) = estar(ix,iy,iz) + (h2m(1)+h2m(2))/2.0_wp*cur2t/(rhot + 1.0e-25_wp)
!                     estar(ix, iy, iz) = estar(ix, iy, iz) + h2m(1) * cur2n / rhon + h2m(2) * cur2p / rhop
!                  End Do
!               End Do
!            End Do
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     estar(ix, iy, iz) = edenst(ix, iy, iz) - edenss(ix, iy, iz)! - estar(ix, iy, iz) !don't subtract tke here
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        plot the excitation energy
!-----------------------------------------------------------------------
            iyy = ncoly / 2
            explt(:ncolx, :ncolz) = estar(:ncolx, iyy, :ncolz)
!            call grid_plot(explt,"Excitation Density",xclx,xclz,ncolx,ncolz,16,2)
!            write(44,*) (((estar(ix,iy,iz),ix=1,ncolx),iy=1,ncoly),iz=1,ncolz)
            Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk,&
           & der1x, der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj,&
           & bmunu, dbmu, upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds,&
           & wxyz, coulplan1, coulplan2, rho2, q, iperc)
!
! end constraint branch
!
         End If
!
         Call tinfo(ncolx, ncoly, ncolz, nx2, ny2, nz2, lpsi, it, it, iprint, it, wx, wy, wz, xclx, xcly, xclz, h2m, t3, x3, alpha, time,&
          & npsi, npmin, vocc, spkine, spenrg, spnorm, spar, ajz, xcm, ycm, zcm, centf, roft, roftq2, rdot, tdot, rho, tau, currnt, ehft,&
          & ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc,&
          & ehfy, mnof, nof, ifixb, vx, vz, xb, tetrs, dt, irest, psi, niter, nfixb, rold, tetold, rscale, q20old, q20tot, xdot,&
          & zdot, ecoul, wcoul, estar, xvcmc, isoedenst, der1x, der1y, der1z, coulplan1, coulplan2, rho2, q, iperc, fermin, fermip, epairn, epairp)
          ! correct if the c.m. is moving with an opposite boost (if turned on in parameter statement)
         If(ifixcm == 1) then
               amc2 = hbc * hbc / (h2m(1)+h2m(2))
               akfx = - amc2 * xvcmc / hbc
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     Do is = 1, 2
                        Do iz = 1, ncolz
                           Do iy = 1, ncoly
                              Do ix = 1, ncolx
                                 x = xclx(ix)
                                 psi(ix, iy, iz, is, nst, iq) = psi(ix, iy, iz, is, nst, iq) * Exp(eye*akfx*x)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
         End If
!-----------------------------------------------------------------------
!        compute densities, currents, potentials etc.
!-----------------------------------------------------------------------
         Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x,&
        & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu,&
        & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1,&
        & coulplan2, rho2, q, iperc)
         Call putpsi(psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
        & wcoul, mnof, tetrs, centf, npsif, npminf, rold, tetold, rscale, q20old, q20tot, xdot, zdot, ecoul)
!
         If(nt > 1) then
            Call write_density(it, time, ncolx, ncoly, ncolz, xclx, xcly, xclz, rho, currnt, clocal, directname)
         End if
      End If
!
      If(iconstr /= 0) Then
         Deallocate(psihalf, vocct, spenrgt, spkinet)
      End If
      Deallocate(estar, explt)
      Deallocate( rhos, taus, clocal, drhos, currnts, rhoplt, taus_int, tausp_int, tautf_int )
!
      Return
End Subroutine tevolv
