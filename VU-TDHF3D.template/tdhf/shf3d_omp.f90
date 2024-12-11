Subroutine shf3d(lpsi, ncolx, ncoly, ncolz, nx2, ny2, nz2, npsi, npmin, itrbx0, mprint, iprint, serr, h2m, t0, t1, t2, t3, t4, &
& t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, upot, bmass, &
& xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, edens, dbmass, dxiq, der1x, der1y, der1z, der2x, der2y, der2z, spenrg, speflu, &
& spkine, vocc, spnorm, spar, ajz, psi, pswk1, pswk2, pswk3, pswk4, cdmpx, cdmpy, cdmpz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, &
& zcm, x0dmp, e0dmp, icoul, irest, mrest, mplot, tpsi, alz, asz, itimrev, nneut, nprot, ipairn, ipairp, fermin, fermip, epairn, &
& epairp, worka, rlam, wxyz, c0, d0, iconstr, q2in, itheta, imode, irwgs, iodds, coulplan1, coulplan2, rho2, q, iperc, v0neut, &
& v0prot, rho0pr, nexadd, nexiter, etheta, ihdiag, isrestart)
!-----------------------------------------------------------------------
!     shf3d   =  the static hartree-fock solution
!-----------------------------------------------------------------------
      Use lapack95
      Use f95_precision
      Implicit None
      Integer, Parameter  :: wp = Kind(1.0D0)
      Integer, Parameter  :: mshmid = 1
      Integer, Parameter  :: mhdiag = 20
      Real(wp), Parameter :: q2eps  = 5.0e-7_wp !alternate converge criteria for q20 constraint
      Real(wp), Parameter :: q2c0   = 0.2_wp    !c0 coefficient of q20 constraint lambda iteration
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: nx2
      Integer, Intent(In) :: ny2
      Integer, Intent(In) :: nz2
      Integer, Intent(In) :: mprint
      Integer, Intent(In) :: iprint
      Integer, Intent(In) :: icoul
      Integer, Intent(In) :: irest
      Integer, Intent(In) :: mrest
      Integer, Intent(In) :: mplot
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: iodds
      Integer, Intent(In) :: ihdiag
      Integer, Intent(In) :: nneut
      Integer, Intent(In) :: nprot
      Integer, Intent(In) :: ipairn
      Integer, Intent(In) :: ipairp
      Integer, Intent(In) :: iconstr
      Integer, Intent(In) :: itheta
      Integer, Intent(In) :: imode
      Integer, Intent(In) :: irwgs
      Integer, Intent(In) :: itrbx0
      Integer, Intent(In) :: iperc
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Integer, Intent(In) :: nexadd
      Integer, Intent(In) :: nexiter
      Integer(8), Intent(In) :: coulplan1, coulplan2
      Real(wp), Intent(In) :: serr
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
      Real(wp), Intent(In) :: q2in
      Real(wp), Intent(In) :: h2m(2)
      Real(wp), Intent(In) :: cdmpx(ncolx, ncolx)
      Real(wp), Intent(In) :: cdmpy(ncoly, ncoly)
      Real(wp), Intent(In) :: cdmpz(ncolz, ncolz)
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
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
      Real(wp), Intent(Inout) :: etheta(3)
      Real(wp), Intent(Out) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp) :: wyuk(ncolx, ncoly, ncolz)
      Real(wp) :: edens(ncolx, ncoly, ncolz)
      Real(wp) :: dbmass(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: dxiq(ncolx, ncoly, ncolz, 2)
      Real(wp) :: spenrg(lpsi, 2)
      Real(wp) :: speflu(lpsi, 2)
      Real(wp) :: spkine(lpsi, 2)
      Real(wp) :: vocc(lpsi, 2)
      Real(wp) :: spnorm(lpsi, 2)
      Real(wp) :: spar(lpsi, 2)
      Real(wp) :: ajz(lpsi, 2)
      Real(wp) :: xcm(2)
      Real(wp) :: ycm(2)
      Real(wp) :: zcm(2)
      Real(wp) :: alz(lpsi, 2)
      Real(wp) :: asz(lpsi, 2)
      Real(wp) :: worka(ncolx, ncoly, ncolz, 2)
      Real(wp) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp) :: rlam(ncolx, ncoly, ncolz, 2)
      Real(wp) :: work(lpsi)
      Real(wp) :: fermin(2), fermip(2), epairn(2), epairp(2)
      Complex(wp) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp) :: psi2(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk2(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk3(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk4(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: tpsi(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: hmat(lpsi, lpsi, 2)
      Complex(wp) :: q(nx2, ny2, nz2), rho2(nx2, ny2, nz2), cnormb
      Logical, Intent(In) :: isrestart
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ix, iy, iz, itrs, iq, is, nst, nst2, ioprint, idp, jfirst, ifirst, itrbx, ireord, iquant!, il
      Real(wp)                                      :: isoedens(ncolx, ncoly, ncolz, 2)
      Real(wp), Dimension(:, :, :), Allocatable     :: damp, taus_int, tausp_int, tautf_int
      Real(wp), Dimension(:, :, :, :), Allocatable  :: rho0, rhoold, rhonew, drho
      Real(wp), Dimension(:, :, :, :), Allocatable  :: rlamold, drlam, rlamcof
      Real(wp), Dimension(:, :), Allocatable        :: rhoplt
      Real(wp), Dimension(:, :), Allocatable        :: erear3, erearc, erearl
      Real(wp), Dimension(:,:,:,:,:,:), Allocatable   :: rhos, taus, clocal
      Real(wp), Dimension(:,:,:,:,:,:,:), Allocatable :: drhos, currnts
      Real(wp) :: efluct, delesum, sumflu, spe, denerg, ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehfc, ehf, ecorc, ecorr, ecorl, &
     & ehfy, gapn, gapp
      Real(wp) :: ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2
      Real(wp) :: q20old, q202old, q20new, q202new, delq, xlam, dellam, x0act
      Real(wp) :: q20(2), q22(2),q20tot, q22tot, qerr, varqsq, dxq, dzq, grsq2, cursq2, tautf
      Real(wp) :: qxold(7), qx2old(7), qxnew(7), qx2new(7), qxin(7), xlamx(7), delqx(7), dellamx(7), varqsqx(7)
      Real(wp) :: v2(lpsi,2), isoscalarenergy, isovectorenergy
      Real(wp) :: dqmu(2,3), dqmu0(2), dqmu0_tf(2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
!
      Save jfirst
      Save ifirst
      Data jfirst / 0 /
      Data ifirst / 0 /
!
! use this if different number of iterations needed for each nuclei
      if(ifirst == 0) then
!        itrbx = 1500
         itrbx = itrbx0
         ifirst = 1
      else
         itrbx = itrbx0
      end if
!
! initialize local arrays for density constraint
      If(iconstr == 2) Then
         Allocate(rho0(ncolx, ncoly, ncolz, 2), rhoold(ncolx, ncoly, ncolz, 2))
         Allocate(drho(ncolx, ncoly, ncolz, 2), rhonew(ncolx, ncoly, ncolz, 2))
         Allocate(rlamold(ncolx, ncoly, ncolz, 2))
         Allocate(drlam(ncolx, ncoly, ncolz, 2), rlamcof(ncolx, ncoly, ncolz, 2))
      End If
      Allocate(damp(ncolx, ncoly, ncolz))
      Allocate(taus_int(2,2,3), tausp_int(2,2,3), tautf_int(2,2,3))
      Allocate(erear3(lpsi, 2), erearc(lpsi, 2), erearl(lpsi, 2))
      Allocate(rhos(ncolx,ncoly,ncolz,2,2,3), taus(ncolx,ncoly,ncolz,2,2,3), clocal(ncolx,ncoly,ncolz,2,2,3))
      Allocate(drhos(ncolx,ncoly,ncolz,3,2,2,3), currnts(ncolx,ncoly,ncolz,3,2,2,3))
! if we want to use the lambda(r) from previous iteration in the next DC run
      If(jfirst == 0) Then
         jfirst = 1
         rlam = 0.0_wp
      End If
! comment next line if you want to activate the above block
      rlam = 0.0_wp
      x0act = x0dmp
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
!        initial density for density constraint
!        For DC-TDHF runs with itimrevs=1 we cannot use psi to build rho
!        for DC initialization (keep in mind if changing below rho0).
!-----------------------------------------------------------------------
      If(iconstr == 2) Then
         rho0 = rho
         rhoold = rho0
      End If
!-----------------------------------------------------------------------
!        initialize some single-particle arrays
!-----------------------------------------------------------------------
      efluct = 1.0_wp
      delesum = 0.0_wp
      sumflu = 0.0_wp
      itrs = 0
!-----------------------------------------------------------------------
!        perform a gram-schmidt orthogonalization of all states
!-----------------------------------------------------------------------
      If(irest == 0 .And. irwgs == 0) Then
         Do iq = 1, 2
            Do nst = npmin(iq), npsi(iq)
               spenrg(nst, iq) = 0.0_wp
               speflu(nst, iq) = 0.0_wp
               spnorm(nst, iq) = 0.0_wp
               sumflu = sumflu + vocc(nst, iq) * efluct
            End Do
         End Do
         Call schmid(lpsi, ncolx, ncoly, ncolz, npsi, npmin, spnorm, psi, tpsi, pswk1, pswk2, itimrev, wxyz)
      Else If(irest /= 0 .And. irwgs == 0) Then
         Call rdtpsi(12, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
        & wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam, itimrev)
         itrs = irest
      Else If(irest == 0 .And. irwgs == 1) Then
         Call rdtpsi(15, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
        & wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam, itimrev)
      Else
         Write(8, '(//,a,//)') "Wrong combination of irest and irwgs entered."
         Stop
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
!
      Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x, &
     & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
     & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1, &
     & coulplan2, rho2, q, iperc)
!
      If(irest == 0) Then
         ioprint = 1
         Do iq = 1, 2
            Do nst = npmin(iq), npsi(iq)
               spe = spenrg(nst, iq)
               Call grstep(ncolx, ncoly, ncolz, iq, mprint, ioprint, itrs, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, &
              & dbmass, dxiq, der1x, der1y, der1z, der2x, der2y, der2z, spe, efluct, denerg, psi(1, 1, 1, 1, nst, iq), pswk1, &
              & pswk2, pswk3, cdmpx, cdmpy, cdmpz, x0act, itrbx, irest, wxyz, xcm, ycm, zcm, xclx, xcly, xclz, xlam, xlamx, damp, 0, idp, &
              & itimrev, cnormb)
               spenrg(nst, iq) = spe
               speflu(nst, iq) = efluct
               sumflu = sumflu + vocc(nst, iq) * efluct
               delesum = delesum + vocc(nst, iq) * denerg
               epairn = 0.0_wp
               epairp = 0.0_wp
               gapn = 0.0_wp
               gapp = 0.0_wp
               fermin = 0.0_wp
               fermip = 0.0_wp
            End Do
         End Do
         If(isrestart) Then
            If(itimrev == 1) Then
               v2 = vocc/2.0_wp
            Else
               v2 = vocc
            End If
            itrs = 100
            Call pairing(ipairn, ipairp, itrs, irest, ncolx, ncoly, ncolz, lpsi, npmin, npsi, nneut, nprot, epairn, epairp,&
                         gapn, gapp, fermin, fermip, spenrg, vocc, v2, wxyz, rho, v0prot, v0neut, rho0pr, psi, nexadd, itimrev)
         End If
         sumflu = sumflu / (nneut+nprot)
      Else
!-----------------------------------------------------------------------
!     perform pairing calculation
!-----------------------------------------------------------------------
         If(imode == 0 .and. irwgs /= 1) Then
            If(itimrev == 1) Then
               v2 = vocc/2.0_wp
            Else
               v2 = vocc
            End If
            Call pairing(ipairn, ipairp, 0, irest, ncolx, ncoly, ncolz, lpsi, npmin, npsi, nneut, nprot, epairn, epairp,&
                        &gapn, gapp, fermin, fermip, spenrg, vocc, v2, wxyz, rho, v0prot, v0neut, rho0pr, psi, nexadd, itimrev)
         End If
      End If
      If(imode == 0 .and. irwgs == 1) Then
         If(itimrev == 1) Then
            v2 = vocc/2.0_wp
         Else
            v2 = vocc
         End If
         itrs = 100
         Call pairing(ipairn, ipairp, itrs, irest, ncolx, ncoly, ncolz, lpsi, npmin, npsi, nneut, nprot, epairn, epairp,&
                     &gapn, gapp, fermin, fermip, spenrg, vocc, v2, wxyz, rho, v0prot, v0neut, rho0pr, psi, nexadd, itimrev)
      End If
!-----------------------------------------------------------------------
!        perform a gram-schmidt orthogonalization of all states
!-----------------------------------------------------------------------
      Call schmid(lpsi, ncolx, ncoly, ncolz, npsi, npmin, spnorm, psi, tpsi, pswk1, pswk2, itimrev, wxyz)
!-----------------------------------------------------------------------
!        print initial quantities
!-----------------------------------------------------------------------
      If(iprint > 0 .And. mprint /= 0) Then
         Call kine(ncolx, ncoly, ncolz, lpsi, npsi, npmin, der2x, der2y, der2z, spkine, psi, pswk1, wxyz)
         Call energy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, der1y, &
        & der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, wyuk, edens, &
        & worka, icoul, ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, &
        & ehfbj2, ehfc, ecorc, ehfy, itheta, itimrev, iodds)
         Call isoenergy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, der1y, &
        & der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, wyuk, &
        & isoedens, worka, icoul, itheta, itimrev, iodds)
         Call sinfo(ncolx, ncoly, ncolz, lpsi, wx, wy, wz, xclx, xcly, xclz, der1x, der1y, der1z, h2m, t3, x3, alpha, itrs, npsi, &
        & npmin, vocc, spkine, spenrg, speflu, spnorm, delesum, spar, ajz, xcm, ycm, zcm, rho, psi, pswk1, pswk2, pswk3, ehft, &
        & ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ehf, ecorc, &
        & ecorr, ecorl, erear3, erearc, erearl, rlam, iconstr, ehfy, alz, asz, epairn, epairp, gapn, gapp, fermin, fermip,&
        & efluct, isoedens, etheta)
      End If
! if we are reading the stored g.s. wfs for dynamic initialization return here
      If(irwgs == 1) Then
         Return
      End If
! initialize delq and qerr so that non-constraint calculations do not accidentally converge
      delq = 1.0_wp
      qerr = 1.0_wp
      If(iconstr == 1 .And. irest == 0) Then
         xlam = 0.0_wp
      End If
! comment the next line to make above block active
      xlam = 0.0_wp ! setting xlam=0 seems better than starting from previous lambda
      xlamx = 0.0_wp
      qxin = 0.0_wp
!-----------------------------------------------------------------------
!        static iteration begins here
!-----------------------------------------------------------------------
      Do itrs = irest + 1, irest + itrbx
!
         If(iconstr == 1 .Or. iconstr == 3 .Or. iconstr == 4) Then
            damp = 0.0_wp
            Call distance(ncolx,ncoly,ncolz,xclx,xcly,xclz,wx,wy,wz,rho,damp)
! can set damp=1 to test constraints without any damping
!           damp = 1.0_wp
! call moments to get the principal axis idp
            Call moments(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, q20, q22, q20tot, q22tot, dxq, &
           & dzq, idp, etheta)
            Call qmomc(idp, ncolx, ncoly, ncolz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, rho, q20old, q202old)
            Call momentqx(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, qxold, qx2old, idp)
         End If
!
         ioprint = 0
         If(mprint /= 0) Then
            If(iprint > 0 .And. (Mod(itrs, mprint) == 0 .Or. itrs == irest+itrbx)) ioprint = 1
         End If
!
         delesum = 0.0_wp
         sumflu = 0.0_wp
!
! ** neutron states **
!
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(nst,nst2,spe,efluct,denerg,cnormb,pswk1,pswk2,pswk3,pswk4) REDUCTION(+:sumflu,delesum)
!$OMP DO
         Do nst = npmin(1), npsi(1)
            spe = spenrg(nst, 1)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
            Call grstep(ncolx, ncoly, ncolz, 1, mprint, iprint, itrs, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, &
           & dxiq, der1x, der1y, der1z, der2x, der2y, der2z, spe, efluct, denerg, pswk4, pswk1, pswk2, pswk3, cdmpx, cdmpy, &
           & cdmpz, x0act, itrbx, irest, wxyz, xcm, ycm, zcm, xclx, xcly, xclz, xlam, xlamx, damp, iconstr, idp, itimrev, cnormb)
            If(ihdiag == 1 .AND. mod(itrs,mhdiag) == 0) Then
               Do nst2 = npmin(1), npsi(1)
                  Call psnorm(ncolx, ncoly, ncolz, psi(:,:,:,:,nst2,1), pswk3, cnormb, 1, wxyz)
                  hmat(nst2, nst, 1) = cnormb
               End Do
               hmat(nst,nst,1) = cnormb + spe
            End If
            psi(:ncolx, :ncoly, :ncolz, :2, nst, 1) = pswk4(:ncolx, :ncoly, :ncolz, :2)
            spenrg(nst, 1) = spe
            speflu(nst, 1) = efluct
            sumflu = sumflu + vocc(nst, 1) * efluct
            delesum = delesum + vocc(nst, 1) * denerg
         End Do
!$OMP END DO NOWAIT
!
! ** proton states **
!
!$OMP DO
         Do nst = npmin(2), npsi(2)
            spe = spenrg(nst, 2)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
            Call grstep(ncolx, ncoly, ncolz, 2, mprint, iprint, itrs, upot, bmass, xiq, hsigma, cq, dcq, bmunu, dbmu, dbmass, &
           & dxiq, der1x, der1y, der1z, der2x, der2y, der2z, spe, efluct, denerg, pswk4, pswk1, pswk2, pswk3, cdmpx, cdmpy, &
           & cdmpz, x0act, itrbx, irest, wxyz, xcm, ycm, zcm, xclx, xcly, xclz, xlam, xlamx, damp, iconstr, idp, itimrev, cnormb)
            If(ihdiag == 1 .AND. mod(itrs,mhdiag) == 0) Then
               Do nst2 = npmin(2), npsi(2)
                  Call psnorm(ncolx, ncoly, ncolz, psi(:,:,:,:,nst2,2), pswk3, cnormb, 1, wxyz)
                  hmat(nst2, nst, 2) = cnormb
               End Do
               hmat(nst,nst,2) = cnormb + spe
            End If
            psi(:ncolx, :ncoly, :ncolz, :2, nst, 2) = pswk4(:ncolx, :ncoly, :ncolz, :2)
            spenrg(nst, 2) = spe
            speflu(nst, 2) = efluct
            sumflu = sumflu + vocc(nst, 2) * efluct
            delesum = delesum + vocc(nst, 2) * denerg
         End Do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!-----------------------------------------------------------------------
!     perform pairing calculation
!-----------------------------------------------------------------------
         If(imode == 0) Then
            If(nexadd == 1 .and. itrs <= nexiter) Then
               ireord = 1
            Else
               ireord = 0
            End If
            Call pairing(ipairn, ipairp, itrs, irest, ncolx, ncoly, ncolz, lpsi, npmin, npsi, nneut, nprot, epairn, epairp,&
                         gapn, gapp, fermin, fermip, spenrg, vocc, v2, wxyz, rho, v0prot, v0neut, rho0pr, psi, ireord, itimrev)
         End If
!-----------------------------------------------------------------------
!        perform a gram-schmidt orthogonalization of all states
!-----------------------------------------------------------------------
         If(mod(itrs,mshmid) == 0 .Or. itrs < 200) then
            Call schmid(lpsi, ncolx, ncoly, ncolz, npsi, npmin, spnorm, psi, tpsi, pswk1, pswk2, itimrev, wxyz)
         End If
!-----------------------------------------------------------------------
!        Diagonalize H-matrix
!-----------------------------------------------------------------------
         If(ihdiag == 1 .AND. mod(itrs,mhdiag) == 0) Then
            Call heev(A=hmat(npmin(1):npsi(1),npmin(1):npsi(1),1), W=work, uplo='L', jobz='V')
            Call heev(A=hmat(npmin(2):npsi(2),npmin(2):npsi(2),2), W=work, uplo='L', jobz='V')
            psi2 = psi
            psi(:,:,:,:,npmin(1):npsi(1),1) = (0.0_wp, 0.0_wp)
            psi(:,:,:,:,npmin(2):npsi(2),2) = (0.0_wp, 0.0_wp)
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(nst,nst2,iq) REDUCTION(+:psi)
!$OMP DO
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                  Do nst2 = npmin(iq), npsi(iq)
                     psi(:,:,:,:,nst,iq) = psi(:,:,:,:,nst,iq) + hmat(nst2,nst,iq)*psi2(:,:,:,:,nst2,iq)
                  End Do
               End Do
            End Do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
         End If
!-----------------------------------------------------------------------
!        program branches to do the density constraint step.
!        compute density, no need to call densit here.
!-----------------------------------------------------------------------
         If(iconstr == 2) Then
!
            rho = 0.0_wp
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                  Do is = 1, 2
                     Do iz = 1, ncolz
                        Do iy = 1, ncoly
                           Do ix = 1, ncolx
                              rho(ix, iy, iz, iq) = rho(ix, iy, iz, iq) + vocc(nst, iq) * psi(ix, iy, iz, is, nst, iq) * &
                             & conjg(psi(ix, iy, iz, is, nst, iq))
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
!
            rhonew = rho
            drho = rhonew - rhoold
            rlamold = rlam
            rlam = rlam + c0 * drho / (2.0_wp*x0dmp*rhoold+d0)
            drho = rhoold - rho0
            drlam = c0 * drho / (2.0_wp*x0dmp*rhoold+d0)
            rlamcof = rlam - rlamold + drlam
!
! ** neutron states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk4) SCHEDULE(STATIC)
            Do nst = npmin(1), npsi(1)
               pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
               Call costep2(1, ncolx, ncoly, ncolz, rlamcof, pswk4, pswk1)
               psi(:ncolx, :ncoly, :ncolz, :2, nst, 1) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1) - x0dmp * pswk1(:ncolx, :ncoly, &
              & :ncolz, :2)
            End Do
!$OMP END PARALLEL DO
!
! ** proton states **
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,pswk1,pswk4) SCHEDULE(STATIC)
            Do nst = npmin(2), npsi(2)
               pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
               Call costep2(2, ncolx, ncoly, ncolz, rlamcof, pswk4, pswk1)
               psi(:ncolx, :ncoly, :ncolz, :2, nst, 2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2) - x0dmp * pswk1(:ncolx, :ncoly, &
              & :ncolz, :2)
            End Do
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------
!        perform a gram-schmidt orthogonalization of all states
!-----------------------------------------------------------------------
            If(mod(itrs,mshmid) == 0 .Or. itrs < 200) then
               Call schmid(lpsi, ncolx, ncoly, ncolz, npsi, npmin, spnorm, psi, tpsi, pswk1, pswk2, itimrev, wxyz)
            End If
!-----------------------------------------------------------------------
!        end density constraint step
!-----------------------------------------------------------------------
         End If
!-----------------------------------------------------------------------
!        program branches to do the quadrupole constraint step.
!        compute the density, no need to call densit here.
!-----------------------------------------------------------------------
         If(iconstr == 1 .Or. iconstr == 3 .Or. iconstr == 4) Then
!
            rho = 0.0_wp
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                  Do is = 1, 2
                     Do iz = 1, ncolz
                        Do iy = 1, ncoly
                           Do ix = 1, ncolx
                              rho(ix, iy, iz, iq) = rho(ix, iy, iz, iq) + vocc(nst, iq) * psi(ix, iy, iz, is, nst, iq) * &
                             & conjg(psi(ix, iy, iz, is, nst, iq))
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
!
            Call qmomc(idp, ncolx, ncoly, ncolz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, rho, q20new, q202new)
            Call momentqx(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, qxnew, qx2new, idp)
!
            delq = q20new - q20old
            varqsq = ABS(q202old-q20old*q20old/(nneut+nprot))
!            dellam = c0*(q20old-q2in) / (2.0_wp*varqsq+d0)
            dellam = c0*(q20new-q2in) / (2.0_wp*varqsq+d0)
! other constraints (dimensions hidden)
            delqx = qxnew - qxold
            varqsqx = ABS(qx2old-qxold*qxold/(nneut+nprot))
            dellamx = c0*(qxold-qxin)/(2.0_wp*varqsqx+d0)
! quadrupole constraint
            xlam = xlam + q2c0*e0dmp/x0dmp * delq / (2.0_wp*varqsq+d0)
! other constraints (hidden array index)
            xlamx = xlamx + q2c0*e0dmp/x0dmp * delqx / (2.0_wp*varqsqx+d0)
!
            qerr = q20new - q2in
            If(abs(qerr) <= q2eps) ioprint = 1
            Write(8, '(/,A,I5)') 'iteration= ', itrs
            Write(8, '(A)') 'Constraints:     xlam          dellam          delq            qin         q2new         qerr'
            Write(8, '(1P,A,6E14.5)') "constraints: ",xlam, dellam, delq, q2in, q20new, qerr
!
            If(iconstr == 1) Then ! just the quadrupole constraint
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                     Call costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm, -dellam, damp, pswk4, &
                                 & psi(1, 1, 1, 1, nst, iq), idp)
                  End Do
               End Do
            ElseIf(iconstr == 3) Then ! quadrupole constraint + on <x>,<y>,<z>, and q21
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                     Call costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm, -dellam, damp, pswk4, &
                                & psi(1, 1, 1, 1, nst, iq), idp)
                  End Do
               End Do
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                    Call costep3(ncolx, ncoly, ncolz, xclx, xcly, xclz, -dellamx, damp, pswk4, &
                               & psi(1, 1, 1, 1, nst, iq), idp)
                  End Do
               End Do
            ElseIf(iconstr == 4) Then ! Everyhing in iconstr=3 + q22 axial constraint
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                     Call costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm,-dellam, damp, pswk4, &
                                & psi(1, 1, 1, 1, nst, iq), idp)
                  End Do
               End Do
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                    Call costep3(ncolx, ncoly, ncolz, xclx, xcly, xclz, -dellamx, damp, pswk4, &
                               & psi(1, 1, 1, 1, nst, iq), idp)
                  End Do
               End Do
               Do iq = 1, 2
                  Do nst = npmin(iq), npsi(iq)
                     pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
                     Call costep4(ncolx, ncoly, ncolz, xclx, xcly, xclz, -dellamx, damp, pswk4, &
                                & psi(1, 1, 1, 1, nst, iq), idp)
                  End Do
               End Do
            End If
!-----------------------------------------------------------------------
!        perform a gram-schmidt orthogonalization of all states
!-----------------------------------------------------------------------
            If(mod(itrs,mshmid) == 0 .Or. itrs < 200) then
               Call schmid(lpsi, ncolx, ncoly, ncolz, npsi, npmin, spnorm, psi, tpsi, pswk1, pswk2, itimrev, wxyz)
            End If
!-----------------------------------------------------------------------
!        end quadrupole constraint step
!-----------------------------------------------------------------------
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
! ** neutron states **
!
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(nst,pswk1,pswk4) REDUCTION(+:rho,tau,currnt,sodens,spinden,kinvden,spincur,rhos,taus,drhos,currnts)
!$OMP DO
         Do nst = npmin(1), npsi(1)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 1)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 1, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
!$OMP END DO NOWAIT
!
! ** proton states **
!
!$OMP DO
         Do nst = npmin(2), npsi(2)
            pswk4(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, 2)
            Call densit(lpsi, ncolx, ncoly, ncolz, nst, 2, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, &
           & der1y, der1z, pswk4, pswk1, rhos, taus, drhos, currnts, itimrev, itheta, iodds)
         End Do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
         Call skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x, &
        & der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, &
        & upot, bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1, &
        & coulplan2, rho2, q, iperc)
!
         If(iconstr == 2) Then
            rhoold = rho
            upot = upot + rlam
         End If
!-----------------------------------------------------------------------
!        for print modulus or last iteration calculate quantities
!-----------------------------------------------------------------------
         If(ioprint == 1) Then
            Call kine(ncolx, ncoly, ncolz, lpsi, npsi, npmin, der2x, der2y, der2z, spkine, psi, pswk1, wxyz)
            Call energy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, der1y, &
           & der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, wyuk, edens, &
           & worka, icoul, ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, &
           & ehfbj2, ehfc, ecorc, ehfy, itheta, itimrev, iodds)
           Call isoenergy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz,der1x,der1y,&
          & der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, wyuk, isoedens,&
          & worka, icoul, itheta, itimrev, iodds)
            Call sinfo(ncolx, ncoly, ncolz, lpsi, wx, wy, wz, xclx, xcly, xclz, der1x, der1y, der1z, h2m, t3, x3, alpha, itrs, &
           & npsi, npmin, vocc, spkine, spenrg, speflu, spnorm, delesum, spar, ajz, xcm, ycm, zcm, rho, psi, pswk1, pswk2, pswk3, &
           & ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, &
           & ehf, ecorc, ecorr, ecorl, erear3, erearc, erearl, rlam, iconstr, ehfy, alz, asz, epairn, epairp, gapn, gapp, fermin, &
           & fermip, efluct, isoedens, etheta)
            sumflu = sumflu / (nneut+nprot)
!-----------------------------------------------------------------------------
!        if converged or end of max iterations go out of static iteration loop
!-----------------------------------------------------------------------------
            If(((efluct < serr .And. itrs > irest+5) .Or. abs(qerr) <= q2eps) .Or. (itrs == irest+itrbx)) Then
              isoscalarenergy = sum(wxyz(:,:,:)*isoedens(:,:,:,1))
              isovectorenergy = sum(wxyz(:,:,:)*isoedens(:,:,:,2))

              Write(8,'(//,A,1P,E15.6,/,A,1P,E15.6,/,A,1P,E15.6,/,A,1P,E15.6)') &
              & ' Coulomb fenergy    (mev)= ', ehfc,&
              & ' Isovector fenergy  (mev)= ', isovectorenergy,&
              & ' Isoscalar fenergy  (mev)= ', isoscalarenergy, &
              & ' IsoTotal fenergy   (mev)= ', isovectorenergy+isoscalarenergy
              flush(8)
! these were for calculating the single-particle potentials to fusion barrier
!            write(99,'(12f10.4)') (spenrg(il,1), il = 1, lpsi)
!            write(99,'(12f10.4)') (spenrg(il,2), il = 1, lpsi)
!            write(99,'(12f10.4)') (h2m(1)*spkine(il,1), il = 1, lpsi)
!            write(99,'(12f10.4)') (h2m(2)*spkine(il,2), il = 1, lpsi)
!            write(99,'(12f10.4)') (erear3(il,1), il = 1, lpsi)
!            write(99,'(12f10.4)') (erear3(il,2), il = 1, lpsi)
!            write(99,'(12f10.4)') (erearc(il,1), il = 1, lpsi)
!            write(99,'(12f10.4)') (erearc(il,2), il = 1, lpsi)
!            write(99,'(12f10.4)') (erearl(il,1), il = 1, lpsi)
!            write(99,'(12f10.4)') (erearl(il,2), il = 1, lpsi)
! these are to write mesh and density when desired for special calculations
!             write(39) (xclx(ix), ix=1,ncolx)
!             write(39) (xcly(iy), iy=1,ncoly)
!             write(39) (xclz(iz), iz=1,ncolz)
!             write(39) (wx(ix), ix=1,ncolx)
!             write(39) (wy(iy), iy=1,ncoly)
!             write(39) (wz(iz), iz=1,ncolz)
!             write(39) ((((rho(ix,iy,iz,1)+rho(ix,iy,iz,2)),ix=1,ncolx),iy=1,ncoly),iz=1,ncolz)
!

               Call wrtpsi(12, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, &
              & spar, wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam)

               If(imode == 0 .And. irwgs == 0) Then
                  Call wrtpsi(15, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, &
                 & spar, wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam)
               End If
               Exit
            End If
!
         End If
!
         If((Mod(itrs, mrest) == 0)) Then
            Call wrtpsi(12, psi, ncolx, ncoly, ncolz, lpsi, npsi, npmin, vocc, xcm, ycm, zcm, spenrg, spnorm, spkine, ajz, spar, &
           & wcoul, epairn, epairp, gapn, gapp, fermin, fermip, xlam)
         End If
! don't plot if we are doing density constraint
         If(mplot /= 0) Then
            If((Mod(itrs, mplot) == 0) .And. (iconstr /= 2)) Then
               Allocate(rhoplt(ncolx, ncolz))
               iy = ncoly / 2
               If(mod(ncoly,2) /= 0) iy = iy + 1
               rhoplt(:ncolx, :ncolz) = rho(:ncolx, iy, :ncolz, 1) + rho(:ncolx, iy, :ncolz, 2)
               Call grid_plot(rhoplt, "Nuclear-Density", xclx, xclz, ncolx, ncolz, 13, 2)
               Deallocate(rhoplt)
            End If
         End If
      End Do
!-----------------------------------------------------------------------------
!        calcuate localization function
!-----------------------------------------------------------------------------
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
                     &          drhos(ix,iy,iz,3,is,iq,iquant)**2
                        cursq2 = currnts(ix,iy,iz,1,is,iq,iquant)**2 + currnts(ix,iy,iz,2,is,iq,iquant)**2 + &
                     &  currnts(ix,iy,iz,3,is,iq,iquant)**2
                        tautf = 3.0_wp/5.0_wp*(6.0_wp*pi**2)**(2.0_wp/3.0_wp)*rhos(ix,iy,iz,is,iq,iquant)**(5.0_wp/3.0_wp)
                        tautf_int(is,iq,iquant) = tautf_int(is,iq,iquant) + wxyz(ix,iy,iz)*(tautf+0.5_wp*grsq2)
                        taus_int(is,iq,iquant) = taus_int(is,iq,iquant) + wxyz(ix,iy,iz)*taus(ix,iy,iz,is,iq,iquant)
                        tausp_int(is,iq,iquant) = tausp_int(is,iq,iquant) + wxyz(ix,iy,iz)*( &
                     &  0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) + cursq2/rhos(ix,iy,iz,is,iq,iquant))
                        clocal(ix,iy,iz,is,iq,iquant) = taus(ix,iy,iz,is,iq,iquant) - &
                     &                                  0.25_wp*grsq2/rhos(ix,iy,iz,is,iq,iquant) - &
                     &                                  cursq2/rhos(ix,iy,iz,is,iq,iquant)
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
!               dqmu(iq,iquant) = dqmu(iq,iquant) + (tautf_int(is,iq,iquant)-tausp_int(is,iq,iquant))
            end do
            !write(*,*) dqmu(iq,iquant)
            write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-s: ", iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0(iq))
!            write(8,'(A,f10.2,2x,2i3,2f18.6)') "Pauli-s: ", roft, iq, iquant, h2m(iq)*dqmu(iq,iquant), h2m(iq)*(dqmu(iq,iquant)-dqmu0_tf(iq))
         end do
      end do
      write(8,*) "======================="
! always plot the final density
      If(iconstr /= 2) Then
         Allocate(rhoplt(ncolx, ncolz))
         iy = ncoly / 2
         If(mod(ncoly,2) /= 0) iy = iy + 1
         rhoplt(:ncolx, :ncolz) = rho(:ncolx, iy, :ncolz, 1) + rho(:ncolx, iy, :ncolz, 2)
         Call grid_plot(rhoplt, "Nuclear-Density", xclx, xclz, ncolx, ncolz, 13, 2)
         Deallocate(rhoplt)
      End If
!      allocate ( rhoplt(ncolx,ncolz) )
!      iy = ncoly/2
!      rhoplt(:ncolx,:ncolz) = clocal(:ncolx,iy,:ncolz,1,2)
!      call grid_plot(rhoplt,"Localization",xclx,xclz,ncolx,ncolz,55,2)
!      rhoplt(:ncolx,:ncolz) = clocal(:ncolx,iy,:ncolz,1,1)
!      call grid_plot(rhoplt,"Localization",xclx,xclz,ncolx,ncolz,56,2)
!      deallocate ( rhoplt )
!
      If(iconstr == 1) Then
         Write(8, '(2/,A,1P,E12.4)') " Quadrupole constraint lambda = ", xlam
      End If
!
      Deallocate(erear3, erearc, erearl)
      Deallocate( rhos, taus, clocal, drhos, currnts )
      If(iconstr == 2) Then
         Deallocate(rho0, rhoold, drho, rhonew, rlamold, drlam, rlamcof)
      Else if(iconstr == 1) Then
         Deallocate(damp)
      End If
!
      Return
End Subroutine shf3d
