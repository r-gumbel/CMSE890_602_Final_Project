Subroutine tinfo(ncolx, ncoly, ncolz, nx2, ny2, nz2, lpsi, mprint, mplot, iprint, it, wx, wy, wz, xclx, xcly, xclz, h2m, t3, x3, alpha, time, &
& npsi, npmin, vocc, spkine, spenrg, spnorm, spar, ajz, xcm, ycm, zcm, centf, roft, roftq20, rdot, tdot, rho, tau, currnt, ehft, ehf0, &
& ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, ehfy, mnof, nof, &
& ifixb, vx, vz, xb, tetrs_s, xdt, irest, psi, niter, nfixb, rold_s, tetold_s, rscale_s, q20old_s, q20tot_s, xdot_s, zdot_s, &
& ecoul, wcoul, estar, xvcmc, isoedens, der1x, der1y, der1z, coulplan1, coulplan2, rho2, q, iperc, fermin, fermip, epairn, epairp)
!
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: nx2
      Integer, Intent(In) :: ny2
      Integer, Intent(In) :: nz2
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: mprint
      Integer, Intent(In) :: mplot
      Integer, Intent(In) :: iprint
      Integer, Intent(In) :: it
      Integer, Intent(In) :: mnof
      Integer, Intent(In) :: nof
      Integer, Intent(In) :: ifixb
      Integer, Intent(In) :: niter
      Integer, Intent(In) :: nfixb
      Integer, Intent(In) :: irest
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Integer, Intent(In) :: iperc
      Integer(8), Intent(In) :: coulplan1, coulplan2
      Real(wp), Intent(In) :: vocc(lpsi, 2)
      Real(wp), Intent(In) :: t3
      Real(wp), Intent(In) :: x3
      Real(wp), Intent(In) :: alpha
      Real(wp), Intent(In) :: time
      Real(wp), Intent(In) :: ehft
      Real(wp), Intent(In) :: ehf0
      Real(wp), Intent(In) :: ehf12
      Real(wp), Intent(In) :: ehf22
      Real(wp), Intent(In) :: ehf3
      Real(wp), Intent(In) :: ehfls
      Real(wp), Intent(In) :: ehf0_odd
      Real(wp), Intent(In) :: ehf12_odd
      Real(wp), Intent(In) :: ehf22_odd
      Real(wp), Intent(In) :: ehf3_odd
      Real(wp), Intent(In) :: ehfls_odd
      Real(wp), Intent(In) :: ehfbj2
      Real(wp), Intent(In) :: ehflj2
      Real(wp), Intent(In) :: ehfc
      Real(wp), Intent(In) :: ehfy
      Real(wp), Intent(In) :: ecorc
      Real(wp), Intent(In) :: xdt
      Real(wp), Intent(In) :: xb
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: h2m(2)
      Real(wp), Intent(In) :: spkine(lpsi, 2)
      Real(wp), Intent(In) :: spenrg(lpsi, 2)
      Real(wp), Intent(In) :: spnorm(lpsi, 2)
      Real(wp), Intent(In) :: spar(lpsi, 2)
      Real(wp), Intent(In) :: ajz(lpsi, 2)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: tau(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: ecoul(mnof)
      Real(wp), Intent(In) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: estar(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: isoedens(ncolx,ncoly,ncolz,2)
      Real(wp), Intent(Inout) :: xcm(2)
      Real(wp), Intent(Inout) :: ycm(2)
      Real(wp), Intent(Inout) :: zcm(2)
      Real(wp), Intent(Inout) :: roft
      Real(wp), Intent(Inout) :: rdot
      Real(wp), Intent(Inout) :: tdot
      Real(wp), Intent(Inout) :: centf(3, mnof)
      Real(wp), Intent(Inout) :: vx
      Real(wp), Intent(Inout) :: vz
      Real(wp), Intent(Inout) :: tetrs_s
      Real(wp), Intent(Inout) :: rold_s
      Real(wp), Intent(Inout) :: tetold_s
      Real(wp), Intent(Inout) :: rscale_s
      Real(wp), Intent(Inout) :: q20old_s
      Real(wp), Intent(Inout) :: q20tot_s
      Real(wp), Intent(Inout) :: xdot_s
      Real(wp), Intent(Inout) :: zdot_s
      Real(wp), Intent(Inout) :: xvcmc
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp)             :: fermin(2), fermip(2), epairn(2), epairp(2)
      Complex(wp) :: q(nx2, ny2, nz2), rho2(nx2, ny2, nz2)
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: hbc = 197.3269885804381_wp
      Real(wp), Parameter :: e2 = 1.439978408596513_wp
      Real(wp), Parameter :: small = 1.0e-25_wp
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iff, iq, iz, iy, ix, iyy, nps, nst, is, npsm
      Integer :: jfirst, lfirst, ipfirst, itcm, ilm, ilx, idp
      Integer :: il, ixcut1, ixcut2, izcut1, izcut2, locx1(1), locx2(1), locz1(1), locz2(1), ifrag, frag(ncolx,ncoly,ncolz)
      Real(wp) :: x1old, x2old, z1old, z2old, x1c, x2c, z1c, z2c, f1c, f2c, pnrq, xcmq, ycmq, zcmq, z, y, x, vol, pnrn, zcmn, &
     & ycmn, xcmn, pnrp, pnrtot1, zcmp, ycmp, xcmp, pnrtot, xcmtot, xcmtot_old, ycmtot, zcmtot, amc2, slope, bb, diff, xmin, &
     & zmin, angle, slopev, bbv, xlmass, xrmass, xlchrg, xrchrg, rhotot, rhopp, xmu, ratio, vx1, vx2, vz1, vz2, tlix, tliz, trix, &
     & triz, tkel, tker, tketot, tetr, tdotc, xlf, xcoul, vcoul, xcent, ecmf, epsf, ttt, tetc, tetf, tets, rmsq, x2mom, y2mom, &
     & z2mom, xsq, ysq, zsq, rmsn, x2mn, y2mn, z2mn, rmsp, x2mp, y2mp, z2mp, rmstot, zcmtot1, rmstot1, x2mtot1, y2mtot1, z2mtot1, &
     & xcmtot1, ycmtot1, x2mtot, y2mtot, z2mtot, ecorr, rhop, rhon, ehf, tke, ehfint, ehfsum, tkerr, xperi, f5, centxl, centxr, &
     & centyl, centyr, centzl, centzr, ehf_even, ehf_odd, tetrs, estarl, estarr, estart, y1old, y2old, weight, x1, &
     & isoscalarenergy, isovectorenergy, momentum(3,0:2), vel(3,0:2), angmom(3,0:2), xcoul2, epairhf, efermi(2), &
     & ctime, ctimetmp
      Real(wp) :: centx(nof), centy(nof), centz(nof), xmass(nof)
      Real(wp) :: centln, centlp, centrn, centrp, ecoultot, ecoul_1, ecoul_2, ekinl, ekinr
      Real(wp) :: roftq20, rold, roldq20, tetold, rscale, q20old, q20(2), q22(2), dxq, dzq
      Real(wp) :: q30q, q40q, q30n, q30p, q40n, q40p
      Real(wp) :: q10xq, q10yq, q10zq, q10xn, q10xp, q10yn, q10yp, q10zn, q10zp
      Real(wp) :: q10totx, q10toty, q10totz, q10tot1x, q10tot1y, q10tot1z
      Real(wp) :: q20tot, q20tot1, q22tot, q22tot1, q30tot, q30tot1, q40tot, q40tot1
      Real(wp) :: xdot, zdot, rdotq20, ucmc, vcmc, cur2q, cur2t, erotx, eroty, erotz
      Real(wp) :: e_coll(2), e_kin(2), e_collt, e_colltt, rhoam1, slate, ecortot, ehfkoop
      Real(wp) :: xlfxl, xlfxr, xlfyl, xlfyr, xlfzl, xlfzr, xlftot, xjltot, vy1, vy2
      Real(wp) :: curtotx, curtoty, curtotz
      Real(wp) :: xi(3), zi(3,3), etheta(3)
      Real(wp) :: rsx, rsy, rsz, sigis
      integer :: ic
      Real(wp), Dimension(:, :), Allocatable :: erear3, erearc
      Real(wp), Dimension(:, :), Allocatable :: rhoplt
      Real(wp), Dimension(:), Allocatable :: xwk
      Real(wp), Dimension(:), Allocatable :: zwk
      Real(wp) :: rhow(ncolx,ncoly,ncolz,mnof), wcoul_1(ncolx,ncoly,ncolz), wcoul_2(ncolx,ncoly,ncolz)
!
      Complex(wp) :: psinq
      complex(wp) :: pswk1(ncolx,ncoly,ncolz,2), pswk2(ncolx,ncoly,ncolz,2), pswk3(ncolx,ncoly,ncolz,2)
!
      Save rold, x1old, x2old, y1old, y2old, z1old, z2old, tetold, tetrs, xmin, zmin, xcmtot_old
      Save x1c, x2c, z1c, z2c, f1c, f2c, ixcut1, ixcut2, izcut1, izcut2, ctime, ctimetmp
      Save jfirst, lfirst, ipfirst, rscale, q20old, q20tot, xdot, zdot
!
      Data jfirst / 0 /
      Data lfirst / 0 /
      Data ipfirst / 0 /
!
      Allocate(erear3(lpsi, 2), erearc(lpsi, 2))
      Allocate(rhoplt(ncolx, ncolz), xwk(ncolx+1), zwk(ncolz+1))
! initialize in case of restart
      If((nof == 2 .And. ifixb == 1) .And. (irest /= 0 .And. jfirst == 0)) Then
         jfirst = 1
         tetrs = tetrs_s
         rold = rold_s
         roldq20 = rold_s
         tetold = tetold_s
         angle = tetold
         rscale = rscale_s
         q20old = q20old_s
         q20tot = q20tot_s
         xdot = xdot_s
         zdot = zdot_s
      End If
!-----------------------------------------------------------------------
!        calculate mass, center of mass, and charge coordinates
!-----------------------------------------------------------------------
      Do iq = 1, 2
!
         pnrq = 0.0_wp
         xcmq = 0.0_wp
         ycmq = 0.0_wp
         zcmq = 0.0_wp
         Do iz = 1, ncolz
            z = xclz(iz)
            Do iy = 1, ncoly
               y = xcly(iy)
               Do ix = 1, ncolx
                  x = xclx(ix)
                  vol = wx(ix)*wy(iy)*wz(iz)*rho(ix, iy, iz, iq)
                  pnrq = vol + pnrq
                  xcmq = vol*x + xcmq
                  ycmq = vol*y + ycmq
                  zcmq = vol*z + zcmq
               End Do
            End Do
         End Do
         If(iq == 1) Then
            pnrn = pnrq
            zcmn = zcmq
            ycmn = ycmq
            xcmn = xcmq
         Else
            pnrp = pnrq
            zcmp = zcmq
            ycmp = ycmq
            xcmp = xcmq
         End If
      End Do
!
      pnrtot = pnrn + pnrp
      xcmtot = (xcmn + xcmp) / pnrtot
      ycmtot = (ycmn + ycmp) / pnrtot
      zcmtot = (zcmn + zcmp) / pnrtot
      xcmtot1 = (xcmn - xcmp) / pnrtot
      ycmtot1 = (ycmn - ycmp) / pnrtot
      zcmtot1 = (zcmn - zcmp) / pnrtot
      xcm(1) = xcmn / pnrn
      xcm(2) = xcmp / pnrp
      ycm(1) = ycmn / pnrn
      ycm(2) = ycmp / pnrp
      zcm(1) = zcmn / pnrn
      zcm(2) = zcmp / pnrp
! quadrupole moments (also used for calculating R)
      Call moments(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, q20, q22, q20tot, q22tot, dxq, dzq, &
     & idp, etheta)
!  initialize R stuff at start once
      If(xdt < 0.0_wp .And. nof == 2 .And. ifixb == 1 .And. jfirst == 0 .And. niter /= 0) Then
         jfirst = 1
         ctime = 0.0
         ctimetmp = 0.0
         x1old = centf(1, 1)
         x2old = centf(1, 2)
         y1old = centf(2, 1)
         y2old = centf(2, 2)
         z1old = centf(3, 1)
         z2old = centf(3, 2)
         rold = sqrt((x2old - x1old)**2 + (z2old - z1old)**2)
         roldq20 = roft
         rscale = roldq20 / Sqrt(Abs(q20tot))
         q20old = q20tot
         xcmtot_old = xcmtot
         vcmc = 0.0_wp
         Call getslope(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcmtot, ycmtot, zcmtot, rho, slope)
         tetold = Atan(slope)
!  calculate the correct angle in different quadrants using the eigenvector of quadrupole tensor
         If(tetold < 0.0_wp .And. abs(slope) > 1.0e-3_wp) Then
            tetold = pi + tetold
         Else If(centf(1,1) > 0.0_wp .And. centf(3,1) > 0.0_wp .And. abs(slope) > 1.0e-3_wp) Then
            tetold = pi + tetold
         Else If(abs(slope) < 1.0e-3_wp) Then
            tetold = 0.0_wp
         End If
         tetrs = tetold
      End If
      If(xdt < 0.0_wp .And. nof <= 2 .And. ifixb == 0) Then
         x1old = centf(1, 1)
         x2old = centf(1, 2)
         y1old = centf(2, 1)
         y2old = centf(2, 2)
         z1old = centf(3, 1)
         z2old = centf(3, 2)
         rold = Sqrt((x2old - x1old)**2 + (z2old - z1old)**2)
      End If
!-----------------------------------------------------------------------
!     calculation of the dividing plane. specialize to x-z reaction
!     plane, which assumes that fragment centers stay around y=0.
!     note: fragment 2 should be elevated along positive z-axis for
!     non-central collisions.
!
!     nof=1 case corresponds to a fissioning nucleus that may end up in
!     two fragments (mnof should be at least 2 here).
!-----------------------------------------------------------------------
      amc2 = hbc*hbc / (h2m(1) + h2m(2))
!
      If(nof == 1) Then
         If(mnof < 2) Then
            Write(8, '(//,a,//)') "Set mnof to at least 2 in case of fission"
            Stop
         End If
         centf = 0.0_wp
         xlmass = 0.0_wp
         xrmass = 0.0_wp
         xlchrg = 0.0_wp
         xrchrg = 0.0_wp
         centln = 0.0_wp
         centlp = 0.0_wp
         centrn = 0.0_wp
         centrp = 0.0_wp
!
         zmin = - 1.0_wp
!
         Do ix = 1, ncolx
            x = xclx(ix)
            Do iy = 1, ncoly
               y = xcly(iy)
               Do iz = 1, ncolz
                  z = xclz(iz)
                  If(z > zmin) Then
                     vol = wx(ix)*wy(iy)*wz(iz)
                     rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                     rhopp = rho(ix, iy, iz, 2)
                     xlmass = xlmass + vol*rhotot
                     xlchrg = xlchrg + vol*rhopp
                     centln = centln + vol*rho(ix, iy, iz, 1)*x
                     centlp = centlp + vol*rho(ix, iy, iz, 2)*x
                     centf(1, 1) = centf(1, 1) + vol*rhotot*x
                     centf(2, 1) = centf(2, 1) + vol*rhotot*y
                     centf(3, 1) = centf(3, 1) + vol*rhotot*z
                  Else
                     vol = wx(ix)*wy(iy)*wz(iz)
                     rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                     rhopp = rho(ix, iy, iz, 2)
                     xrmass = xrmass + vol*rhotot
                     xrchrg = xrchrg + vol*rhopp
                     centrn = centrn + vol*rho(ix, iy, iz, 1)*x
                     centrp = centrp + vol*rho(ix, iy, iz, 2)*x
                     centf(1, 2) = centf(1, 2) + vol*rhotot*x
                     centf(2, 2) = centf(2, 2) + vol*rhotot*y
                     centf(3, 2) = centf(3, 2) + vol*rhotot*z
                  End If
               End Do
            End Do
         End Do
         centln = centln / (xlmass-xlchrg)
         centlp = centlp / xlchrg
         centrn = centrn / (xrmass-xrchrg)
         centrp = centrp / xrchrg
         centf(1, 1) = centf(1, 1) / xlmass
         centf(2, 1) = centf(2, 1) / xlmass
         centf(3, 1) = centf(3, 1) / xlmass
         centf(1, 2) = centf(1, 2) / xrmass
         centf(2, 2) = centf(2, 2) / xrmass
         centf(3, 2) = centf(3, 2) / xrmass
         centx(1) = centf(1, 1)
         centx(2) = centf(1, 2)
         centz(1) = centf(3, 1)
         centz(2) = centf(3, 2)
         xmass(1) = xlmass
         xmass(2) = xrmass
         If(xdt > 1.0e-5_wp) Then
            vx1 = (centf(1, 1) - x1old) / xdt
            vx2 = (centf(1, 2) - x2old) / xdt
            vy1 = (centf(2, 1) - y1old) / xdt
            vy2 = (centf(2, 2) - y2old) / xdt
            vz1 = (centf(3, 1) - z1old) / xdt
            vz2 = (centf(3, 2) - z2old) / xdt
            tlix = 0.5_wp*xlmass*amc2*vx1**2
            tliz = 0.5_wp*xlmass*amc2*vz1**2
            trix = 0.5_wp*xrmass*amc2*vx2**2
            triz = 0.5_wp*xrmass*amc2*vz2**2
            tkel = tlix + tliz
            tker = trix + triz
         End If
         x1old = centf(1, 1)
         x2old = centf(1, 2)
         y1old = centf(2, 1)
         y2old = centf(2, 2)
         z1old = centf(3, 1)
         z2old = centf(3, 2)
!
      Else If(nof == 2 .And. ifixb == 1) Then
!
         Do itcm = 1, 10
            If(itcm > 1) Then
               centx = centf(1, 1:2)
               centz = centf(3, 1:2)
               slope = (centz(2) - centz(1)) / (centx(2) - centx(1) + small)
               bb = centz(1) - slope*centx(1)
            Else
               centx = centf(1, 1:2)
               centz = centf(3, 1:2)
! first time: get direction of axis of largest quadrupole moment to calculate intercept bb
               Call getslope(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcmtot, ycmtot, zcmtot, rho, slope)
               bb = zcmtot - slope*xcmtot
            End If
!
            If(lfirst == 0) Then
               locx1 = minloc(Abs(xclx - centf(1, 1)))
               locx2 = minloc(Abs(xclx - centf(1, 2)))
               locz1 = minloc(Abs(xclz - centf(3, 1)))
               locz2 = minloc(Abs(xclz - centf(3, 2)))
               ixcut1 = locx1(1)
               ixcut2 = locx2(1)
               izcut1 = locz1(1)
               izcut2 = locz2(1)
            End If
            If(niter <= nfixb) Then
               xmin = (centx(1) + centx(2))/2.0_wp
               zmin = (centz(1) + centz(2))/2.0_wp
            Else
               Call divpoint(ncolx, ncoly, ncolz, xclx, xclz, slope, bb, ixcut1, ixcut2, izcut1, izcut2, rho, &
           &                 xmin, zmin, ctimetmp)
            End If
            angle = Atan(slope)
            slopev = Tan(angle + pi/2.0_wp)
! if you want to manually set these do it here (for defining dividing plane and restarting).
!         xmin = -9.0_wp
!         zmin = 0.0_wp
!        centx(1) = 0.5_wp
!        centx(2) = 0.0_wp
!        centz(1) = -15.0_wp
!        centz(2) = 5.0_wp
            bbv = zmin - slopev*xmin
            rhow = 0.0_wp
            centf = 0.0_wp
            xlmass = 0.0_wp
            xrmass = 0.0_wp
            xlchrg = 0.0_wp
            xrchrg = 0.0_wp
            centln = 0.0_wp
            centlp = 0.0_wp
            centrn = 0.0_wp
            centrp = 0.0_wp
            estarl = 0.0_wp
            estarr = 0.0_wp
            estart = 0.0_wp
            ekinl = 0.0_wp
            ekinr = 0.0_wp
            momentum = 0.0_wp
            Do iz = 1, ncolz
               z = xclz(iz)
               Do ix = 1, ncolx
                  x = xclx(ix)
                  diff = z - slopev*x - bbv
                  If(abs(slope) <= 1.0e-3_wp) diff = x - xmin
                  If(diff <= 0.0_wp) Then
                     Do iy = 1, ncoly
                        y = xcly(iy)
                        vol = wx(ix)*wy(iy)*wz(iz)
                        rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                        rhow(ix,iy,iz,1) = rho(ix,iy,iz,2)
                        rhopp = rho(ix, iy, iz, 2)
                        curtotx = currnt(ix,iy,iz,1,1) + currnt(ix,iy,iz,1,2)
                        curtoty = currnt(ix,iy,iz,2,1) + currnt(ix,iy,iz,2,2)
                        curtotz = currnt(ix,iy,iz,3,1) + currnt(ix,iy,iz,3,2)
                        momentum(:,1) = momentum(:,1) + currnt(ix,iy,iz,:,1) + currnt(ix,iy,iz,:,2)
                        momentum(:,0) = momentum(:,0) + currnt(ix,iy,iz,:,1) + currnt(ix,iy,iz,:,2)
                        xlmass = xlmass + vol*rhotot
                        xlchrg = xlchrg + vol*rhopp
                        centln = centln + vol*rho(ix, iy, iz, 1)*x
                        centlp = centlp + vol*rho(ix, iy, iz, 2)*x
                        centf(1, 1) = centf(1, 1) + vol*rhotot*x
                        centf(2, 1) = centf(2, 1) + vol*rhotot*y
                        centf(3, 1) = centf(3, 1) + vol*rhotot*z
                        ekinl = (h2m(1) + h2m(2))/2.0_wp*(curtotx**2 + curtoty**2 + curtotz**2)/(rhotot + small)
                        estarl = estarl + vol*(estar(ix, iy, iz) - ekinl)
                        frag(ix,iy,iz) = 1
                     End Do
                  Else
                     Do iy = 1, ncoly
                        y = xcly(iy)
                        vol = wx(ix)*wy(iy)*wz(iz)
                        rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                        rhow(ix,iy,iz,2) = rho(ix,iy,iz,2)
                        rhopp = rho(ix, iy, iz, 2)
                        curtotx = currnt(ix,iy,iz,1,1) + currnt(ix,iy,iz,1,2)
                        curtoty = currnt(ix,iy,iz,2,1) + currnt(ix,iy,iz,2,2)
                        curtotz = currnt(ix,iy,iz,3,1) + currnt(ix,iy,iz,3,2)
                        momentum(:,2) = momentum(:,2) + currnt(ix,iy,iz,:,1) + currnt(ix,iy,iz,:,2)
                        momentum(:,0) = momentum(:,0) + currnt(ix,iy,iz,:,1) + currnt(ix,iy,iz,:,2)
                        xrmass = xrmass + vol*rhotot
                        xrchrg = xrchrg + vol*rhopp
                        centrn = centrn + vol*rho(ix, iy, iz, 1)*x
                        centrp = centrp + vol*rho(ix, iy, iz, 2)*x
                        centf(1, 2) = centf(1, 2) + vol*rhotot*x
                        centf(2, 2) = centf(2, 2) + vol*rhotot*y
                        centf(3, 2) = centf(3, 2) + vol*rhotot*z
                        ekinr = (h2m(1) + h2m(2))/2.0_wp*(curtotx**2 + curtoty**2 + curtotz**2)/(rhotot + small)
                        estarr = estarr + vol*(estar(ix, iy, iz) - ekinr)
                        frag(ix,iy,iz) = 2
                     End Do
                  End If
               End Do
            End Do
            estart = estarl + estarr
            centln = centln / (xlmass - xlchrg + small)
            centlp = centlp / (xlchrg + small)
            centrn = centrn / (xrmass - xrchrg + small)
            centrp = centrp / (xrchrg + small)
            centf(1, 1) = centf(1, 1) / (xlmass + small)
            centf(2, 1) = centf(2, 1) / (xlmass + small)
            centf(3, 1) = centf(3, 1) / (xlmass + small)
            centf(1, 2) = centf(1, 2) / (xrmass + small)
            centf(2, 2) = centf(2, 2) / (xrmass + small)
            centf(3, 2) = centf(3, 2) / (xrmass + small)
            xmass(1) = xlmass
            xmass(2) = xrmass
            vel(:,1) = momentum(:,1) / xmass(1)
            vel(:,2) = momentum(:,2) / xmass(2)
            vel(:,0) = momentum(:,0) / (xmass(1) + xmass(2))
            angmom = 0.0_wp
            Do ifrag = 1, 2
               Do iz = 1, ncolz
                  z = xclz(iz) - centf(3, ifrag)
                  Do iy = 1, ncoly
                     y = xcly(iy) - centf(2, ifrag)
                     Do ix = 1, ncolx
                        x = xclx(ix) - centf(1, ifrag)
                        If(frag(ix,iy,iz) /= ifrag) Cycle
                        vol = wx(ix)*wy(iy)*wz(iz)
                        rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                        curtotx = currnt(ix,iy,iz,1,1) + currnt(ix,iy,iz,1,2) - rhotot*vel(1,ifrag)
                        curtoty = currnt(ix,iy,iz,2,1) + currnt(ix,iy,iz,2,2) - rhotot*vel(2,ifrag)
                        curtotz = currnt(ix,iy,iz,3,1) + currnt(ix,iy,iz,3,2) - rhotot*vel(3,ifrag)
                        angmom(1,ifrag) = angmom(1,ifrag) + vol*(y*curtotz-z*curtoty)
                        angmom(2,ifrag) = angmom(2,ifrag) - vol*(x*curtotz-z*curtotx)
                        angmom(3,ifrag) = angmom(3,ifrag) + vol*(x*curtoty-y*curtotx)
                     End Do
                  End Do
               End Do
            End Do
            Do iz = 1, ncolz
               z = xclz(iz) - zcmtot
               Do iy = 1, ncoly
                  y = xcly(iy) - ycmtot
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcmtot
                     vol = wx(ix)*wy(iy)*wz(iz)
                     rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                     curtotx = currnt(ix,iy,iz,1,1) + currnt(ix,iy,iz,1,2) - rhotot*vel(1,0)
                     curtoty = currnt(ix,iy,iz,2,1) + currnt(ix,iy,iz,2,2) - rhotot*vel(2,0)
                     curtotz = currnt(ix,iy,iz,3,1) + currnt(ix,iy,iz,3,2) - rhotot*vel(3,0)
                     angmom(1,0) = angmom(1,0) + vol*(y*curtotz-z*curtoty)
                     angmom(2,0) = angmom(2,0) - vol*(x*curtotz-z*curtotx)
                     angmom(3,0) = angmom(3,0) + vol*(x*curtoty-y*curtotx)
                  End Do
               End Do
            End Do
            If(Abs(Max(MAXVAL(centx - centf(1, 1:2)), MAXVAL(centz - centf(3, 1:2)))) < 1.0e-05_wp) Exit
         End Do
         If(lfirst == 0) Then
            lfirst = 1
            x1c = centf(1, 1)
            x2c = centf(1, 2)
            z1c = centf(3, 1)
            z2c = centf(3, 2)
            f1c = xmass(2) / (xmass(1) + xmass(2))
            f2c = xmass(1) / (xmass(1) + xmass(2))
         End If
!
! calculate moments
         Call moments(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, q20, q22, q20tot, q22tot, dxq, dzq, &
        & idp, etheta)
! do these for full time-step
         If(xdt > 1.0e-5_wp) Then
            vx1 = (centf(1, 1) - x1old) / xdt
            If(Abs(vx1) > 1.0_wp) vx1 = 0.0_wp
            vx2 = (centf(1, 2) - x2old) / xdt
            If(Abs(vx2) > 1.0_wp) vx2 = 0.0_wp
            vy1 = (centf(2, 1) - y1old) / xdt
            If(Abs(vy1) > 1.0_wp) vy1 = 0.0_wp
            vy2 = (centf(2, 2) - y2old) / xdt
            If(Abs(vy2) > 1.0_wp) vy2 = 0.0_wp
            vz1 = (centf(3, 1) - z1old) / xdt
            If(Abs(vz1) > 1.0_wp) vz1 = 0.0_wp
            vz2 = (centf(3, 2) - z2old) / xdt
            If(Abs(vz2) > 1.0_wp) vz2 = 0.0_wp
            tlix = 0.5_wp*xlmass*amc2*vx1**2
            tliz = 0.5_wp*xlmass*amc2*vz1**2
            trix = 0.5_wp*xrmass*amc2*vx2**2
            triz = 0.5_wp*xrmass*amc2*vz2**2
            tkel = tlix + tliz
            tker = trix + triz
            roftq20 = rscale*Sqrt(Abs(q20tot))
            roft = Sqrt((centf(1, 2) - centf(1, 1))**2 + (centf(3, 2) - centf(3, 1))**2)
            If(niter > nfixb) Then
               roldq20 = rscale*Sqrt(Abs(q20old))
               q20old = q20tot
               x1old = centf(1, 1)
               x2old = centf(1, 2)
               y1old = centf(2, 1)
               y2old = centf(2, 2)
               z1old = centf(3, 1)
               z2old = centf(3, 2)
            End If
            rdot = (roft - rold)/xdt
            rdotq20 = (roftq20 - roldq20)/xdt
            If(Abs(rdot) > 1.0_wp) rdot = 0.0_wp
            If(niter > nfixb) rold = roft
!-------------------------------------------------------------------------------------
! calculate the velocity of the c.m., if moving, correct with opposite boost
! we average the velocity every mprint iterations
!-------------------------------------------------------------------------------------
            ucmc = (xcmtot - xcmtot_old) / xdt
            vcmc = vcmc + ucmc
            xvcmc = 0.0_wp
            If(Mod(it, mprint) == 0) Then
               xvcmc = vcmc / mprint
               vcmc = 0.0_wp
            End If
            xcmtot_old = xcmtot
!
            tetr = angle
!-----------------------------------------------------------------------------
!  correct tetr for different quadrants
!-----------------------------------------------------------------------------
            If(tetr < 0.0_wp .And. abs(slope) > 1.0e-3_wp) Then
               tetr = pi + tetr
            Else If(centf(1,1) > 0.0_wp .And. centf(3,1) > 0.0_wp .And. abs(slope) > 1.0e-3_wp) Then
               tetr = pi + tetr
            Else If(abs(slope) < 1.0e-3_wp .And. rdot < 0.0_wp) Then
               tetr = 0.0_wp
               tetc = pi - tetr
            Else If(abs(slope) < 1.0e-3_wp .And. rdot > 0.0_wp) Then
               tetr = pi
               tetc = pi - tetr
            End If
!-----------------------------------------------------------------------------
!  this sums over the angle since we can have multiple rotations.
!-----------------------------------------------------------------------------
            If(niter > nfixb) Then
               tetrs = tetrs + (tetr - tetold)
            End If
            tdotc = (tetr - tetold) / xdt
            If(niter <= nfixb) Then
               tdotc = tdot
               tetold = tetr
            End If
            If(time == xdt) tdotc = tdot
            If(Abs(tdotc) > pi/180.0_wp) tdotc = 0.0_wp
! calculate the actual coulomb interaction potential using asymptotic values in the
! entrance channel. in exit channel use the point coulomb formula (not good for deformed systems).
            Call pois(iperc, coulplan1, coulplan2, ncolx, ncoly, ncolz, nx2, ny2, nz2, wx(1), wy(1), wz(1),&
                   rhow(1, 1, 1, 1), rho2, q, wcoul_1)
            Call pois(iperc, coulplan1, coulplan2, ncolx, ncoly, ncolz, nx2, ny2, nz2, wx(1), wy(1), wz(1),&
                   rhow(1, 1, 1, 2), rho2, q, wcoul_2)
            ecoultot = 0.0_wp
            ecoul_1 = 0.0_wp
            ecoul_2 = 0.0_wp
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     ecoultot = ecoultot + 0.5_wp*wx(ix)*wy(iy)*wz(iz)*wcoul(ix, iy, iz)*rho(ix, iy, iz, 2)
                     ecoul_1  = ecoul_1  + 0.5_wp*wx(ix)*wy(iy)*wz(iz)*wcoul_1(ix, iy, iz)*rhow(ix, iy, iz, 1)
                     ecoul_2  = ecoul_2  + 0.5_wp*wx(ix)*wy(iy)*wz(iz)*wcoul_2(ix, iy, iz)*rhow(ix, iy, iz, 2)
                  End Do
               End Do
            End Do
            xcoul = ecoultot - ecoul_1 - ecoul_2
            xcoul2 = ecoultot - ecoul(1) - ecoul(2)
            vcoul = xlchrg*xrchrg*e2 / (roft + small)
            Do iq = 1, 2
               e_coll(iq) = 0.0_wp
               Do iz = 1, ncolz
                  Do iy = 1, ncoly
                     Do ix = 1, ncolx
                        vol = wx(ix)*wy(iy)*wz(iz)
                        cur2q = currnt(ix, iy, iz, 1, iq)**2 + currnt(ix, iy, iz, 2, iq)**2 + currnt(ix, iy, iz, 3, iq)**2
                        e_coll(iq) = e_coll(iq) + h2m(iq)*vol*cur2q / (rho(ix, iy, iz, iq) + small)
                     End Do
                  End Do
               End Do
            End Do
            e_collt = e_coll(1) + e_coll(2)
!
            centxl = centf(1, 1)
            centxr = centf(1, 2)
            centyl = centf(2, 1)
            centyr = centf(2, 2)
            centzl = centf(3, 1)
            centzr = centf(3, 2)
! cross-product for orbital angular momentum
            xlfxl = (centyl*momentum(3,1)-centzl*momentum(2,1))
            xlfxr = (centyr*momentum(3,2)-centzr*momentum(2,2))
            xlfyl = (centzl*momentum(1,1)-centxl*momentum(3,1)) ! minus sign used here
            xlfyr = (centzr*momentum(1,2)-centxr*momentum(3,2)) ! minus sign used here
            xlfzl = (centxl*momentum(2,1)-centyl*momentum(1,1))
            xlfzr = (centxr*momentum(2,2)-centyr*momentum(1,2))
            xlftot = sqrt((xlfxl+xlfxr)**2+(xlfyl+xlfyr)**2+(xlfzl+xlfzr)**2)
!
            xmu = xlmass*xrmass / (xlmass + xrmass + small)*amc2
            ratio = xmu / hbc**2
            tketot = 0.5_wp*xmu*rdot**2
!            xlf = tdotc*xmu / hbc*roft**2
            xlf = xlftot
            xlf = abs(xlf)
            xcent = xlf**2 / (2.0_wp*ratio*roft**2 + small)
            If(rdot <= 0.0_wp) Then
               ecmf = e_collt + xcoul
            Else
               ecmf = tketot + xcoul + xcent
            End If
            If(ecmf < 0.0_wp) ecmf = 0.0_wp
            epsf = Sqrt(1.0_wp + 2.0_wp*ecmf*xlf**2/(ratio*(xlchrg*xrchrg*e2)**2 + small))
            ttt = xlf**2 / (ratio*xlchrg*xrchrg*e2*roft)
            tetf = Acos(Min(1.0_wp/(epsf + small), 1.0_wp)) - Acos(Min(1.0_wp/(epsf + small)*(1.0_wp + ttt), 1.0_wp))
            xlf = tdotc*xmu / hbc*roft**2
!-----------------------------------------------------------------------------
!  calculate the correct angle in different quadrants using the information
!-----------------------------------------------------------------------------
            if(xlmass > xrmass) x1 = centf(1, 2)
            if(xlmass <= xrmass) x1 = centf(1, 1)
            if(angle < 0.0_wp .and. x1 >= 0.0) then
               tetc = -angle
               tets = tetc - tetf
            else if(angle < 0.0_wp .and. x1 <= 0.0 .and. xb /= 0.0) then
               tetc = -(pi + angle)
               tets = tetc - tetf
            else if(angle < 0.0_wp .and. x1 <= 0.0 .and. xb == 0.0) then
               tetc = -(pi + angle)
               tets = tetc + tetf
            else if(angle > 0.0_wp .and. x1 >= 0.0) then
               tetc = -angle
               tets = tetc - tetf
            else if(angle > 0.0_wp .and. x1 <= 0.0 .and. xb /= 0.0) then
               tetc = pi - angle
               tets = tetc - tetf
            else if(angle > 0.0_wp .and. x1 <= 0.0 .and. xb == 0.0) then
               tetc = pi - angle
               tets = tetc - tetf
            end if
!
            vx = rdot*Cos(tetr) - roft*Sin(tetr)*tdotc
            xdot = vx
            vz = rdot*Sin(tetr) + roft*Cos(tetr)*tdotc
            zdot = vz
            x1c = x1c - f1c*xdot*xdt
            x2c = x2c + f2c*xdot*xdt
            z1c = z1c - f1c*zdot*xdt
            z2c = z2c + f2c*zdot*xdt
            locx1 = minloc(Abs(xclx - x1c))
            locx2 = minloc(Abs(xclx - x2c))
            locz1 = minloc(Abs(xclz - z1c))
            locz2 = minloc(Abs(xclz - z2c))
            ixcut1 = locx1(1)
            ixcut2 = locx2(1)
            izcut1 = locz1(1)
            izcut2 = locz2(1)
            If(niter > nfixb) Then
               tetold = tetr
               tetrs_s = tetrs
               rold_s = rold
               tetold_s = tetold
               rscale_s = rscale
               q20old_s = q20old
               q20tot_s = q20tot
               xdot_s = xdot
               zdot_s = zdot
            End If
         End If
!
      Else If(nof >= 2 .And. ifixb == 0) Then
         Do iff = 1, nof
            centx(iff) = centf(1, iff)
            centy(iff) = centf(2, iff)
            centz(iff) = centf(3, iff)
         End Do
!
         If(nof == 2) Then
!-----------------------------------------------------------------------
!     in case of 2 fragments but in LAB frame we can still try to define
!     a L-R dividing planes. Repeat the relevant chunks from above that
!     are not specific to c.m. frame.
!-----------------------------------------------------------------------
            xmin = (centx(1)+centx(2)) / 2.0_wp
            zmin = (centz(1)+centz(2)) / 2.0_wp
!
            slope = (centz(2) - centz(1)) / (centx(2) - centx(1) + small)
            angle = Atan(slope)
            slopev = Tan(angle + pi/2.0_wp)
            bbv = zmin - slopev*xmin
!
            centf = 0.0_wp
            xlmass = 0.0_wp
            xrmass = 0.0_wp
            xlchrg = 0.0_wp
            xrchrg = 0.0_wp
            centln = 0.0_wp
            centlp = 0.0_wp
            centrn = 0.0_wp
            centrp = 0.0_wp
            Do iz = 1, ncolz
               z = xclz(iz)
               Do ix = 1, ncolx
                  x = xclx(ix)
                  diff = z - slopev*x - bbv
                  If(diff <= 0.0_wp) Then
                     Do iy = 1, ncoly
                        vol = wx(ix)*wy(iy)*wz(iz)
                        rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                        rhopp = rho(ix, iy, iz, 2)
                        xlmass = xlmass + vol*rhotot
                        xlchrg = xlchrg + vol*rhopp
                        centln = centln + vol*rho(ix, iy, iz, 1)*x
                        centlp = centlp + vol*rho(ix, iy, iz, 2)*x
                        centf(1, 1) = centf(1, 1) + vol*rhotot*x
                        centf(2, 1) = centf(2, 1) + vol*rhotot*y
                        centf(3, 1) = centf(3, 1) + vol*rhotot*z
                     End Do
                  Else
                     Do iy = 1, ncoly
                        vol = wx(ix)*wy(iy)*wz(iz)
                        rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                        rhopp = rho(ix, iy, iz, 2)
                        xrmass = xrmass + vol*rhotot
                        xrchrg = xrchrg + vol*rhopp
                        centrn = centrn + vol*rho(ix, iy, iz, 1)*x
                        centrp = centrp + vol*rho(ix, iy, iz, 2)*x
                        centf(1, 2) = centf(1, 2) + vol*rhotot*x
                        centf(2, 2) = centf(2, 2) + vol*rhotot*y
                        centf(3, 2) = centf(3, 2) + vol*rhotot*z
                     End Do
                  End If
               End Do
            End Do
!
            centln = centln / (xlmass - xlchrg + small)
            centlp = centlp / (xlchrg + small)
            centrn = centrn / (xrmass - xrchrg + small)
            centrp = centrp / (xrchrg + small)
            centf(1, 1) = centf(1, 1) / (xlmass + small)
            centf(2, 1) = centf(2, 1) / (xlmass + small)
            centf(3, 1) = centf(3, 1) / (xlmass + small)
            centf(1, 2) = centf(1, 2) / (xrmass + small)
            centf(2, 2) = centf(2, 2) / (xrmass + small)
            centf(3, 2) = centf(3, 2) / (xrmass + small)
            centx(1) = centf(1, 1)
            centx(2) = centf(1, 2)
            centz(1) = centf(3, 1)
            centz(2) = centf(3, 2)
!
            If(xdt > 1.0e-5_wp) Then
               vx1 = (centf(1, 1) - x1old) / xdt
               vx2 = (centf(1, 2) - x2old) / xdt
               vy1 = (centf(2, 1) - y1old) / xdt
               vy2 = (centf(2, 2) - y2old) / xdt
               vz1 = (centf(3, 1) - z1old) / xdt
               vz2 = (centf(3, 2) - z2old) / xdt
!
               roft = Sqrt((centf(1, 2) - centf(1, 1))**2 + (centf(3, 2) - centf(3, 1))**2)
               rdot = (roft - rold) / xdt
            End If
         End If
!
      End If
!-----------------------------------------------------------------------
!        calculate radii, moments, etc. (in c.m. frame)
!-----------------------------------------------------------------------
! Accumulate ctime
      If(xdt>1e-5) Call divpoint(ncolx, ncoly, ncolz, xclx, xclz, slope, bb, ixcut1, ixcut2, izcut1, izcut2, rho, &
      &            xmin, zmin, ctime)
      If(mprint /= 0) Then
         If(iprint > 0 .And. (Mod(it,mprint) == 0 .Or. it == 1)) Then
!
            Call moments(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, q20, q22, q20tot, q22tot, &
                         dxq, dzq, idp, etheta)
            Call inertia(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcmtot, ycmtot, zcmtot, rho, xi, zi)
!
            Do iq = 1, 2
               rmsq = 0.0_wp
               q30q = 0.0_wp
               q40q = 0.0_wp
               x2mom = 0.0_wp
               y2mom = 0.0_wp
               z2mom = 0.0_wp
               Do iz = 1, ncolz
                  z = xclz(iz) - zcm(iq)
                  zsq = z*z
                  Do iy = 1, ncoly
                     y = xcly(iy) - ycm(iq)
                     ysq = y*y
                     Do ix = 1, ncolx
                        x = xclx(ix) - xcm(iq)
                        xsq = x*x
                        vol = wx(ix)*wy(iy)*wz(iz)*rho(ix, iy, iz, iq)
                        rmsq = vol*(xsq + ysq + zsq) + rmsq
                        q30q = vol*(xsq*x - 1.5_wp*x*(zsq + ysq)) + q30q
                        q40q = vol*(xsq*xsq - 3.0_wp*xsq*(zsq + ysq) + 0.375_wp*(zsq + ysq)**2) + q40q
                        x2mom = vol*xsq + x2mom
                        y2mom = vol*ysq + y2mom
                        z2mom = vol*zsq + z2mom
                     End Do
                  End Do
               End Do
               q10xq = xcm(iq)
               q10yq = ycm(iq)
               q10zq = zcm(iq)
   !
               If(iq == 1) Then
                  rmsn = rmsq
                  q10xn = q10xq
                  q10yn = q10yq
                  q10zn = q10zq
                  q30n = q30q
                  q40n = q40q
                  x2mn = x2mom
                  y2mn = y2mom
                  z2mn = z2mom
               Else
                  rmsp = rmsq
                  q10xp = q10xq
                  q10yp = q10yq
                  q10zp = q10zq
                  q30p = q30q
                  q40p = q40q
                  x2mp = x2mom
                  y2mp = y2mom
                  z2mp = z2mom
               End If
            End Do
            q10xn = Sqrt(3.0_wp/(4.0_wp*pi))*q10xn
            q10xp = Sqrt(3.0_wp/(4.0_wp*pi))*q10xp
            q10yn = Sqrt(3.0_wp/(4.0_wp*pi))*q10yn
            q10yp = Sqrt(3.0_wp/(4.0_wp*pi))*q10yp
            q10zn = Sqrt(3.0_wp/(4.0_wp*pi))*q10zn
            q10zp = Sqrt(3.0_wp/(4.0_wp*pi))*q10zp
            q30n = Sqrt(7.0_wp/(4.0_wp*pi))*q30n
            q30p = Sqrt(7.0_wp/(4.0_wp*pi))*q30p
            q40n = Sqrt(9.0_wp/(4.0_wp*pi))*q40n
            q40p = Sqrt(9.0_wp/(4.0_wp*pi))*q40p
   !
            pnrtot1 = pnrn - pnrp
            rmstot = (rmsn+rmsp) / pnrtot
            rmstot1 = (rmsn-rmsp) / pnrtot
            rmstot = Sqrt(rmstot)
            rmstot1 = Sqrt(Abs(rmstot1))
            rmsp = Sqrt(rmsp/pnrp)
            rmsn = Sqrt(rmsn/pnrn)
            q10totx = q10xp + q10xn
            q10toty = q10yp + q10yn
            q10totz = q10zp + q10zn
            q10tot1x = q10xp - q10xn
            q10tot1y = q10yp - q10yn
            q10tot1z = q10zp - q10zn
            q20tot1 = q20(1) - q20(2)
            q22tot1 = q22(1) - q22(2)
            q30tot = q30n + q30p
            q30tot1 = q30n - q30p
            q40tot = q40n + q40p
            q40tot1 = q40n - q40p
            x2mtot = (x2mn+x2mp) / pnrtot
            y2mtot = (y2mn+y2mp) / pnrtot
            z2mtot = (z2mn+z2mp) / pnrtot
            x2mtot1 = (x2mn-x2mp) / pnrtot
            y2mtot1 = (y2mn-y2mp) / pnrtot
            z2mtot1 = (z2mn-z2mp) / pnrtot
            x2mn = x2mn / pnrn
            x2mp = x2mp / pnrp
            y2mn = y2mn / pnrn
            y2mp = y2mp / pnrp
            z2mn = z2mn / pnrn
            z2mp = z2mp / pnrp
!-----------------------------------------------------------------------
!     calculate the correction term to the hf energy due to the 3-body
!     force. this could be replaced but the answer is general for the
!     forms of skyrme force in use
!-----------------------------------------------------------------------
            ecorr = 0.0_wp
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     vol = wx(ix)*wy(iy)*wz(iz)
                     rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                     rhon = rho(ix, iy, iz, 1)
                     rhop = rho(ix, iy, iz, 2)
                     ecorr = ecorr + vol*rhotot**alpha*((1.0_wp + x3/2.0_wp)*rhotot**2 - (0.50_wp + x3)*(rhop**2 + rhon**2))
                  End Do
               End Do
            End Do
            ecorr = - alpha*t3 / 24.0_wp*ecorr
! consruct the correction terms for each single-particle state
            slate = - 3.0_wp / 4.0_wp*(3.0_wp/pi)**(1.0_wp/3.0_wp)*e2
            erear3 = 0.0_wp
            erearc = 0.0_wp
            ecortot = 0.0_wp
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
!                  if(spenrg(nst,iq) > 0.0_wp) exit
                  Do is = 1, 2
                     Do ix = 1, ncolx
                        Do iy = 1, ncoly
                           Do iz = 1, ncolz
                              vol = wx(ix)*wy(iy)*wz(iz)
                              rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                              rhon = rho(ix, iy, iz, 1)
                              rhop = rho(ix, iy, iz, 2)
                              rhoam1 = rhotot**(alpha-1.0_wp+1.0e-25_wp)
                              psinq = psi(ix, iy, iz, is, nst, iq)
                              erear3(nst, iq) = erear3(nst, iq) - alpha*t3 / 24.0_wp*vol*vocc(nst, iq)*conjg(psinq)*psinq &
                             &*rhoam1*((1.0_wp+x3/2.0_wp)*rhotot**2-(0.50_wp+x3)*(rhop**2+rhon**2))
                              erearc(nst, iq) = erearc(nst, iq) + vol*slate / 3.0_wp*vocc(nst, iq)*conjg(psinq)*psinq*rhop &
                             &**(1.0_wp/3.0_wp)
                           End Do
                        End Do
                     End Do
                  End Do
                  ecortot = ecortot + erear3(nst, iq)
                  If(iq == 2) ecortot = ecortot + erearc(nst, iq)
               End Do
            End Do
            erearc(:, 1) = 0.0_wp ! slater correction only for protons
!-----------------------------------------------------------------------
!     instead of brute force integration of the energy density we
!     calculate the hf energy from the single particle energies, the
!     kinetic energies, and the correction term
!
!     note: there are two ways to calculate the hf energy: 1) using
!           the koopman + correction terms. 2) integrating the energy
!           density. the integration of tau and the sum of the sp
!           kinetic energies give different results for coarse grids.
!           they approach one another for finer meshes or higher order
!           splines. the koopman expression also has a correction term
!           for the coulomb exchange force. with this koopman gives
!           same result if we use the sum over sp kinetic energies as
!           the total kinetic energy in both methods.
!-----------------------------------------------------------------------
            ehf = 0.0_wp
            tke = 0.0_wp
            Do iq = 1, 2
               nps = npsi(iq)
               npsm = npmin(iq)
               Do il = npsm, nps
!                  if(spenrg(il,iq) > 0.0_wp) cycle
                  ehf = ehf + vocc(il, iq)*h2m(iq)*spkine(il, iq) + vocc(il, iq)*spenrg(il, iq)
                  tke = tke + vocc(il, iq)*h2m(iq)*spkine(il, iq)
               End Do
            End Do
            ehf = ehf / 2.0_wp
            ehfkoop = ehf + ecortot
! same as ehfkoop but correction terms calculated by integration
            ehf = ehf + ecorr + ecorc
            ehfint = ehft + ehf0 + ehf12 + ehf22 + ehf3 + ehfls + ehfc + ehfy + ehf0_odd + ehf12_odd + ehf22_odd + ehf3_odd + &
           & ehfls_odd + ehflj2 + ehfbj2
            ehfsum = tke + ehf0 + ehf12 + ehf22 + ehf3 + ehfls + ehfc + ehfy + ehf0_odd + ehf12_odd + ehf22_odd + ehf3_odd + &
           & ehfls_odd + ehflj2 + ehfbj2
            ehf_odd = ehf0_odd + ehf12_odd + ehf22_odd + ehf3_odd + ehfls_odd + ehflj2
            ehf_even = ehfint - ehf_odd
            tkerr = ehft - tke
            epairhf = epairn(2) + epairp(2)
            ehfsum = ehfsum - epairhf
            ehfint = ehfint - epairhf
            ehf = ehf - epairhf
            efermi(1) = fermin(1)
            efermi(2) = fermip(1)
!-----------------------------------------------------------------------
!        calculate components of the excitation energy
!-----------------------------------------------------------------------
            erotx = 0.0_wp
            eroty = 0.0_wp
            erotz = 0.0_wp
            Do iq = 1, 2
               e_coll(iq) = 0.0_wp
               e_kin(iq) = 0.0_wp
               Do iz = 1, ncolz
                  Do iy = 1, ncoly
                     Do ix = 1, ncolx
                        vol = wx(ix)*wy(iy)*wz(iz)
                        cur2q = currnt(ix, iy, iz, 1, iq)**2 + currnt(ix, iy, iz, 2, iq)**2 + currnt(ix, iy, iz, 3, iq)**2
                        erotx = erotx + h2m(iq)*vol*currnt(ix, iy, iz, 1, iq)**2 / (rho(ix, iy, iz, iq) + small)
                        eroty = eroty + h2m(iq)*vol*currnt(ix, iy, iz, 2, iq)**2 / (rho(ix, iy, iz, iq) + small)
                        erotz = erotz + h2m(iq)*vol*currnt(ix, iy, iz, 3, iq)**2 / (rho(ix, iy, iz, iq) + small)
                        e_coll(iq) = e_coll(iq) + h2m(iq)*vol*cur2q / (rho(ix, iy, iz, iq) + small)
                        e_kin(iq) = e_kin(iq) + h2m(iq)*vol*tau(ix, iy, iz, iq)
                     End Do
                  End Do
               End Do
            End Do
            e_collt = e_coll(1) + e_coll(2)
            e_colltt = 0.0_wp
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     vol = wx(ix)*wy(iy)*wz(iz)
                     cur2t = (currnt(ix, iy, iz, 1, 1) + currnt(ix, iy, iz, 1, 2))**2 + (currnt(ix, iy, iz, 2, 1)&
                           + currnt(ix, iy, iz, 2, 2))**2 + (currnt(ix, iy, iz, 3, 1) + currnt(ix, iy, iz, 3, 2))**2
                     rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                     e_colltt = e_colltt + (h2m(1) + h2m(2)) / 2.0_wp*vol*cur2t / (rhotot + small)
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        calculate spin vector
!-----------------------------------------------------------------------
            rsz = 0.0_wp
            rsy = 0.0_wp
            rsx = 0.0_wp
            Do iq = 1, 2
               Do nst = npmin(iq), npsi(iq)
                  pswk3(:, :, :, :) = psi(:, :, :, :, nst, iq)
                  Call cmulz(ncolx, ncoly, ncolz, der1z, pswk3, pswk1, 0)
                  Call cmuly(ncolx, ncoly, ncolz, der1y, pswk3, pswk2, 0)
                  Do is = 1, 2
                     ic = 3 - is
                     sigis = 3 - 2 * is
                     Do iz = 1, ncolz
                        z = xclz(iz) - zcm(iq)
                        Do iy = 1, ncoly
                           y = xcly(iy) - ycm(iq)
                           Do ix = 1, ncolx
                              x = xclx(ix) - xcm(iq)
                              weight = wx(ix) * wy(iy) * wz(iz)
                              rsx = rsx + weight * conjg(pswk3(ix, iy, iz, is)) * pswk3(ix, iy, iz, ic)
                           End Do
                        End Do
                     End Do
                  End Do
                  Call cmulx(ncolx, ncoly, ncolz, der1x, pswk3, pswk1, 0)
                  Call cmulz(ncolx, ncoly, ncolz, der1z, pswk3, pswk2, 0)
                  Do is = 1, 2
                     ic = 3 - is
                     sigis = 3 - 2 * is
                     Do iz = 1, ncolz
                        z = xclz(iz) - zcm(iq)
                        Do iy = 1, ncoly
                           y = xcly(iy) - ycm(iq)
                           Do ix = 1, ncolx
                              x = xclx(ix) - xcm(iq)
                              weight = wx(ix) * wy(iy) * wz(iz)
                              rsy = rsy + 2.0_wp * sigis * weight * aimag(conjg(pswk3(ix, iy, iz, is))*pswk3(ix, iy, iz, ic))
                           End Do
                        End Do
                     End Do
                  End Do
                  Call cmulx(ncolx, ncoly, ncolz, der1x, pswk3, pswk1, 0)
                  Call cmuly(ncolx, ncoly, ncolz, der1y, pswk3, pswk2, 0)
                  Do is = 1, 2
                     ic = 3 - is
                     sigis = 3 - 2 * is
                     Do iz = 1, ncolz
                        z = xclz(iz) - zcm(iq)
                        Do iy = 1, ncoly
                           y = xcly(iy) - ycm(iq)
                           Do ix = 1, ncolx
                              x = xclx(ix) - xcm(iq)
                              weight = wx(ix) * wy(iy) * wz(iz)
                              rsz = rsz + sigis * weight * conjg(pswk3(ix, iy, iz, is)) * pswk3(ix, iy, iz, is)
                           End Do
                        End Do
                     End Do
                  End Do
!
               End Do
            End Do
            xjltot = sqrt((angmom(1,0)+0.5*rsx)**2+(angmom(2,0)+0.5*rsy)**2+(angmom(3,0)+0.5*rsz)**2)
!
            If(iprint > 0) Then
               If(ipfirst == 0) Then
                  Write(78, '(a)') "     time, pnrtot, rmstot, q10totx, q10toty, q10totz, q20tot, q22tot, q30tot, q40tot"
                  Write(79, '(a)') "     time, pnrtot1, rmstot1, q10x, q10y, q10z, q20tot1, q22tot1, q30tot1, q40tot1"
                  Write(80, '(a)') "     time, q21R, q21I, q22R, q22I, q31R, q31I, q32R, q32I, q33R, q33I"
                  ipfirst = 1
               End If
               Write(78, '(3F10.4,7E12.4)') time, pnrtot, rmstot, q10totx, q10toty, q10totz, q20tot, q22tot, q30tot, q40tot
               Write(79, '(3F10.4,7E12.4)') time, pnrtot1, rmstot1, q10tot1x, q10tot1y, q10tot1z, q20tot1, q22tot1, q30tot1, q40tot1
               Call moment_c(time, xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, idp)
            End If
!
!         If(iprint > 0 .And. (Mod(it, mprint) == 0 .Or. it == 1)) Then
            Write(8, '(3/,A,F12.4,/,A,F12.4,/,A,F12.4,/,A,F12.4)') ' time       (fm/c) = ', time, ' energy     (mev)  = ', &
           & ehfint, ' t.k.e.     (mev)  = ', ehft, ' tke error  (mev)  = ', tkerr
!-----------------------------------------------------------------------
!        plot the nuclear density
!-----------------------------------------------------------------------
            If(mplot /= 0) Then
               If(Mod(it, mplot) == 0 .Or. it == 1) Then
                  iyy = ncoly / 2
                  rhoplt(:ncolx, :ncolz) = rho(:ncolx, iyy, :ncolz, 2) + rho(:ncolx, iyy, :ncolz, 1)
                  xperi = (xclx(ncolx) - xclx(1)) / 12.0_wp
                  f5 = 0.140_wp
                  Call pagplt(rhoplt, xclx, ncolx, xclz, ncolz, f5, xperi, xwk, zwk)
                  Call grid_plot(rhoplt, "Nuclear-Density", xclx, xclz, ncolx, ncolz, 13, 2)
               End If
            End If
!
            If(nof == 2 .And. ifixb == 1) Then
!
               Write(8, '(2/,A,/,4(A,F12.4,A,F12.4,A,F12.4,A,F12.4,/),4(A,F12.4))') ' relative motion/coulomb kinematics:', '    &
              &  red. mass=', xmu / amc2, '          l/hbar=', xlf, '      ecmf(mev)=', ecmf, '        roft(fm)=', roft, '       &
              &  rdot/c=', rdot, '       td/c(deg)=', tdotc*180.0_wp / pi, '      trke(mev)=', tketot, '       tetr(deg)=', &
              & tetr*180.0_wp / pi, '     tetrs(deg)=', tetrs*180.0_wp / pi, '       tetf(deg)=', tetf*180.0_wp / pi, '    &
              &  tetc(deg)=', tetc*180.0_wp / pi, '       tets(deg)=', tets*180.0_wp / pi, '  vcoul(f)(mev)=', xcoul, '   vco&
              &ul(p)(mev)=', vcoul, '     vcent(mev)=', xcent, '        xcut(fm)=', xmin, '       zcut(fm)=', zmin, &
              & '      roftq2(fm)=', roftq20, '       rdotq2/c=', rdotq20, '        e*(mev)=', max(0.0_wp,estart)
!
               Write(8, '(2/,A,/,2A)') ' collision kinematics:', '                    mass       charge      <x>   ', '     <y>  &
              &       <z>         <vx>        <vz>        tke         E*'
               Write(8, '(A,2F12.4,1P,7E12.4)') ' fragment 1: ', xlmass, xlchrg, centxl, centyl, centzl, vx1, vz1, tkel, &
               max(0.0_wp,estarl)
               Write(8, '(A,2F12.4,1P,7E12.4)') ' fragment 2: ', xrmass, xrchrg, centxr, centyr, centzr, vx2, vz2, tker, &
               max(0.0_wp,estarr)
               Write(8, '(/,A,F12.4)') ' ctime = ', ctime*xdt
               Write(8, '(/,A,F12.4)') ' ctimetmp = ', ctimetmp*xdt
            End If
            If(nof >= 1 .And. ifixb == 0) Then
               If(nof == 1) Then
                  centxl = centf(1, 1)
                  centxr = centf(1, 2)
                  centyl = centf(2, 1)
                  centyr = centf(2, 2)
                  centzl = centf(3, 1)
                  centzr = centf(3, 2)
                  Write(8, '(2/,A,/,2A)') ' collision kinematics:', '                    mass       charge      <x>   ', '     <y&
                 &>         <z>         <vx>        <vz>        tke'
                  Write(8, '(A,2F12.4,1P,6E12.4)') ' fragment 1: ', xlmass, xlchrg, centxl, centyl, centzl, vx1, vz1, tkel
                  Write(8, '(A,2F12.4,1P,6E12.4)') ' fragment 2: ', xrmass, xrchrg, centxr, centyr, centzr, vx2, vz2, tker
               Else If(nof > 2) Then
                  Write(8, '(2/,A,/,2A)') ' fragment properties:', '                    mass       <x>   ', '     <y>         <z>&
                 &'
                  Do iff = 1, nof
                     Write(8, '(A,i2,A,F12.4,1P,3E12.4)') ' fragment ', iff, ': ', xmass(iff), centf(1, iff), centf(2, iff), &
                    & centf(3, iff)
                  End Do
               Else If(nof == 2) Then
                  centxl = centf(1, 1)
                  centxr = centf(1, 2)
                  centyl = centf(2, 1)
                  centyr = centf(2, 2)
                  centzl = centf(3, 1)
                  centzr = centf(3, 2)
                  Write(8, '(2/,A,/,2A)') ' collision kinematics:', '                    mass       charge      <x>   ', '     <y&
                 &>         <z>         <vx>        <vz>    '
                  Write(8, '(A,2F12.4,1P,6E12.4)') ' fragment 1: ', xlmass, xlchrg, centxl, centyl, centzl, vx1, vz1
                  Write(8, '(A,2F12.4,1P,6E12.4)') ' fragment 2: ', xrmass, xrchrg, centxr, centyr, centzr, vx2, vz2
!
                  Write(8, '(/,A,F12.4,A,F12.4)') ' roft = ', roft, '    rdot = ', rdot
                  Write(8, '(/,A,F12.4)') ' ctimetmp = ', ctimetmp*xdt
                  Write(8, '(/,A,F12.4)') ' ctime = ', ctime*xdt
               End If
            End If
!
            isoscalarenergy = 0.0_wp
            isovectorenergy = 0.0_wp
            Do ix = 1, ncolx
              Do iy = 1, ncoly
                Do iz = 1, ncolz
                   weight = wx(ix)*wy(iy)*wz(iz)
                   isoscalarenergy = isoscalarenergy + weight*isoedens(ix,iy,iz,1)
                   isovectorenergy = isovectorenergy + weight*isoedens(ix,iy,iz,2)
                End Do
              End Do
            End Do
            Write(8, '(/,A)') '                             total             time-even           time-odd'
            Write(8, '(A,1P,E15.6,4X,E15.6,4X,E15.6,/,A,E15.6,4X,E15.6,4X,E15.6,/,A,   E15.6,4X,E15.6,4X,E15.6,/,A,E15.6,4X,E15.6&
           &,4X,E15.6,/,A,   E15.6,4X,E15.6,4X,E15.6,/,A,   E15.6,23X,E15.6,/,A,   E15.6,4X,E15.6,/,A,E15.6,/,A,E15.6,/,A,E15.6,/,A,E15.6)'&
           & ) ' intg. energy    (mev)= ', ehfint, ehf_even, ehf_odd, ' t0 contribution (mev)= ', ehf0 + ehf0_odd, ehf0, &
           & ehf0_odd, ' t1,t2 terms     (mev)= ', ehf12 + ehf22 + ehf12_odd + ehf22_odd + ehflj2, ehf12 + ehf22, ehf12_odd + &
           & ehf22_odd + ehflj2, ' t3 contribution (mev)= ', ehf3 + ehf3_odd, ehf3, ehf3_odd, ' t4 contribution (mev)= ', ehfls + &
           & ehfls_odd, ehfls, ehfls_odd, ' j**2 term       (mev)= ', ehflj2, ehflj2, ' J**2 term       (mev)= ', ehfbj2, ehfbj2, &
           & ' yukawa energy   (mev)= ', ehfy, ' coulomb energy  (mev)= ', ehfc, ' Koopman energy  (mev)= ', ehfkoop, &
           & ' pairing energy  (mev)= ', -epairhf
            Write(8,'(/,A,1P,E15.6,/,A,1P,E15.6,/,A,1P,E15.6)') ' Isovector tenergy  (mev)= ', isovectorenergy, &
           & ' Isoscalar tenergy  (mev)= ', isoscalarenergy, &
           & ' IsoTotal tenergy   (mev)= ', isovectorenergy+isoscalarenergy
            Write(8, '(2/,A,/,3A)') ' excitation energies:', '                 e_coll     e_colltt     e_kin       erotx', '     &
           &  eroty       erotz         E*'
            Write(8, '(A,7F12.4)') ' total(ex):', e_collt, e_colltt, e_kin(1) + e_kin(2), erotx, eroty, erotz, max(0.0_wp,estart)
            Write(8, '(A,F12.4,12X,F12.4)') ' neutrons: ', e_coll(1), e_kin(1)
            Write(8, '(A,F12.4,12X,F12.4)') ' protons:  ', e_coll(2), e_kin(2)
!
            Write(8, '(2/,A,/,3A)') ' orbital angular momenta and spin:',&
            &'                    total         x           y           z'
            Write(8, '(A,4F12.4)') ' Relative L/hbar:',  xlftot, xlfxl + xlfxr, xlfyl + xlfyr, xlfzl + xlfzr
            Write(8, '(A,4F12.4)') ' Left L/hbar:    ', sqrt(sum(angmom(:,1)**2)), angmom(1,1), angmom(2,1), angmom(3,1)
            Write(8, '(A,4F12.4)') ' Right L/hbar:   ', sqrt(sum(angmom(:,2)**2)), angmom(1,2), angmom(2,2), angmom(3,2)
            Write(8, '(A,4F12.4)') ' Total L/hbar:   ', sqrt(sum(angmom(:,0)**2)), angmom(1,0), angmom(2,0), angmom(3,0)
            Write(8, '(A,4F12.4)') ' s/hbar:         ', 0.5_wp*sqrt(rsx**2+rsy**2+rsz**2), 0.5_wp*rsx, 0.5_wp*rsy, 0.5_wp*rsz
            Write(8, '(A,F12.4)') ' J/hbar:         ', xjltot
            Write(8, '(/,A,/,A)') ' inertia/m', '              I(1)        I(2)        I(3)'
            Write(8, '(A,3F12.4,/)')  '         ', xi(1), xi(2), xi(3)
            Write(8,*) 'Cartesian inertia principal values and axes:'
            Write(8,'(3(F10.2,''('',3F8.4,'')''))') (xi(ix),zi(:,ix),ix=1,3)
!-----------------------------------------------------------------------
!     the following block writes the single-particle quantities
!-----------------------------------------------------------------------
            Write(8, '(/,12x,a,29x,a)') '  neutron single particle states', ' proton single particle states'
            Write(8, '(  12x,a,29x,a)') '  ------------------------------', ' -----------------------------'
            Write(8, '(a,5x,a)') '  nr    jz      par   weight   norm     ekin     energy', '  nr    jz      par   weight   norm &
           &    ekin     energy'
            ilm = Min(npsi(1), npsi(2))
            ilx = Max(npsi(1), npsi(2))
            Do il = 1, ilm
               Write(8, '(1x,i3,4(f8.4),2(f10.4),5x,i3,4(f8.4),2(f10.4))') il, ajz(il, 1), spar(il, 1), vocc(il, 1), spnorm(il, &
              & 1), h2m(1)*spkine(il, 1), spenrg(il, 1), il, ajz(il, 2), spar(il, 2), vocc(il, 2), spnorm(il, 2), h2m(2)*&
              & spkine(il, 2), spenrg(il, 2)
            End Do
            Do il = ilm + 1, ilx
               If(npsi(1) > npsi(2)) Then
                  Write(8, '(1x,i3,4(f8.4),2(f10.4))') il, ajz(il, 1), spar(il, 1), vocc(il, 1), spnorm(il, 1), h2m(1)*&
                 & spkine(il, 1), spenrg(il, 1)
               Else
                  Write(8, '(61x,i3,4(f8.4),2(f10.4))') il, ajz(il, 2), spar(il, 2), vocc(il, 2), spnorm(il, 2), h2m(2)*&
                 & spkine(il, 2), spenrg(il, 2)
               End If
            End Do
!
            Write(8, '(2/,A,i2,A,/,2A)') ' c.m. etc. (principal axis =', 1, ' ):', '              <x**2>      <y**2>      <z**2> &
           &  ', '     <x>         <y>         <z>    '
            Write(8, '(A,1P,6E12.4)') '   total: ', x2mtot, y2mtot, z2mtot, xcmtot, ycmtot, zcmtot
            Write(8, '(A,1P,6E12.4)') ' neutron: ', x2mn, y2mn, z2mn, xcm(1), ycm(1), zcm(1)
            Write(8, '(A,1P,6E12.4)') '  proton: ', x2mp, y2mp, z2mp, xcm(2), ycm(2), zcm(2)
!
            Write(8, '(2/,A,i2,A,/,3A)') ' the moments of the density (principal axis =', 1, ' ):', '                part.num.   &
           & rms       q10x    ', '    q20         q22         q30         q40'
            Write(8, '(A,2F10.4,1P,5E12.4)') ' total (T=0): ', pnrtot, rmstot, q10totx, q20tot, q22tot, q30tot, q40tot
            Write(8, '(A,2F10.4,1P,5E12.4)') ' total (T=1): ', pnrtot1, rmstot1, q10tot1x, q20tot1, q22tot1, q30tot1, q40tot1
            Write(8, '(A,2F10.4,1P,5E12.4)') '     neutron: ', pnrn, rmsn, q10xn, q20(1), q22(1), q30n, q40n
            Write(8, '(A,2F10.4,1P,5E12.4)') '      proton: ', pnrp, rmsp, q10xp, q20(2), q22(2), q30p, q40p
         End If
      End If
!
      Deallocate(erear3, erearc, rhoplt, xwk, zwk)
      Call flush(8)
      Call flush(78)
      Call flush(79)
      Call flush(80)
      Return
End Subroutine tinfo
