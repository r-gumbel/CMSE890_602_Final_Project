Subroutine sinfo(ncolx, ncoly, ncolz, lpsi, wx, wy, wz, xclx, xcly, xclz, der1x, der1y, der1z, h2m, t3, x3, alpha, itrs, npsi, &
& npmin, vocc, spkine, spenrg, speflu, spnorm, delesum, spar, ajz, xcm, ycm, zcm, rho, psi, pswk1, pswk2, pswk3, ehft, ehf0, &
& ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ehf, ecorc, ecorr, ecorl, &
& erear3, erearc, erearl, rlam, iconstr, ehfy, alz, asz, epairn, epairp, gapn, gapp, fermin, fermip, efluct, isoedens, etheta)
!-----------------------------------------------------------------------
!     calculates static observables for printout
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: itrs
      Integer, Intent(In) :: npsi(2)
      Integer, Intent(In) :: npmin(2)
      Integer, Intent(In) :: iconstr
      Real(wp), Intent(In) :: t3
      Real(wp), Intent(In) :: x3
      Real(wp), Intent(In) :: alpha
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
      Real(wp), Intent(In) :: ecorc
      Real(wp), Intent(Out) :: ecorr
      Real(wp), Intent(Out) :: ecorl
      Real(wp), Intent(Out) :: ehf
      Real(wp), Intent(In) :: ehfc
      Real(wp), Intent(In) :: ehfy
      Real(wp), Intent(In) :: gapn
      Real(wp), Intent(In) :: gapp
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp), Intent(In) :: h2m(2)
      Real(wp), Intent(In) :: vocc(lpsi, 2)
      Real(wp), Intent(In) :: spkine(lpsi, 2)
      Real(wp), Intent(In) :: spenrg(lpsi, 2)
      Real(wp), Intent(In) :: speflu(lpsi, 2)
      Real(wp), Intent(In) :: spnorm(lpsi, 2)
      Real(wp), Intent(Out) :: erear3(lpsi, 2)
      Real(wp), Intent(Out) :: erearc(lpsi, 2)
      Real(wp), Intent(Out) :: erearl(lpsi, 2)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: rlam(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: epairn(2)
      Real(wp), Intent(In) :: epairp(2)
      Real(wp), Intent(In) :: fermin(2)
      Real(wp), Intent(In) :: fermip(2)
      Real(wp), Intent(Inout) :: delesum
      Real(wp), Intent(Inout) :: efluct
      Real(wp), Intent(Inout) :: spar(lpsi, 2)
      Real(wp), Intent(Inout) :: ajz(lpsi, 2)
      Real(wp), Intent(Inout) :: xcm(2)
      Real(wp), Intent(Inout) :: ycm(2)
      Real(wp), Intent(Inout) :: zcm(2)
      Real(wp), Intent(Inout) :: alz(lpsi, 2)
      Real(wp), Intent(Inout) :: asz(lpsi, 2)
      Real(wp), Intent(Inout) :: etheta(3)
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk2(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk3(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: e2 = 1.439978408596513_wp
      Integer, Parameter :: jz_exact = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, iz, iy, ix, nps, npsm, il, nst, is, ic, izz, iyy, ixx, npart, idp, nstp, lpsiq
      Real(wp), Allocatable, Dimension(:, :) :: q2sp, hjz, hlz, hsz, xja
      Real(wp) :: pnrq, xcmq, ycmq, zcmq, z, y, x, vol, pnrn, zcmn, ycmn, xcmn, pnrp, zcmp, ycmp, xcmp, pnrtot, xcmtot, ycmtot, &
     & zcmtot, rmsq, q10q, q30q, q40q, x2mom, y2mom, z2mom, xymom, xzmom, yzmom, zsq, ysq, xsq, rmsn, q10n, q30n, q40n, x2mn, &
     & y2mn, z2mn, xymomn, xzmomn, yzmomn, rmsp, q10p, q30p, q40p, x2mp, y2mp, z2mp, xymomp, xzmomp, yzmomp, rmstot, q10tot, &
     & q30tot, q40tot, x2mtot, y2mtot, z2mtot, xymomt, xzmomt, yzmomt, rhotot, rhop, rhon, tke, ecortot, ehfkoop, ehfint, &
     & ehfsum, tkerr, epair, epairhf, weight, psi2, cli, rsz, sigis, xpar, ehf_even, ehf_odd, q20(2), q22(2), q20tot, q22tot, &
     & dxq, dzq, rhoam1, slate, xjz(lpsi), isoscalarenergy, isovectorenergy, isoedens(ncolx,ncoly,ncolz,2)
      Real(wp) :: xlmass, xrmass, xlchrg, xrchrg, rhopp
      Real(wp) :: xi(3), zi(3,3)
      Complex(wp) :: psinq
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
                  vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                  pnrq = vol + pnrq
                  xcmq = vol * x + xcmq
                  ycmq = vol * y + ycmq
                  zcmq = vol * z + zcmq
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
!
      End Do
!
      pnrtot = pnrn + pnrp
      xcmtot = (xcmn+xcmp) / pnrtot
      ycmtot = (ycmn+ycmp) / pnrtot
      zcmtot = (zcmn+zcmp) / pnrtot
      xcm(1) = xcmn / pnrn
      xcm(2) = xcmp / pnrp
      ycm(1) = ycmn / pnrn
      ycm(2) = ycmp / pnrp
      zcm(1) = zcmn / pnrn
      zcm(2) = zcmp / pnrp
! calculate left/right mass and charge
      xlmass = 0.0_wp
      xrmass = 0.0_wp
      xlchrg = 0.0_wp
      xrchrg = 0.0_wp
!
      Do ix = 1, ncolx
         x = xclx(ix)
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               If(x < 0.0_wp) Then
                  vol = wx(ix)*wy(iy)*wz(iz)
                  rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                  rhopp = rho(ix, iy, iz, 2)
                  xlmass = xlmass + vol*rhotot
                  xlchrg = xlchrg + vol*rhopp
               Else
                  vol = wx(ix)*wy(iy)*wz(iz)
                  rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                  rhopp = rho(ix, iy, iz, 2)
                  xrmass = xrmass + vol*rhotot
                  xrchrg = xrchrg + vol*rhopp
               End If
            End Do
         End Do
      End Do
!      Write(*,*) xlmass, xlchrg, xrmass, xrchrg
      Call inertia(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcmtot, ycmtot, zcmtot, rho, xi, zi)
!      Write(8, '(/,A,/,A)') ' inertia/m', '         I(1)        I(2) I(3)'
!      Write(8, '(A,3F12.4)')  '         ', xi(1), xi(2), xi(3)
!-----------------------------------------------------------------------
!        calculate radii, moments, etc.
!
!        q10 = r*y10 = sqrt(3/4pi)*z  for isovector (rhop-rhon)
!        q20 = r**2*y20 = sqrt(5/16pi)*(2*z**2-x**2-y**2)
!        q22 = r**2*(y22+y2-2) = sqrt(15/8pi)*(x**2-y**2)
!        q30 = r**3*y30 = sqrt(7/4pi)*(z**3-1.5*z*(x**2+y**2))
!        q40 = r**4*y40 = sqrt(9/4pi)*(z**4-3*z**2*(x**2+y**2)
!                                         +3/8*(x**2+y**2)**2)
!
! first diagonalize the quadrupole tensor to obtain q2x and principal axis
! the moments are calculated based on principal axis (assuming x,y, or z).
!-----------------------------------------------------------------------
      Call moments(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, q20, q22, q20tot, q22tot, dxq, dzq, &
     & idp, etheta)
! for spherical systems set idp=3 to use z as axis of quantization
      If(Abs(q20tot) <= 1.e-1_wp) idp = 3
!
      Do iq = 1, 2
!
         rmsq = 0.0_wp
         q30q = 0.0_wp
         q40q = 0.0_wp
         x2mom = 0.0_wp
         y2mom = 0.0_wp
         z2mom = 0.0_wp
         xymom = 0.0_wp
         xzmom = 0.0_wp
         yzmom = 0.0_wp
         If(idp == 1) Then
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               zsq = z * z
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  ysq = y * y
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     xsq = x * x
                     vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                     rmsq = vol * (xsq + ysq + zsq) + rmsq
                     q30q = vol * (xsq*x - 1.5_wp*x*(zsq + ysq)) + q30q
                     q40q = vol * (xsq*xsq - 3.0_wp*xsq*(zsq + ysq) + 0.375_wp*(zsq + ysq)**2) + q40q
                     x2mom = vol * xsq + x2mom
                     y2mom = vol * ysq + y2mom
                     z2mom = vol * zsq + z2mom
                     xymom = vol * x * y + xymom
                     xzmom = vol * x * z + xzmom
                     yzmom = vol * y * z + yzmom
                  End Do
               End Do
            End Do
            q10q = xcm(iq)
         Else If(idp == 2) Then
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               zsq = z * z
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  ysq = y * y
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     xsq = x * x
                     vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                     rmsq = vol * (ysq + xsq + zsq) + rmsq
                     q10q = vol * y + q10q
                     q30q = vol * (ysq*y - 1.5_wp*y*(zsq + xsq)) + q30q
                     q40q = vol * (ysq*ysq - 3.0_wp*ysq*(zsq + xsq) + 0.375_wp*(zsq + xsq)**2) + q40q
                     x2mom = vol * xsq + x2mom
                     y2mom = vol * ysq + y2mom
                     z2mom = vol * zsq + z2mom
                     xymom = vol * x * y + xymom
                     xzmom = vol * x * z + xzmom
                     yzmom = vol * y * z + yzmom
                  End Do
               End Do
            End Do
            q10q = ycm(iq)
         Else If(idp == 3) Then
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               zsq = z * z
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  ysq = y * y
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     xsq = x * x
                     vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                     rmsq = vol * (xsq+ysq+zsq) + rmsq
                     q10q = vol * z + q10q
                     q30q = vol * (zsq*z - 1.5_wp*z*(xsq + ysq)) + q30q
                     q40q = vol * (zsq*zsq - 3.0_wp*zsq*(xsq + ysq) + 0.375_wp*(xsq + ysq)**2) + q40q
                     x2mom = vol * xsq + x2mom
                     y2mom = vol * ysq + y2mom
                     z2mom = vol * zsq + z2mom
                     xymom = vol * x * y + xymom
                     xzmom = vol * x * z + xzmom
                     yzmom = vol * y * z + yzmom
                  End Do
               End Do
            End Do
            q10q = zcm(iq)
         End If
!
         If(iq == 1) Then
            rmsn = rmsq
            q10n = q10q
            q30n = q30q
            q40n = q40q
            x2mn = x2mom
            y2mn = y2mom
            z2mn = z2mom
            xymomn = xymom
            xzmomn = xzmom
            yzmomn = yzmom
         Else
            rmsp = rmsq
            q10p = q10q
            q30p = q30q
            q40p = q40q
            x2mp = x2mom
            y2mp = y2mom
            z2mp = z2mom
            xymomp = xymom
            xzmomp = xzmom
            yzmomp = yzmom
         End If
      End Do
      q10n = Sqrt(3.0_wp/(4.0_wp*pi)) * q10n
      q10p = Sqrt(3.0_wp/(4.0_wp*pi)) * q10p
      q30n = Sqrt(7.0_wp/(4.0_wp*pi)) * q30n
      q30p = Sqrt(7.0_wp/(4.0_wp*pi)) * q30p
      q40n = Sqrt(9.0_wp/(4.0_wp*pi)) * q40n
      q40p = Sqrt(9.0_wp/(4.0_wp*pi)) * q40p
!
      rmstot = (rmsn+rmsp) / pnrtot
      rmstot = Sqrt(rmstot)
      rmsp = Sqrt(rmsp/pnrp)
      rmsn = Sqrt(rmsn/pnrn)
      q10tot = q10p - q10n
      q30tot = q30n + q30p
      q40tot = q40n + q40p
      x2mtot = (x2mn+x2mp) / pnrtot
      y2mtot = (y2mn+y2mp) / pnrtot
      z2mtot = (z2mn+z2mp) / pnrtot
      x2mn = x2mn / pnrn
      x2mp = x2mp / pnrp
      y2mn = y2mn / pnrn
      y2mp = y2mp / pnrp
      z2mn = z2mn / pnrn
      z2mp = z2mp / pnrp
      xymomn = xymomn / pnrn
      xzmomn = xzmomn / pnrn
      yzmomn = yzmomn / pnrn
      xymomp = xymomp / pnrp
      xzmomp = xzmomp / pnrp
      yzmomp = yzmomp / pnrp
      xymomt = xymomn + xymomp
      xzmomt = xzmomn + xzmomp
      yzmomt = yzmomn + yzmomp
!-----------------------------------------------------------------------
!     calculate the correction term to the hf energy due to the 3-body
!     force. this could be replaced but the answer is general for the
!     forms of skyrme force in use
!-----------------------------------------------------------------------
      ecorr = 0.0_wp
      ecorl = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               vol = wx(ix) * wy(iy) * wz(iz)
               rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               rhon = rho(ix, iy, iz, 1)
               rhop = rho(ix, iy, iz, 2)
               ecorr = ecorr + vol * rhotot ** alpha * ((1.0_wp+x3/2.0_wp)*rhotot**2-(0.50_wp+x3)*(rhop**2+rhon**2))
               ecorl = ecorl + vol * (rlam(ix, iy, iz, 1)*rhon+rlam(ix, iy, iz, 2)*rhop)
            End Do
         End Do
      End Do
      ecorr = - alpha * t3 / 24.0_wp * ecorr
!
      slate = - 3.0_wp / 4.0_wp * (3.0_wp/pi) ** (1.0_wp/3.0_wp) * e2
      erear3 = 0.0_wp
      erearc = 0.0_wp
      erearl = 0.0_wp
      ecortot = 0.0_wp
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            If(spenrg(nst, iq) > 0.0_wp) Cycle
            Do is = 1, 2
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        vol = wx(ix) * wy(iy) * wz(iz)
                        rhotot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                        rhon = rho(ix, iy, iz, 1)
                        rhop = rho(ix, iy, iz, 2)
                        rhoam1 = rhotot ** (alpha-1.0_wp+1.0e-25_wp)
                        psinq = psi(ix, iy, iz, is, nst, iq)
                        erear3(nst, iq) = erear3(nst, iq) - alpha * t3 / 24.0_wp * vol * vocc(nst, iq) * conjg(psinq) * psinq * &
                       & rhoam1 * ((1.0_wp+x3/2.0_wp)*rhotot**2-(0.50_wp+x3)*(rhop**2+rhon**2))
                        erearc(nst, iq) = erearc(nst, iq) + vol * slate / 3.0_wp * vocc(nst, iq) * conjg(psinq) * psinq * rhop ** &
                       & (1.0_wp/3.0_wp)
                        erearl(nst, iq) = erearl(nst, iq) - 0.5_wp * vol * rlam(ix, iy, iz, iq) * vocc(nst, iq) * conjg(psinq) * &
                       & psinq
                     End Do
                  End Do
               End Do
            End Do
            ecortot = ecortot + erear3(nst, iq)
            If(iq == 2) ecortot = ecortot + erearc(nst, iq)
            If(iconstr == 2) ecortot = ecortot + erearl(nst, iq)
         End Do
      End Do
      erearc(:, 1) = 0.0_wp ! slater correction only for protons
!---------------------------------------------------------------------------
!     calculate the hf energy and energy fluctuations.
!     instead of brute force integration of the energy density we
!     calculate the hf energy from the single particle energies, the
!     kinetic energies, and the correction term
!
!     note: there are two ways to calculate the hf energy: 1) using
!           the koopman + correction terms. 2) integrating the energy
!           density. The koopman expression has correction terms for
!           the coulomb exchange, three-body, and density constraint.
!----------------------------------------------------------------------------
      ehf = 0.0_wp
      tke = 0.0_wp
      efluct = 0.0_wp

      Do iq = 1, 2
         nps = npsi(iq)
         npsm = npmin(iq)
         Do il = npsm, nps
            If(spenrg(il, iq) > 0.0_wp) Cycle
            ehf = ehf + vocc(il, iq) * h2m(iq) * spkine(il, iq) + vocc(il, iq) * spenrg(il, iq)
            tke = tke + vocc(il, iq) * h2m(iq) * spkine(il, iq)
            efluct = efluct + vocc(il, iq) * speflu(il, iq)
         End Do
      End Do
      ehf = ehf / 2.0_wp
      ehfkoop = ehf + ecortot
!
      efluct = efluct / pnrtot
      delesum = delesum / pnrtot
! same as ehfkoop but correction terms calculated by integration
      ehf = ehf + ecorr + ecorc
      If(iconstr == 2) ehf = ehf - 0.5_wp*ecorl
      ehfint = ehft + ehf0 + ehf12 + ehf22 + ehf3 + ehfls + ehfc + ehfy + ehf0_odd + ehf12_odd + ehf22_odd + ehf3_odd + ehfls_odd &
     & + ehflj2 + ehfbj2
      ehfsum = tke + ehf0 + ehf12 + ehf22 + ehf3 + ehfls + ehfc + ehfy + ehf0_odd + ehf12_odd + ehf22_odd + ehf3_odd + ehfls_odd &
     & + ehflj2 + ehfbj2
      tkerr = ehft - tke
      ehf_odd = ehf0_odd + ehf12_odd + ehf22_odd + ehf3_odd + ehfls_odd + ehflj2
      ehf_even = ehfint - ehf_odd
!-----------------------------------------------------------------------
!         energy - total energy, ehf + epair(1)
!         t.k.e. - total kinetic energy
!-----------------------------------------------------------------------
      epair = epairn(1) + epairp(1) ! not used, leftover from L-N
      epairhf = epairn(2) + epairp(2)
      ehfsum = ehfsum - epairhf
      ehfint = ehfint - epairhf
      ehf = ehf - epairhf
!-----------------------------------------------------------------------
!     calculate the z-component of total angular momentum and parity
!     this calculation of the parity only makes sense for meshes that
!     are symmetric about the origin (not neccesarily equally spaced).
!-----------------------------------------------------------------------
      npart = Max(npsi(1), npsi(2))
      Allocate(q2sp(npart, 2))
!
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
!
            Do is = 1, 2
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        pswk3(ix, iy, iz, is) = psi(ix, iy, iz, is, nst, iq)
                     End Do
                  End Do
               End Do
            End Do
!
            q2sp(nst, iq) = 0.0
            Do is = 1, 2
               Do iz = 1, ncolz
                  z = xclz(iz) - zcm(iq)
                  zsq = z * z
                  Do iy = 1, ncoly
                     y = xcly(iy) - ycm(iq)
                     ysq = y * y
                     Do ix = 1, ncolx
                        x = xclx(ix) - xcm(iq)
                        xsq = x * x
                        If(idp == 1) Then
                           weight = wx(ix) * wy(iy) * wz(iz) * (xsq+xsq-zsq-ysq)
                        Else If(idp == 2) Then
                           weight = wx(ix) * wy(iy) * wz(iz) * (ysq+ysq-zsq-xsq)
                        Else If(idp == 3) Then
                           weight = wx(ix) * wy(iy) * wz(iz) * (zsq+zsq-xsq-ysq)
                        End If
                        psi2 = conjg(pswk3(ix, iy, iz, is)) * pswk3(ix, iy, iz, is)
                        q2sp(nst, iq) = q2sp(nst, iq) + weight * vocc(nst, iq) * psi2
                     End Do
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        calculate jz using idp (principal axis) as quantization axis
!-----------------------------------------------------------------------
            If(idp == 1) Then
               Call cmulz(ncolx, ncoly, ncolz, der1z, pswk3, pswk1, 0)
               Call cmuly(ncolx, ncoly, ncolz, der1y, pswk3, pswk2, 0)
            Else If(idp == 2) Then
               Call cmulx(ncolx, ncoly, ncolz, der1x, pswk3, pswk1, 0)
               Call cmulz(ncolx, ncoly, ncolz, der1z, pswk3, pswk2, 0)
            Else If(idp == 3) Then
               Call cmulx(ncolx, ncoly, ncolz, der1x, pswk3, pswk1, 0)
               Call cmuly(ncolx, ncoly, ncolz, der1y, pswk3, pswk2, 0)
            End If
!-----------------------------------------------------------------------
!        lz which requires combination of first derivatives
!-----------------------------------------------------------------------
            cli = 0.0_wp
            rsz = 0.0_wp
            Do is = 1, 2
               ic = 3 - is
               sigis = 3 - 2 * is
!
               Do iz = 1, ncolz
                  z = xclz(iz) - zcm(iq)
                  Do iy = 1, ncoly
                     y = xcly(iy) - ycm(iq)
                     Do ix = 1, ncolx
                        x = xclx(ix) - xcm(iq)
                        weight = wx(ix) * wy(iy) * wz(iz)
                        If(idp == 1) Then
                           cli = cli + weight * aimag((y*conjg(pswk3(ix, iy, iz, is))*pswk1(ix, iy, iz, is)-z*conjg(pswk3(ix, iy, &
                          & iz, is))*pswk2(ix, iy, iz, is)))
                           rsz = rsz + weight * conjg(pswk3(ix, iy, iz, is)) * pswk3(ix, iy, iz, ic)
                        Else If(idp == 2) Then
                           cli = cli + weight * aimag((z*conjg(pswk3(ix, iy, iz, is))*pswk1(ix, iy, iz, is)-x*conjg(pswk3(ix, iy, &
                          & iz, is))*pswk2(ix, iy, iz, is)))
                           rsz = rsz + 2.0_wp * sigis * weight * aimag(conjg(pswk3(ix, iy, iz, is))*pswk3(ix, iy, iz, ic))
                        Else If(idp == 3) Then
                           cli = cli + weight * aimag(x*conjg(pswk3(ix, iy, iz, is))*pswk2(ix, iy, iz, is)-y*conjg(pswk3(ix, iy, &
                          & iz, is))*pswk1(ix, iy, iz, is))
                           rsz = rsz + sigis * weight * conjg(pswk3(ix, iy, iz, is)) * pswk3(ix, iy, iz, is)
                        End If
                     End Do
                  End Do
               End Do
!
            End Do
!
            ajz(nst, iq) = cli + 0.5_wp * rsz
            alz(nst, iq) = cli
            asz(nst, iq) = 0.5_wp * rsz
!-----------------------------------------------------------------------
!        calculation of parity
!        NOTE: This only works for symmetric meshes, nucleus center at origin
!-----------------------------------------------------------------------
            xpar = 0.0_wp
            Do is = 1, 2
               Do iz = 1, ncolz
                  izz = ncolz - iz + 1
                  Do iy = 1, ncoly
                     iyy = ncoly - iy + 1
                     Do ix = 1, ncolx
                        ixx = ncolx - ix + 1
                        weight = wx(ix) * wy(iy) * wz(iz)
                        xpar = xpar + weight * conjg(pswk3(ix, iy, iz, is)) * pswk3(ixx, iyy, izz, is)
                     End Do
                  End Do
               End Do
            End Do
            spar(nst, iq) = xpar
!
         End Do
      End Do
! calculate jz exactly
      If(jz_exact == 1) Then
         Do iq = 1, 2
            lpsiq = npsi(iq) - npmin(iq) + 1
            Allocate(hjz(lpsiq, lpsiq), hlz(lpsiq, lpsiq), hsz(lpsiq, lpsiq), xja(lpsiq, lpsiq))
            Do nst = npmin(iq), npsi(iq)
               Do nstp = npmin(iq), npsi(iq)
                  pswk3(:, :, :, :) = psi(:, :, :, :, nstp, iq)
                  If(idp == 1) Then
                     Call cmulz(ncolx, ncoly, ncolz, der1z, pswk3, pswk1, 0)
                     Call cmuly(ncolx, ncoly, ncolz, der1y, pswk3, pswk2, 0)
                  Else If(idp == 2) Then
                     Call cmulx(ncolx, ncoly, ncolz, der1x, pswk3, pswk1, 0)
                     Call cmulz(ncolx, ncoly, ncolz, der1z, pswk3, pswk2, 0)
                  Else If(idp == 3) Then
                     Call cmulx(ncolx, ncoly, ncolz, der1x, pswk3, pswk1, 0)
                     Call cmuly(ncolx, ncoly, ncolz, der1y, pswk3, pswk2, 0)
                  End If
                  cli = 0.0_wp
                  rsz = 0.0_wp
                  Do is = 1, 2
                     sigis = 3 - 2 * is
                     ic = 3 - is
                     Do iz = 1, ncolz
                        z = xclz(iz) - zcm(iq)
                        Do iy = 1, ncoly
                           y = xcly(iy) - ycm(iq)
                           Do ix = 1, ncolx
                              x = xclx(ix) - xcm(iq)
                              weight = wx(ix) * wy(iy) * wz(iz)
                              If(idp == 1) Then
                                 cli = cli + weight * aimag((y*conjg(psi(ix, iy, iz, is, nst, iq))*pswk1(ix, iy, iz, &
                                & is)-z*conjg(psi(ix, iy, iz, is, nst, iq))*pswk2(ix, iy, iz, is)))
                                 rsz = rsz + weight * conjg(psi(ix, iy, iz, is, nst, iq)) * pswk3(ix, iy, iz, ic)
                              Else If(idp == 2) Then
                                 cli = cli + weight * aimag((z*conjg(psi(ix, iy, iz, is, nst, iq))*pswk1(ix, iy, iz, &
                                & is)-x*conjg(psi(ix, iy, iz, is, nst, iq))*pswk2(ix, iy, iz, is)))
                                 rsz = rsz + 2.0_wp * sigis * weight * aimag(conjg(psi(ix, iy, iz, is, nst, iq))*pswk3(ix, iy, &
                                & iz, ic))
                              Else If(idp == 3) Then
                                 cli = cli + weight * aimag(x*conjg(psi(ix, iy, iz, is, nst, iq))*pswk2(ix, iy, iz, &
                                & is)-y*conjg(psi(ix, iy, iz, is, nst, iq))*pswk1(ix, iy, iz, is))
                                 rsz = rsz + sigis * weight * conjg(psi(ix, iy, iz, is, nst, iq)) * pswk3(ix, iy, iz, is)
                              End If
                           End Do
                        End Do
                     End Do
                  End Do
                  hjz(nst-npmin(iq)+1, nstp-npmin(iq)+1) = cli + 0.5_wp * rsz
                  hlz(nst-npmin(iq)+1, nstp-npmin(iq)+1) = cli
                  hsz(nst-npmin(iq)+1, nstp-npmin(iq)+1) = 0.5_wp * rsz
               End Do
            End Do
            Call jacobi(hjz, lpsiq, lpsiq, xjz, xja)
            Do is = npmin(iq), npsi(iq)
               ajz(is, iq) = xjz(is-npmin(iq)+1)
            End Do
            Call jacobi(hlz, lpsiq, lpsiq, xjz, xja)
            Do is = npmin(iq), npsi(iq)
               alz(is, iq) = xjz(is-npmin(iq)+1)
            End Do
            Call jacobi(hsz, lpsiq, lpsiq, xjz, xja)
            Do is = npmin(iq), npsi(iq)
               asz(is, iq) = xjz(is-npmin(iq)+1)
            End Do
            Deallocate(hjz, hlz, hsz, xja)
         End Do
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

      Write(8, '(3/,A,I10,/,A,F12.4,/,A,F12.4,/,A,1P,E12.4,/,A,E12.4,/,A,E12.4)') ' iteration         = ', itrs, ' energy      (m&
     &ev) = ', ehfint, ' t.k.e.      (mev) = ', ehft, ' de/e              = ', delesum, ' ener. fluc. (mev) = ', efluct, ' tke er&
     &ror   (mev) = ', tkerr
!
      Write(8, '(/,A)') '                             total             time-even           time-odd'
      Write(8, '(A,1P,E15.6,4X,E15.6,4X,E15.6,/,A,E15.6,4X,E15.6,4X,E15.6,/,A,   E15.6,4X,E15.6,4X,E15.6,/,A,E15.6,4X,E15.6,4X,E1&
     &5.6,/,A,   E15.6,4X,E15.6,4X,E15.6,/,A,   E15.6,23X,E15.6,/,A,   E15.6,4X,E15.6,/,A,E15.6,/,A,E15.6,/,A,E15.6,2/     5(A,E1&
     &5.6,1X,E15.6/),2/)') ' intg. energy    (mev)= ', ehfint, ehf_even, ehf_odd, ' t0 contribution (mev)= ', ehf0 + ehf0_odd, &
     & ehf0, ehf0_odd, ' t1,t2 terms     (mev)= ', ehf12 + ehf22 + ehf12_odd + ehf22_odd, ehf12 + ehf22, ehf12_odd + ehf22_odd, '&
     & t3 contribution (mev)= ', ehf3 + ehf3_odd, ehf3, ehf3_odd, ' t4 contribution (mev)= ', ehfls + ehfls_odd, ehfls, &
     & ehfls_odd, ' j**2 term       (mev)= ', ehflj2, ehflj2, ' J**2 term       (mev)= ', ehfbj2, ehfbj2, ' yukawa energy   (mev)&
     &= ', ehfy, ' coulomb energy  (mev)= ', ehfc, ' Koopman energy  (mev)= ', ehfkoop, &
     &  ' pairing correction (n,p) (mev) = ', epairn(2), epairp(2), &
     &  ' fermi energy (n,p)       (mev) = ', fermin(1), fermip(1), &
     &  ' av. gap (n,p)            (mev) = ', fermin(2), fermip(2), &
     &  ' av. delta (n,p)          (mev) = ', gapn, gapp
     If(iconstr == 2) Write(8,'(A,1P,E15.6)') ' Constraint energy (mev)= ',ecorl
     Write(8,'(A,1P,E15.6,/,A,1P,E15.6,/,A,1P,E15.6)') ' Isovector energy  (mev)= ', isovectorenergy, &
     & ' Isoscalar energy  (mev)= ', isoscalarenergy, &
     & ' IsoTotal energy   (mev)= ', isovectorenergy+isoscalarenergy
!-----------------------------------------------------------------------
!     the following block Writes the single-particle quantities
!-----------------------------------------------------------------------
      Write(8, '(/,a)') ' neutron single particle states:'
      Write(8, '(a)') '  nr    jz      lz      sz     par    weight    efluct     norm      q20         ekin    energy'
      Do il = npmin(1), npsi(1)
         Write(8, '(1x,i3,5(f8.4),1pe12.4,0p,f8.4,1p,e12.4,0p,2(f10.4))') il, ajz(il, 1), alz(il, 1), asz(il, 1), spar(il, 1), &
        & vocc(il, 1), speflu(il, 1), spnorm(il, 1), q2sp(il, 1), h2m(1) * spkine(il, 1), spenrg(il, 1)
      End Do
      Write(8, '(/,a)') ' proton single particle states:'
      Write(8, '(a)') '  nr    jz      lz      sz     par    weight    efluct     norm      q20         ekin    energy'
      Do il = npmin(2), npsi(2)
         Write(8, '(1x,i3,5(f8.4),1pe12.4,0p,f8.4,1p,e12.4,0p,2(f10.4))') il, ajz(il, 2), alz(il, 2), asz(il, 2), spar(il, 2), &
        & vocc(il, 2), speflu(il, 2), spnorm(il, 2), q2sp(il, 2), h2m(1) * spkine(il, 2), spenrg(il, 2)
      End Do
!
      Write(8, '(2/,A,i2,A,/,3A)') ' the moments of the density (principal axis =', idp, ' ):', '               part.num.   rms-r&
     &adius', '   <x**2>      <y**2>      <z**2>   ', '  <x>         <y>         <z>    '
      Write(8, '(A,2F12.4,1P,6E12.4)') '    total: ', pnrtot, rmstot, x2mtot, y2mtot, z2mtot, xcmtot, ycmtot, zcmtot
      Write(8, '(A,2F12.4,1P,6E12.4)') '  neutron: ', pnrn, rmsn, x2mn, y2mn, z2mn, xcm(1), ycm(1), zcm(1)
      Write(8, '(A,2F12.4,1P,6E12.4)') '   proton: ', pnrp, rmsp, x2mp, y2mp, z2mp, xcm(2), ycm(2), zcm(2)
!
      Write(8, '(2/,A,i2,A,/,3A)') ' the moments of the density continued (principal axis =', idp, ' ):', '              q10(t=1)&
     &    q20(t=0)', '    q22(t=0)    q30(t=0)    q40(t=0)', '       x*y         x*z         y*z'
      Write(8, '(A,1P,8E12.4)') '    total: ', q10tot, q20tot, q22tot, q30tot, q40tot, xymomt, xzmomt, yzmomt
      Write(8, '(A,12X,1P,7E12.4)') '  neutron: ', q20(1), q22(1), q30n, q40n, xymomn, xzmomn, yzmomn
      Write(8, '(A,12X,1P,7E12.4)') '   proton: ', q20(2), q22(2), q30p, q40p, xymomp, xzmomp, yzmomp
!
      Deallocate(q2sp)
!
      Call flush(8)
      Return
End Subroutine sinfo
