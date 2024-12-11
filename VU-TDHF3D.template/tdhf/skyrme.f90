Subroutine skyrme(ncolx, ncoly, ncolz, nx2, ny2, nz2, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, der1x, &
& der1y, der1z, der2x, der2y, der2z, wx, wy, wz, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, bmunu, dbmu, upot, &
& bmass, xiq, hsigma, cq, dcq, tmpvec, wcoul, wyuk, dbmass, dxiq, icoul, itheta, itimrev, iodds, wxyz, coulplan1, coulplan2, &
& rho2, q, iperc)
!
!-----------------------------------------------------------------------
!        calculates the Skyrme mean-field potential U_q as well as all
!        quantities needed in calculating h_HF on psi. Time-odd parts
!        are only calculated if there is no time-reversal invariance
!        (itimrev=0).
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: nx2
      Integer, Intent(In) :: ny2
      Integer, Intent(In) :: nz2
      Integer, Intent(In) :: icoul, itimrev, itheta, iodds
      Integer, Intent(In) :: iperc
      Integer(8), Intent(In) :: coulplan1, coulplan2
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
      Real(wp), Intent(In) :: h2m(2)
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: tau(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: sodens(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: spinden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: kinvden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: spincur(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp) :: bmass(ncolx, ncoly, ncolz, 2)
      Real(wp) :: xiq(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: hsigma(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: cq(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: dcq(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: tmpvec(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp) :: wyuk(ncolx, ncoly, ncolz)
      Real(wp) :: dbmass(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: dxiq(ncolx, ncoly, ncolz, 2)
      Real(wp) :: delxj(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: bmunu(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: dbmu(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: upot(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: q(nx2, ny2, nz2), rho2(nx2, ny2, nz2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Integer, Parameter :: ioddls = 1 ! should be 1 except special circumstances
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: e2 = 1.439978408596513_wp
      Real(wp), Parameter :: small = 1.0e-25_wp
      Real(wp), Parameter :: zero = 0.0_wp
      Real(wp), Parameter :: half = 0.5_wp
      Real(wp), Parameter :: one = 1.0_wp
      Real(wp), Parameter :: two = 2.0_wp
      Real(wp), Parameter :: three = 3.0_wp
      Real(wp), Parameter :: four = 4.0_wp
      Real(wp), Parameter :: eight = 8.0_wp
      Real(wp), Parameter :: twelve = 12.0_wp
      Real(wp), Parameter :: sixteen = 16.0_wp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, iz, iy, ix, ic, icomp, imu, inu
      Real(wp) :: t3a, t3b, rhoq, rhon, rhop, rhot, djq, djtot, slate, t0a, t0b, rhotot, t1a, t1b, tauq, tautot, t1c, t1d, &
     & d2rhoq, d2rhot, t4a, t4b, dxrhot, dyrhot, dzrhot, scq, sctot, s2t, s2n, s2p, st, sq, d2st, d2sq, tkt, tkq, dxj1, dxj2
!-----------------------------------------------------------------------
!        three-body contribution to mean-field U_q (t3 terms)
!-----------------------------------------------------------------------
      t3a = t3 / twelve * (one+x3/two)
      t3b = t3 / twelve * (one/two+x3)
      Do iq = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  rhoq = rho(ix, iy, iz, iq)
                  rhon = rho(ix, iy, iz, 1)
                  rhop = rho(ix, iy, iz, 2)
                  rhot = rhon + rhop
                  upot(ix, iy, iz, iq) = t3a * (alpha+two) * rhot ** (alpha+one) - t3b * alpha * rhot ** (alpha-one+small) * &
                 & (rhon*rhon+rhop*rhop) - t3b * two * rhot ** alpha * rhoq
               End Do
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        three-body contribution to mean-field spin-density (t3 terms) (time-odd)
!-----------------------------------------------------------------------
      If(itimrev == 0 .And. iodds == 1) Then
!  square the spin-density and use bmass as temporary storage
         Do iq = 1, 2
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     bmass(ix, iy, iz, iq) = spinden(ix, iy, iz, 1, iq) ** 2 + spinden(ix, iy, iz, 2, iq) ** 2 + spinden(ix, iy, &
                    & iz, 3, iq) ** 2
                  End Do
               End Do
            End Do
         End Do
!  now add to the mean-field
         t3a = t3 / (twelve*two) * alpha
         Do iq = 1, 2
            ic = 3 - iq
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                     s2t = (spinden(ix, iy, iz, 1, 1)+spinden(ix, iy, iz, 1, 2)) ** 2 + (spinden(ix, iy, iz, 2, 1)+spinden(ix, &
                    & iy, iz, 2, 2)) ** 2 + (spinden(ix, iy, iz, 3, 1)+spinden(ix, iy, iz, 3, 2)) ** 2
!                     s2t  = bmass(ix,iy,iz,1) + bmass(ix,iy,iz,2)&
!                          + 2.0_wp*(spinden(ix,iy,iz,1,iq)*spinden(ix,iy,iz,1,ic)&
!                          +         spinden(ix,iy,iz,2,iq)*spinden(ix,iy,iz,2,ic)&
!                          +         spinden(ix,iy,iz,3,iq)*spinden(ix,iy,iz,3,ic))
                     s2n = bmass(ix, iy, iz, 1)
                     s2p = bmass(ix, iy, iz, 2)
                     upot(ix, iy, iz, iq) = upot(ix, iy, iz, iq) + t3a * x3 * rhot ** (alpha-one+small) * s2t - t3a * rhot ** &
                    & (alpha-one+small) * (s2n+s2p)
                  End Do
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        contribution of the spin-orbit potential to mean-field U_q (t4 terms)
!        first calculate the divergence of the spin-orbit current.
!        Use bmass as temporary storage
!-----------------------------------------------------------------------
      Do iq = 1, 2
         Call divergence(ncolx, ncoly, ncolz, der1x, der1y, der1z, sodens(1, 1, 1, 1, iq), bmass(1, 1, 1, iq))
      End Do
!
      Do iq = 1, 2
         ic = 3 - iq
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  djq = bmass(ix, iy, iz, iq)
                  djtot = djq + bmass(ix, iy, iz, ic)
                  upot(ix, iy, iz, iq) = upot(ix, iy, iz, iq) - t4 / two * djtot - t4p / two * djq
               End Do
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        solve poisson equation for the coulomb potential (wcoul).
!        exchange part is the slater approximation
!        all other arrays used as workspace.
!-----------------------------------------------------------------------
      If(icoul /= 0) Then
         Call pois(iperc, coulplan1, coulplan2, ncolx, ncoly, ncolz, nx2, ny2, nz2, wx(1), wy(1), wz(1),&
                   rho(1, 1, 1, 2), rho2, q, wcoul)
!
         slate = (three/pi) ** (one/three) * e2
         upot(:, :, :, 2) = upot(:, :, :, 2) + wcoul - slate * rho(:, :, :, 2) ** (one/three)
      End If
!-------------------------------------------------------------------------
!        solve helmholtz equation for the Yukawa potential (wyuk) for BKN
!        all other arrays used as workspace.
!-------------------------------------------------------------------------
      If(ayuk > 1.0e-5_wp) Then
         Call yukrun(rho, wyuk, dxiq(1, 1, 1, 1), dxiq(1, 1, 1, 2), bmass(1, 1, 1, 1), bmass(1, 1, 1, 2), der2x, der2y, der2z, &
        & wxyz, ayuk, ncolx, ncoly, ncolz)
         upot(:, :, :, 1) = upot(:, :, :, 1) + vyuk * wyuk
         upot(:, :, :, 2) = upot(:, :, :, 2) + vyuk * wyuk
      End If
!-----------------------------------------------------------------------
!        store the gradient of rho in tmpvec
!        these will be used in constructing the B_q term
!        also calculate the laplacian of rho (store in dxiq)
!-----------------------------------------------------------------------
      Do iq = 1, 2
         Call gradient(ncolx, ncoly, ncolz, der1x, der1y, der1z, rho(1, 1, 1, iq), tmpvec(1, 1, 1, 1, iq))
!
         Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, rho(1, 1, 1, iq), dxiq(1, 1, 1, iq))
      End Do
!-----------------------------------------------------------------------
!        construct the anti-symmetric drhobar (in delxj) tensor (used in B_q)
!-----------------------------------------------------------------------
      t4a = t4 / two
      t4b = t4p / two
      delxj = zero
      Do iq = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  dxrhot = tmpvec(ix, iy, iz, 1, 1) + tmpvec(ix, iy, iz, 1, 2)
                  dyrhot = tmpvec(ix, iy, iz, 2, 1) + tmpvec(ix, iy, iz, 2, 2)
                  dzrhot = tmpvec(ix, iy, iz, 3, 1) + tmpvec(ix, iy, iz, 3, 2)
                  delxj(ix, iy, iz, 1, 2, iq) = t4a * dzrhot + t4b * tmpvec(ix, iy, iz, 3, iq)
                  delxj(ix, iy, iz, 1, 3, iq) = - (t4a*dyrhot+t4b*tmpvec(ix, iy, iz, 2, iq))
                  delxj(ix, iy, iz, 2, 1, iq) = - delxj(ix, iy, iz, 1, 2, iq)
                  delxj(ix, iy, iz, 2, 3, iq) = t4a * dxrhot + t4b * tmpvec(ix, iy, iz, 1, iq)
                  delxj(ix, iy, iz, 3, 1, iq) = - delxj(ix, iy, iz, 1, 3, iq)
                  delxj(ix, iy, iz, 3, 2, iq) = - delxj(ix, iy, iz, 2, 3, iq)
               End Do
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        construct the B_mu_nu tensor (not neccessarily odd, e.g. SLy5)
!-----------------------------------------------------------------------
      bmunu = 0.0_wp
      If(itheta == 1) Then
         t1a = - (t2-t1) / four
         t1b = - (t1*x1+t2*x2) / four
         Do iq = 1, 2
            ic = 3 - iq
            Do imu = 1, 3
               Do inu = 1, 3
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           scq = spincur(ix, iy, iz, imu, inu, iq)
                           sctot = scq + spincur(ix, iy, iz, imu, inu, ic)
                           bmunu(ix, iy, iz, imu, inu, iq) = t1a * scq + t1b * sctot
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If
      Do iq = 1, 2
         ic = 3 - iq
         Do imu = 1, 3
            Do inu = 1, 3
               Do iz = 1, ncolz
                  Do iy = 1, ncoly
                     Do ix = 1, ncolx
                        bmunu(ix, iy, iz, imu, inu, iq) = bmunu(ix, iy, iz, imu, inu, iq) + delxj(ix, iy, iz, imu, inu, iq)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        construct the dBmu vector (db_mu = d_nu*b_nu_mu)
!-----------------------------------------------------------------------
      Do iq = 1, 2
         Do imu = 1, 3
            Call divergence(ncolx, ncoly, ncolz, der1x, der1y, der1z, bmunu(1, 1, 1, 1, imu, iq), dbmu(1, 1, 1, imu, iq))
         End Do
      End Do
!-----------------------------------------------------------------------
!        construct the C_q term depending on spin-density
!-----------------------------------------------------------------------
      cq = zero
      dcq = zero
      If(itheta == 1 .And. iodds == 1) Then
         t1a = (t2-t1) / eight
         t1b = (t1*x1+t2*x2) / eight
         Do iq = 1, 2
            ic = 3 - iq
            Do imu = 1, 3
               Do iz = 1, ncolz
                  Do iy = 1, ncoly
                     Do ix = 1, ncolx
                        scq = spinden(ix, iy, iz, imu, iq)
                        sctot = scq + spinden(ix, iy, iz, imu, ic)
                        cq(ix, iy, iz, imu, iq) = t1a * scq + t1b * sctot
                     End Do
                  End Do
               End Do
            End Do
         End Do
!-----------------------------------------------------------------------
!        construct the gradient tensor of dcq_mu_nu = d_mu C_q_nu term
!-----------------------------------------------------------------------
         Do iq = 1, 2
            Do inu = 1, 3
               Call gradient(ncolx, ncoly, ncolz, der1x, der1y, der1z, cq(1, 1, 1, inu, iq), dcq(1, 1, 1, 1, inu, iq))
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        add all two-body terms
!-----------------------------------------------------------------------
      xiq = zero
      bmass = zero
! calculate the gradient of the spin-current as a tensor to use in cross product in I_q term
! this part gives zero since in del.I_q we have del.delXs=0
      If(itimrev == 0 .And. ioddls == 1) Then
         Do iq = 1, 2
            Do imu = 1, 3
               Call gradient(ncolx, ncoly, ncolz, der1x, der1y, der1z, spinden(1, 1, 1, imu, iq), delxj(1, 1, 1, 1, imu, iq))
            End Do
         End Do
      End If
!
      Do iq = 1, 2
         ic = 3 - iq
!-----------------------------------------------------------------------
!        t0-dependent part
!-----------------------------------------------------------------------
         t0a = t0 * (one+x0/two)
         t0b = t0 * (half+x0)
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  rhoq = rho(ix, iy, iz, iq)
                  rhotot = rhoq + rho(ix, iy, iz, ic)
                  upot(ix, iy, iz, iq) = upot(ix, iy, iz, iq) + t0a * rhotot - t0b * rhoq
               End Do
            End Do
         End Do
!-----------------------------------------------------------------------
!        t1,t2 and tau-dependent part
!-----------------------------------------------------------------------
         t1a = (t1*(one+x1/two)+t2*(one+x2/two)) / four
         t1b = (t2*(x2+half)-t1*(x1+half)) / four
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  tauq = tau(ix, iy, iz, iq)
                  tautot = tauq + tau(ix, iy, iz, ic)
                  upot(ix, iy, iz, iq) = upot(ix, iy, iz, iq) + t1a * tautot + t1b * tauq
               End Do
            End Do
         End Do
!-----------------------------------------------------------------------
!        two-body laplacian*rho-dependent part
!        laplacian of rho was stored in dxiq
!-----------------------------------------------------------------------
         t1c = (three*t1*(one+x1/two)-t2*(one+x2/two)) / eight
         t1d = (three*t1*(x1+half)+t2*(x2+half)) / eight
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  d2rhoq = dxiq(ix, iy, iz, iq)
                  d2rhot = d2rhoq + dxiq(ix, iy, iz, ic)
                  upot(ix, iy, iz, iq) = upot(ix, iy, iz, iq) - t1c * d2rhot + t1d * d2rhoq
               End Do
            End Do
         End Do
!-----------------------------------------------------------------------
!        store the effective mass in bmass
!-----------------------------------------------------------------------
         t1a = (t1*(one+x1/two)+t2*(one+x2/two)) / four
         t1b = (t2*(half+x2)-t1*(half+x1)) / four
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  rhoq = rho(ix, iy, iz, iq)
                  rhotot = rhoq + rho(ix, iy, iz, ic)
                  bmass(ix, iy, iz, iq) = bmass(ix, iy, iz, iq) + h2m(iq) + t1a * rhotot + t1b * rhoq
               End Do
            End Do
         End Do
!-----------------------------------------------------------------------
!        begin constructing the iq-vector
!        these will be used in constructing the (del.iq+iq.del) term
!        they are stored in xiq
!        These are zero for the time-reversal invariant case
!-----------------------------------------------------------------------
         If(itimrev == 0) Then
            t1a = (t1*(one+x1/two)+t2*(one+x2/two)) / two
            t1b = (t2*(half+x2)-t1*(half+x1)) / two
            Do icomp = 1, 3
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        djq = currnt(ix, iy, iz, icomp, iq)
                        djtot = djq + currnt(ix, iy, iz, icomp, ic)
                        xiq(ix, iy, iz, icomp, iq) = - t1a * djtot - t1b * djq
                     End Do
                  End Do
               End Do
            End Do
!  add the spin-density cross-product term Del x (s+s_q) SHOULD GIVE ZERO in Del.Iq part!
            If(ioddls == 1) Then
               t4a = - t4 / two
               t4b = - t4p / two
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        dxj1 = (t4a+t4b) * delxj(ix, iy, iz, 2, 3, iq) + t4a * delxj(ix, iy, iz, 2, 3, ic)
                        dxj2 = (t4a+t4b) * delxj(ix, iy, iz, 3, 2, iq) + t4a * delxj(ix, iy, iz, 3, 2, ic)
                        xiq(ix, iy, iz, 1, iq) = xiq(ix, iy, iz, 1, iq) + (dxj1-dxj2)
                        dxj1 = (t4a+t4b) * delxj(ix, iy, iz, 1, 3, iq) + t4a * delxj(ix, iy, iz, 1, 3, ic)
                        dxj2 = (t4a+t4b) * delxj(ix, iy, iz, 3, 1, iq) + t4a * delxj(ix, iy, iz, 3, 1, ic)
                        xiq(ix, iy, iz, 2, iq) = xiq(ix, iy, iz, 2, iq) - (dxj1-dxj2)
                        dxj1 = (t4a+t4b) * delxj(ix, iy, iz, 1, 2, iq) + t4a * delxj(ix, iy, iz, 1, 2, ic)
                        dxj2 = (t4a+t4b) * delxj(ix, iy, iz, 2, 1, iq) + t4a * delxj(ix, iy, iz, 2, 1, ic)
                        xiq(ix, iy, iz, 3, iq) = xiq(ix, iy, iz, 3, iq) + (dxj1-dxj2)
                     End Do
                  End Do
               End Do
            End If
         End If
      End Do
!-----------------------------------------------------------------------
!        calculate the divergence of xiq and store in dxiq
!-----------------------------------------------------------------------
      If(itimrev == 0) Then
         Do iq = 1, 2
            Call divergence(ncolx, ncoly, ncolz, der1x, der1y, der1z, xiq(1, 1, 1, 1, iq), dxiq(1, 1, 1, iq))
         End Do
      End If
!-----------------------------------------------------------------------
!        calculate the gradient of the effective mass and store in
!        dbmass. this is used as part of the kinetic energy term in
!        hpsi
!-----------------------------------------------------------------------
      Do iq = 1, 2
         Call gradient(ncolx, ncoly, ncolz, der1x, der1y, der1z, bmass(1, 1, 1, iq), dbmass(1, 1, 1, 1, iq))
      End Do
!-----------------------------------------------------------------------
!        start building the Sigma_q term (time-odd)
!-----------------------------------------------------------------------
      hsigma = zero
      If(itimrev == 0) Then
!  calculate the laplacian of the spin-density vector and store in tmpvec
         Do iq = 1, 2
            Do icomp = 1, 3
               Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, spinden(1, 1, 1, icomp, iq), tmpvec(1, 1, 1, icomp, iq))
            End Do
         End Do
!  calculate the gradient of the current as a tensor to use in cross product
!  and store in delxj to use in calculating the total energy
         Do iq = 1, 2
            Do imu = 1, 3
               Call gradient(ncolx, ncoly, ncolz, der1x, der1y, der1z, currnt(1, 1, 1, imu, iq), delxj(1, 1, 1, 1, imu, iq))
            End Do
         End Do
!  add the spin-orbit cross-product term Del x (j+j_q)
         If(ioddls == 1) Then
            t4a = - t4 / two
            t4b = - t4p / two
            Do iq = 1, 2
               ic = 3 - iq
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        dxj1 = (t4a+t4b) * delxj(ix, iy, iz, 2, 3, iq) + t4a * delxj(ix, iy, iz, 2, 3, ic)
                        dxj2 = (t4a+t4b) * delxj(ix, iy, iz, 3, 2, iq) + t4a * delxj(ix, iy, iz, 3, 2, ic)
                        hsigma(ix, iy, iz, 1, iq) = hsigma(ix, iy, iz, 1, iq) + (dxj1-dxj2)
                        dxj1 = (t4a+t4b) * delxj(ix, iy, iz, 1, 3, iq) + t4a * delxj(ix, iy, iz, 1, 3, ic)
                        dxj2 = (t4a+t4b) * delxj(ix, iy, iz, 3, 1, iq) + t4a * delxj(ix, iy, iz, 3, 1, ic)
                        hsigma(ix, iy, iz, 2, iq) = hsigma(ix, iy, iz, 2, iq) - (dxj1-dxj2)
                        dxj1 = (t4a+t4b) * delxj(ix, iy, iz, 1, 2, iq) + t4a * delxj(ix, iy, iz, 1, 2, ic)
                        dxj2 = (t4a+t4b) * delxj(ix, iy, iz, 2, 1, iq) + t4a * delxj(ix, iy, iz, 2, 1, ic)
                        hsigma(ix, iy, iz, 3, iq) = hsigma(ix, iy, iz, 3, iq) + (dxj1-dxj2)
                     End Do
                  End Do
               End Do
            End Do
         End If
         If(iodds == 1) Then
            t1a = (t2+three*t1) / sixteen
            t1b = (t2*x2-three*t1*x1) / sixteen
            If(itheta == 1) Then
               t1c = (t2-t1) / eight
               t1d = (t1*x1+t2*x2) / eight
            Else
               t1c = zero
               t1d = zero
            End If
            Do iq = 1, 2
               Do icomp = 1, 3
                  Do ix = 1, ncolx
                     Do iy = 1, ncoly
                        Do iz = 1, ncolz
                           st = spinden(ix, iy, iz, icomp, 1) + spinden(ix, iy, iz, icomp, 2)
                           sq = spinden(ix, iy, iz, icomp, iq)
                           d2st = tmpvec(ix, iy, iz, icomp, 1) + tmpvec(ix, iy, iz, icomp, 2)
                           d2sq = tmpvec(ix, iy, iz, icomp, iq)
                           rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                           tkt = kinvden(ix, iy, iz, icomp, 1) + kinvden(ix, iy, iz, icomp, 2)
                           tkq = kinvden(ix, iy, iz, icomp, iq)
                           hsigma(ix, iy, iz, icomp, iq) = hsigma(ix, iy, iz, icomp, iq) + t0 / two * (x0*st-sq) + t3 / twelve * &
                          & rhot ** (alpha) * (x3*st-sq) + t1c * tkq + t1d * tkt
!                                                     + t1a*d2sq + t1b*d2st - OMITTED DUE TO INSTABILITY
!
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End If
      End If
!
      Return
End Subroutine skyrme
