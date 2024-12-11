Subroutine isoenergy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, der1y, &
& der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, wyuk, isoedens, worka, &
& icoul, itheta, itimrev, iodds)
!
!-----------------------------------------------------------------------
!        calculates the hf energy by integrating the energy functional
!        parameters t0...x3+coulomb+spin-orbit
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: icoul
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: itheta
      Integer, Intent(In) :: iodds
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
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: tau(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: sodens(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: spinden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: kinvden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(In) :: spincur(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp), Intent(In) :: delxj(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp), Intent(In) :: wcoul(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: wyuk(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp), Intent(Out) :: isoedens(ncolx, ncoly, ncolz, 2)
      Real(wp) :: cdens(ncolx, ncoly, ncolz)
      Real(wp) :: worka(ncolx, ncoly, ncolz, 2)
      Real(wp) :: workb(ncolx, ncoly, ncolz, 2)
      Real(wp) :: workc(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: tmpvec(ncolx, ncoly, ncolz, 3, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Integer, Parameter :: ioddls = 1 ! should be 1 except special circumstances
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: e2 = 1.439978408596513_wp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iz, iy, ix, iq, imu, inu
      Real(wp) :: rho2vec, delrhovec, ctauvec, jmunu2vec, deljvec
      Real(wp) :: weight, rhot, rhop, d2rho, taut, taun, taup, xj2, &
     &            dj, slate, wcr, s2t, stt, bjt, isorho, ehfc, ecorc, ehft
      Real(wp) :: c0r0, c0r1, c0s0, c0s1, c3r0, c3r1, c3s0, c3s1, cdr0, cdr1, ctau0, ctau1, cds0, cds1, ct0, ct1, cdj0, cdj1
!      real(wp) ::  sd2st, sd2sn, sd2sp
!
      rho2vec = 0.0_wp
      delrhovec = 0.0_wp
      ctauvec = 0.0_wp
      jmunu2vec = 0.0_wp
      deljvec = 0.0_wp
!----------------------------------------------------------------------
! Setting the coupling constants (t0 & t3)
!----------------------------------------------------------------------
      c0r0 = 3.0_wp/8.0_wp * t0
      c3r0 = 1.0_wp/16.0_wp * t3
      c0r1 = -1.0_wp/8.0_wp * t0 * (1.0_wp + 2.0_wp * x0)
      c3r1 = -1.0_wp/48.0_wp * t3 * (1.0_wp + 2.0_wp * x3)
      c0s0 = 1.0_wp/8.0_wp * t0 * (2.0_wp * x0 - 1.0_wp)
      c3s0 = 1.0_wp/48.0_wp * t3 * (2.0_wp * x3 - 1.0_wp)
      c0s1 = -1.0_wp/8.0_wp * t0
      c3s1 = -1.0_wp/48.0_wp * t3
!----------------------------------------------------------------------
! Setting the coupling constants (t1, t2, etc.)
!----------------------------------------------------------------------
      cdr0 = 1.0_wp/64.0_wp * (5.0_wp * t2 + 4.0_wp * t2 * x2 - 9.0_wp * t1)
      cdr1 = 1.0_wp/64.0_wp * (3.0_wp * t1 * (1.0_wp + 2.0_wp * x1) + t2 * (1.0_wp + 2.0_wp * x2))
      ctau0 = 1.0_wp/16.0_wp * (3.0_wp * t1 + t2 * (5.0_wp + 4.0_wp * x2))
      ctau1 = 1.0_wp/16.0_wp * (t2 * (1.0_wp + 2.0_wp * x2) - t1 * (1.0_wp + 2.0_wp * x1))
      cds0 = 1.0_wp/64.0_wp * (3.0_wp * t1 * (1.0_wp - 2.0_wp * x1) + t2 * (1.0_wp + 2.0_wp * x2))
      cds1 = 1.0_wp/64.0_wp * (3.0_wp * t1 + t2)
      ct0 = 1.0_wp/16.0_wp * (-t1 * (1.0_wp - 2.0_wp * x1) + t2 * (1.0_wp + 2.0_wp * x2))
      ct1 = 1.0_wp/16.0_wp * (t2 - t1)
      cdj0 = -0.5_wp * t4 - 1.0_wp/4.0_wp * t4p  !-3.0_wp/4.0_wp * t4
      cdj1 = -1.0_wp/4.0_wp * t4p
      !write(*,*) c0r0, c0r1, c0s0, c0s1, c3r0, c3r1, c3s0, c3s1, cdr0, cdr1, ctau0, ctau1, cds0, cds1, ct0, ct1, cdj0, cdj1
!-----------------------------------------------------------------------
!        t0 and t3 terms
!-----------------------------------------------------------------------
      isoedens = 0.0_wp
! First isocalar part
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + (c0r0 + c3r0 * rhot**alpha) * rhot**2
            End Do
         End Do
      End Do
! Now the isovector part
      rho2vec = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               isorho = rho(ix, iy, iz, 1) - rho(ix, iy, iz, 2)
               isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + (c0r1 + c3r1 * rhot**alpha) * isorho**2
               rho2vec = rho2vec + wxyz(ix,iy,iz) * (c0r1 + c3r1 * rhot**alpha) * isorho**2
            End Do
         End Do
      End Do

! add the time-odd part
      If(itimrev == 0 .And. iodds == 1) Then
! isoscalar spin density contribution
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                  s2t = (spinden(ix, iy, iz, 1, 1)+spinden(ix, iy, iz, 1, 2)) ** 2 + (spinden(ix, iy, iz, 2, 1)+spinden(ix, iy, &
                 & iz, 2, 2)) ** 2 + (spinden(ix, iy, iz, 3, 1)+spinden(ix, iy, iz, 3, 2)) ** 2
                  isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + (c0s0 + c3s0*rhot**alpha) * s2t
               End Do
            End Do
         End Do
! Isovector part
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                  s2t = (spinden(ix, iy, iz, 1, 1)-spinden(ix, iy, iz, 1, 2)) ** 2 + (spinden(ix, iy, iz, 2, 1)-spinden(ix, iy, &
                 & iz, 2, 2)) ** 2 + (spinden(ix, iy, iz, 3, 1)-spinden(ix, iy, iz, 3, 2)) ** 2
                  isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + (c0s1 + c3s1*rhot**alpha) * s2t
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms depending on the laplacian of density (t1,t2)
!-----------------------------------------------------------------------
      workb(:, :, :, 1) = rho(:, :, :, 1) + rho(:, :, :, 2)
      workb(:, :, :, 2) = rho(:, :, :, 1) - rho(:, :, :, 2)
      Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, workb(:, :, :, 1), worka(:, :, :, 1))
      Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, workb(:, :, :, 2), worka(:, :, :, 2))
!
      delrhovec = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
! Isoscalar Part
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               d2rho = worka(ix, iy, iz, 1)
               isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + cdr0 * rhot * d2rho
! Isovector Part
               rhot = rho(ix, iy, iz, 1) - rho(ix, iy, iz, 2)
               d2rho = worka(ix, iy, iz, 2)
               isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + cdr1 * rhot * d2rho
               delrhovec = delrhovec + wxyz(ix,iy,iz) * cdr1 * rhot * d2rho
            End Do
         End Do
      End Do
!
! add the time-odd spin-density terms  - OMITTED DUE TO INSTABILITY
!      if(itimrev == 0 .and. iodds == 1) then
! calculate the laplacian of the spin-density vector (store in tmpvec)
      If(-1 > 0) Then ! this is to avoid tmpvec being flagged as unused by compiler
         Do iq = 1, 2
            Do imu = 1, 3
               Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, spinden(1, 1, 1, imu, iq), tmpvec(1, 1, 1, imu, iq))
            End Do
         End Do
      End If
! calculate the dot product s_q . del2 s_q and store in worka
!         worka = 0.0_wp
!         do iq = 1, 2
!            do imu = 1, 3
!               do iz = 1, ncolz
!                  do iy = 1, ncoly
!                     do ix = 1, ncolx
!                        worka(ix,iy,iz,iq) = worka(ix,iy,iz,iq) + spinden(ix,iy,iz,imu,iq)&
!                                                                * tmpvec(ix,iy,iz,imu,iq)
!                     end do
!                  end do
!               end do
!            end do
!         end do
!
!         t1a = (t2 + 3.0_wp*t1)/32.0_wp
!         t1b = (t2*x2 - 3.0_wp*t1*x1)/32.0_wp
!         do iz = 1, ncolz
!            do iy = 1, ncoly
!               do ix = 1, ncolx
!                  weight = wxyz(ix,iy,iz)
!                  sd2st = worka(ix,iy,iz,1) + worka(ix,iy,iz,2)&
!                        + spinden(ix,iy,iz,1,1)*tmpvec(ix,iy,iz,1,2)&
!                        + spinden(ix,iy,iz,2,1)*tmpvec(ix,iy,iz,2,2)&
!                        + spinden(ix,iy,iz,3,1)*tmpvec(ix,iy,iz,3,2)&
!                        + spinden(ix,iy,iz,1,2)*tmpvec(ix,iy,iz,1,1)&
!                        + spinden(ix,iy,iz,2,2)*tmpvec(ix,iy,iz,2,1)&
!                        + spinden(ix,iy,iz,3,2)*tmpvec(ix,iy,iz,3,1)
!                  sd2sn = worka(ix,iy,iz,1)
!                  sd2sp = worka(ix,iy,iz,2)
!                  ehf12_odd = ehf12_odd + weight*(t1a*(sd2sn + sd2sp) + t1b*sd2st)
!                  edens(ix,iy,iz) = edens(ix,iy,iz) + (t1a*(sd2sn + sd2sp) + t1b*sd2st)
!               end do
!            end do
!         end do
!      end if
!

      If(itimrev == itimrev) Then
!-------------------------------------------------------------------------------
!        square of the currents j_n and j_p (dot product of the current vector).
!        calculate the other t1,t2 term (rho*tau-j**2).
!        Note that lowercase j is not the spin-orbit current.
!-------------------------------------------------------------------------------
! subrtact -j**2 from rho*tau term (time-odd, traditionally always present)
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                  taut = tau(ix, iy, iz, 1) + tau(ix, iy, iz, 2)
                  xj2 = (currnt(ix, iy, iz, 1, 1)+currnt(ix, iy, iz, 1, 2)) ** 2 + (currnt(ix, iy, iz, 2, 1)+currnt(ix, iy, iz, &
                 & 2, 2)) ** 2 + (currnt(ix, iy, iz, 3, 1)+currnt(ix, iy, iz, 3, 2)) ** 2
                  isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + ctau0 * (rhot*taut - xj2)
               End Do
            End Do
         End Do
         ctauvec = 0.0_wp
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  rhot = rho(ix, iy, iz, 1) - rho(ix, iy, iz, 2)
                  taut = tau(ix, iy, iz, 1) - tau(ix, iy, iz, 2)
                  xj2 = (currnt(ix, iy, iz, 1, 1)-currnt(ix, iy, iz, 1, 2)) ** 2 + (currnt(ix, iy, iz, 2, 1)-currnt(ix, iy, iz, &
                 & 2, 2)) ** 2 + (currnt(ix, iy, iz, 3, 1)-currnt(ix, iy, iz, 3, 2)) ** 2
                  isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + ctau1 * (rhot*taut - xj2)
                  ctauvec = ctauvec + wxyz(ix,iy,iz) * ctau1 * (rhot*taut - xj2)
               End Do
            End Do
         End Do
      End If
!
      If(itheta == 1 .And. iodds == 1) Then
         If(itimrev == 0) Then
! construct the dot product s_q.T_q and store in worka
            worka = 0.0_wp
               Do imu = 1, 3
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           worka(ix, iy, iz, 1) = worka(ix, iy, iz, 1)+(spinden(ix, iy, iz, imu, 1)+spinden(ix, iy, iz, imu, 2))&
                           & * (kinvden(ix, iy, iz, imu, 1) + kinvden(ix, iy, iz, imu, 2))
                           worka(ix, iy, iz, 2) = worka(ix, iy, iz, 2)+(spinden(ix, iy, iz, imu, 1)-spinden(ix, iy, iz, imu, 2))&
                           & * (kinvden(ix, iy, iz, imu, 1) - kinvden(ix, iy, iz, imu, 2))
                        End Do
                     End Do
                  End Do
               End Do
! due to lack of work space do the s.T-J**2 contribution in two steps, first s.T

            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     stt = worka(ix, iy, iz, 1)
                     isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + ct0 * stt
                  End Do
               End Do
            End Do

            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     stt = worka(ix, iy, iz, 2)
                     isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + ct1 * stt
                  End Do
               End Do
            End Do
         End If
      End If
!
      If(itheta == 1) Then
! now construct the square of the J_mu_nu tensor (q'=q part) (not time-odd, flagged by itheta!)
         worka = 0.0_wp
            Do imu = 1, 3
               Do inu = 1, 3
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           worka(ix, iy, iz, 1) = worka(ix, iy, iz, 1) + &
                                (spincur(ix, iy, iz, imu, inu, 1) + spincur(ix, iy, iz, imu, inu, 2)) * &
                           &    (spincur(ix, iy, iz, imu, inu, 1) + spincur(ix, iy, iz, imu, inu, 2))
                           worka(ix, iy, iz, 2) = worka(ix, iy, iz, 2) + &
                                (spincur(ix, iy, iz, imu, inu, 1) - spincur(ix, iy, iz, imu, inu, 2)) * &
                           &    (spincur(ix, iy, iz, imu, inu, 1) - spincur(ix, iy, iz, imu, inu, 2))
                        End Do
                     End Do
                  End Do
               End Do
            End Do
! do the -J_q**2 contribution (q=q' part)
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  bjt = worka(ix, iy, iz, 1)
                  isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) - ct0 * bjt
               End Do
            End Do
         End Do
         jmunu2vec = 0.0_wp
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  bjt = worka(ix, iy, iz, 2)
                  isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) - ct1 * bjt
                  jmunu2vec = jmunu2vec - wxyz(ix,iy,iz) * ct1 * bjt
               End Do
            End Do
         End Do
         If(abs(jmunu2vec) < 1E-20) jmunu2vec = 0.0_wp
       End If

!-----------------------------------------------------------------------
!        spin-orbit contribution to hf energy
!-----------------------------------------------------------------------

      Call divergence(ncolx, ncoly, ncolz, der1x, der1y, der1z, &
      & (sodens(:, :, :, :, 1)+sodens(:,:,:,:, 2)), worka(:,:,:,1))
      Call divergence(ncolx, ncoly, ncolz, der1x, der1y, der1z, &
      & (sodens(:, :, :, :, 1)-sodens(:,:,:,:, 2)), worka(:,:,:,2))
!
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + cdj0 * rhot * worka(ix,iy,iz,1)
            End Do
         End Do
      End Do
      deljvec = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               rhot = rho(ix, iy, iz, 1) - rho(ix, iy, iz, 2)
               isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + cdj1 * rhot * worka(ix,iy,iz,2)
               deljvec = deljvec + wxyz(ix,iy,iz) * cdj1 * rhot * worka(ix,iy,iz,2)
            End Do
         End Do
      End Do
! time-odd term s_q. del x j_q, store the dot product in worka (q=q' part)
      If(itimrev == 0 .And. ioddls == 1) Then
         worka = 0.0_wp
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     workc(ix, iy, iz, 1, 1) = (delxj(ix, iy, iz, 2, 3, 1)-delxj(ix, iy, iz, 3, 2, 1)) &
                     & + (delxj(ix, iy, iz, 2, 3, 2)-delxj(ix, iy, iz, 3, 2, 2))
                     workc(ix, iy, iz, 2, 1) = - (delxj(ix, iy, iz, 1, 3, 1)-delxj(ix, iy, iz, 3, 1, 1)) &
                     & - (delxj(ix, iy, iz, 1, 3, 2)-delxj(ix, iy, iz, 3, 1, 2))
                     workc(ix, iy, iz, 3, 1) = (delxj(ix, iy, iz, 1, 2, 1)-delxj(ix, iy, iz, 2, 1, 1)) &
                     & + (delxj(ix, iy, iz, 1, 2, 2)-delxj(ix, iy, iz, 2, 1, 2))
                     workc(ix, iy, iz, 1, 2) = (delxj(ix, iy, iz, 2, 3, 1)-delxj(ix, iy, iz, 3, 2, 1)) &
                     & - (delxj(ix, iy, iz, 2, 3, 2)-delxj(ix, iy, iz, 3, 2, 2))
                     workc(ix, iy, iz, 2, 2) = - (delxj(ix, iy, iz, 1, 3, 1)-delxj(ix, iy, iz, 3, 1, 1)) &
                     & + (delxj(ix, iy, iz, 1, 3, 2)-delxj(ix, iy, iz, 3, 1, 2))
                     workc(ix, iy, iz, 3, 2) = (delxj(ix, iy, iz, 1, 2, 1)-delxj(ix, iy, iz, 2, 1, 1)) &
                     & - (delxj(ix, iy, iz, 1, 2, 2)-delxj(ix, iy, iz, 2, 1, 2))
                     worka(ix,iy,iz,1) = sum((spinden(ix,iy,iz,:,1)+spinden(ix,iy,iz,:,2))*workc(ix,iy,iz,:,1))
                     worka(ix,iy,iz,2) = sum((spinden(ix,iy,iz,:,1)-spinden(ix,iy,iz,:,2))*workc(ix,iy,iz,:,2))
                  End Do
               End Do
            End Do

! integrate
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  dj = worka(ix, iy, iz, 1)
                  isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + cdj0 * dj
               End Do
            End Do
         End Do
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  dj = worka(ix, iy, iz, 2)
                  isoedens(ix, iy, iz, 2) = isoedens(ix, iy, iz, 2) + cdj1 * dj
               End Do
            End Do
         End Do
       End If
!-----------------------------------------------------------------------
!        coulomb contribution with slater exchange
!        here we also calculate the correction term to koompan formula
!        for the density dependent coulomb exchange term.
!-----------------------------------------------------------------------
      ehfc = 0.0_wp
      ecorc = 0.0_wp
      cdens = 0.0_wp
      If(icoul /= 0) Then
         slate = - 3.0_wp / 4.0_wp * (3.0_wp/pi) ** (1.0_wp/3.0_wp) * e2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  rhop = rho(ix, iy, iz, 2)
                  wcr = wcoul(ix, iy, iz)
                  ehfc = ehfc + weight * (0.5_wp*rhop*wcr + slate*rhop**(4.0_wp/3.0_wp))
                  ecorc = ecorc + weight * slate / 3.0_wp * rhop**(4.0_wp/3.0_wp)
                  cdens(ix, iy, iz) = cdens(ix, iy, iz) + (0.5_wp*rhop*wcr + slate*rhop**(4.0_wp/3.0_wp))
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        Yukawa energy for BKN force
!-----------------------------------------------------------------------
            If(ayuk > 1.0e-5_wp) Then
                cdens(:, :, :) = cdens(:, :, :) + vyuk / 2.0_wp * wyuk(:, :, :)*(rho(:, :, :, 1) + rho(:, :, :, 2))
            End If
!-----------------------------------------------------------------------
!        kinetic energy contribution
!-----------------------------------------------------------------------
      ehft = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               weight = wxyz(ix, iy, iz)
               taun = tau(ix, iy, iz, 1)
               taup = tau(ix, iy, iz, 2)
               ehft = ehft + weight * (h2m(1)*taun+h2m(2)*taup)
               isoedens(ix, iy, iz, 1) = isoedens(ix, iy, iz, 1) + (h2m(1)*taun+h2m(2)*taup) + cdens(ix,iy,iz)
            End Do
         End Do
      End Do
!      write(8,'(a,3E15.6)') "rho*tau,j**2,diff = ",ehf22,ehflj2,(ehf22+ehflj2)
!      write(8,'(a,3E15.6)') "S.T,J**2,diff = ",ehf22_odd,ehfbj2,(ehf22_odd+ehfbj2)
!
      Return
End Subroutine isoenergy
