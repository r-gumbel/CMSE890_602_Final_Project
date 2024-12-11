Subroutine energy(ncolx, ncoly, ncolz, h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, wxyz, der1x, der1y, &
& der1z, der2x, der2y, der2z, rho, tau, currnt, sodens, spinden, kinvden, spincur, delxj, tmpvec, wcoul, wyuk, edens, worka, &
& icoul, ehft, ehf0, ehf12, ehf22, ehf3, ehfls, ehf0_odd, ehf12_odd, ehf22_odd, ehf3_odd, ehfls_odd, ehflj2, ehfbj2, ehfc, ecorc, &
& ehfy, itheta, itimrev, iodds)
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
      Real(wp), Intent(Out) :: edens(ncolx, ncoly, ncolz)
      Real(wp), Intent(Out) :: ehft
      Real(wp), Intent(Out) :: ehf0
      Real(wp), Intent(Out) :: ehf12
      Real(wp), Intent(Out) :: ehf22
      Real(wp), Intent(Out) :: ehf3
      Real(wp), Intent(Out) :: ehfls
      Real(wp), Intent(Out) :: ehfc
      Real(wp), Intent(Out) :: ecorc
      Real(wp), Intent(Out) :: ehfy
      Real(wp), Intent(Out) :: ehf0_odd
      Real(wp), Intent(Out) :: ehf12_odd
      Real(wp), Intent(Out) :: ehf22_odd
      Real(wp), Intent(Out) :: ehf3_odd
      Real(wp), Intent(Out) :: ehfls_odd
      Real(wp), Intent(Out) :: ehfbj2
      Real(wp), Intent(Out) :: ehflj2
      Real(wp) :: worka(ncolx, ncoly, ncolz, 2)
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
      Integer :: iz, iy, ix, iq, imu, inu, ic
      Real(wp) :: t0a, t0b, t3a, t3b, weight, rhot, rhon, rhop, t1a, t1b, d2rho, d2rhon, d2rhop, taut, taun, taup, xj2, xj2n, &
     & xj2p, t4a, t4b, dj, djn, djp, slate, wcr, s2t, s2n, s2p, stt, stp, stn, bjt, bjn, bjp
!      real(wp) ::  sd2st, sd2sn, sd2sp
!
      edens = 0.0_wp
!-----------------------------------------------------------------------
!        t0 and t3 terms
!-----------------------------------------------------------------------
      t0a = t0 / 2.0_wp * (1.0_wp+0.5_wp*x0)
      t0b = t0 / 2.0_wp * (0.5_wp+x0)
      t3a = t3 / 12.0_wp * (1.0_wp+0.5_wp*x3)
      t3b = t3 / 12.0_wp * (0.5_wp+x3)
!
      ehf0 = 0.0_wp
      ehf3 = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               weight = wxyz(ix, iy, iz)
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               rhon = rho(ix, iy, iz, 1)
               rhop = rho(ix, iy, iz, 2)
               ehf0 = ehf0 + weight * (t0a*rhot**2-t0b*(rhop**2+rhon**2))
               ehf3 = ehf3 + weight * rhot ** alpha * (t3a*rhot**2-t3b*(rhop**2+rhon**2))
               edens(ix, iy, iz) = edens(ix, iy, iz) + (t0a*rhot**2-t0b*(rhop**2+rhon**2)) + rhot ** alpha * &
              & (t3a*rhot**2-t3b*(rhop**2+rhon**2))
            End Do
         End Do
      End Do
! add the time-odd part
      ehf0_odd = 0.0_wp
      ehf3_odd = 0.0_wp
      If(itimrev == 0 .And. iodds == 1) Then
! square the spin-density and use worka as temporary storage
         Do iq = 1, 2
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     worka(ix, iy, iz, iq) = spinden(ix, iy, iz, 1, iq) ** 2 + spinden(ix, iy, iz, 2, iq) ** 2 + spinden(ix, iy, &
                    & iz, 3, iq) ** 2
                  End Do
               End Do
            End Do
         End Do
!
         t0a = t0 * x0 / 4.0_wp
         t0b = t0 / 4.0_wp
         t3a = t3 * x3 / 24.0_wp
         t3b = t3 / 24.0_wp
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
                  s2t = (spinden(ix, iy, iz, 1, 1)+spinden(ix, iy, iz, 1, 2)) ** 2 + (spinden(ix, iy, iz, 2, 1)+spinden(ix, iy, &
                 & iz, 2, 2)) ** 2 + (spinden(ix, iy, iz, 3, 1)+spinden(ix, iy, iz, 3, 2)) ** 2
                  s2n = worka(ix, iy, iz, 1)
                  s2p = worka(ix, iy, iz, 2)
                  ehf0_odd = ehf0_odd + weight * (t0a*s2t-t0b*(s2n+s2p))
                  ehf3_odd = ehf3_odd + weight * rhot ** alpha * (t3a*s2t-t3b*(s2n+s2p))
                  edens(ix, iy, iz) = edens(ix, iy, iz) + (t0a*s2t-t0b*(s2n+s2p)) + rhot ** alpha * (t3a*s2t-t3b*(s2n+s2p))
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms depending on the laplacian of density (t1,t2)
!-----------------------------------------------------------------------
      Do iq = 1, 2
         Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, rho(1, 1, 1, iq), worka(1, 1, 1, iq))
      End Do
!
      t1a = (t2*(1.0_wp+x2/2.0_wp)-3.0_wp*t1*(1.0_wp+x1/2.0_wp)) / 16.0_wp
      t1b = (3.0_wp*t1*(0.5_wp+x1)+t2*(0.5_wp+x2)) / 16.0_wp
      ehf12 = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               weight = wxyz(ix, iy, iz)
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               rhon = rho(ix, iy, iz, 1)
               rhop = rho(ix, iy, iz, 2)
               d2rho = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2)
               d2rhon = worka(ix, iy, iz, 1)
               d2rhop = worka(ix, iy, iz, 2)
               ehf12 = ehf12 + weight * (t1a*rhot*d2rho+t1b*(rhop*d2rhop+rhon*d2rhon))
               edens(ix, iy, iz) = edens(ix, iy, iz) + (t1a*rhot*d2rho+t1b*(rhop*d2rhop+rhon*d2rhon))
            End Do
         End Do
      End Do
! add the time-odd spin-density terms  - OMITTED DUE TO INSTABILITY
      ehf12_odd = 0.0_wp
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
      t1a = (t1*(1.0_wp+x1/2.0_wp)+t2*(1.0_wp+x2/2.0_wp)) / 4.0_wp
      t1b = (t2*(x2+0.5_wp)-t1*(x1+0.5_wp)) / 4.0_wp
      ehf22 = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               weight = wxyz(ix, iy, iz)
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               rhon = rho(ix, iy, iz, 1)
               rhop = rho(ix, iy, iz, 2)
               taut = tau(ix, iy, iz, 1) + tau(ix, iy, iz, 2)
               taun = tau(ix, iy, iz, 1)
               taup = tau(ix, iy, iz, 2)
               ehf22 = ehf22 + weight * (t1a*rhot*taut+t1b*(rhon*taun+rhop*taup))
               edens(ix, iy, iz) = edens(ix, iy, iz) + (t1a*rhot*taut+t1b*(rhon*taun+rhop*taup))
            End Do
         End Do
      End Do
!
      ehflj2 = 0.0_wp
      If(itimrev == 0) Then
!-------------------------------------------------------------------------------
!        square of the currents j_n and j_p (dot product of the current vector).
!        calculate the other t1,t2 term (rho*tau-j**2).
!        Note that lowercase j is not the spin-orbit current.
!-------------------------------------------------------------------------------
         Do iq = 1, 2
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     worka(ix, iy, iz, iq) = currnt(ix, iy, iz, 1, iq) ** 2 + currnt(ix, iy, iz, 2, iq) ** 2 + currnt(ix, iy, iz, &
                    & 3, iq) ** 2
                  End Do
               End Do
            End Do
         End Do
! subrtact -j**2 from rho*tau term (time-odd, traditionally always present)
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  xj2 = (currnt(ix, iy, iz, 1, 1)+currnt(ix, iy, iz, 1, 2)) ** 2 + (currnt(ix, iy, iz, 2, 1)+currnt(ix, iy, iz, &
                 & 2, 2)) ** 2 + (currnt(ix, iy, iz, 3, 1)+currnt(ix, iy, iz, 3, 2)) ** 2
                  xj2n = worka(ix, iy, iz, 1)
                  xj2p = worka(ix, iy, iz, 2)
                  ehflj2 = ehflj2 + weight * (-t1a*xj2-t1b*(xj2n+xj2p))
                  edens(ix, iy, iz) = edens(ix, iy, iz) + (-t1a*xj2-t1b*(xj2n+xj2p))
               End Do
            End Do
         End Do
      End If
!
      ehf22_odd = 0.0_wp
      If(itheta == 1 .And. iodds == 1) Then
         If(itimrev == 0) Then
! construct the dot product s_q.T_q and store in worka
            worka = 0.0_wp
            Do iq = 1, 2
               Do imu = 1, 3
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           worka(ix, iy, iz, iq) = worka(ix, iy, iz, iq) + spinden(ix, iy, iz, imu, iq) * kinvden(ix, iy, iz, &
                          & imu, iq)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
! due to lack of work space do the s.T-J**2 contribution in two steps, first s.T
            t1a = (t1*x1+t2*x2) / 8.0_wp
            t1b = (t2-t1) / 8.0_wp
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     weight = wxyz(ix, iy, iz)
                     stt = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2) + spinden(ix, iy, iz, 1, 1) * kinvden(ix, iy, iz, 1, 2) + &
                    & spinden(ix, iy, iz, 2, 1) * kinvden(ix, iy, iz, 2, 2) + spinden(ix, iy, iz, 3, 1) * kinvden(ix, iy, iz, 3, &
                    & 2) + spinden(ix, iy, iz, 1, 2) * kinvden(ix, iy, iz, 1, 1) + spinden(ix, iy, iz, 2, 2) * kinvden(ix, iy, &
                    & iz, 2, 1) + spinden(ix, iy, iz, 3, 2) * kinvden(ix, iy, iz, 3, 1)
                     stn = worka(ix, iy, iz, 1)
                     stp = worka(ix, iy, iz, 2)
                     ehf22_odd = ehf22_odd + weight * (t1a*stt+t1b*(stn+stp))
                     edens(ix, iy, iz) = edens(ix, iy, iz) + (t1a*stt+t1b*(stn+stp))
                  End Do
               End Do
            End Do
         End If
      End If
!
      ehfbj2 = 0.0_wp
      If(itheta == 1) Then
! now construct the square of the J_mu_nu tensor (q'=q part) (not time-odd, flagged by itheta!)
         worka = 0.0_wp
         Do iq = 1, 2
            Do imu = 1, 3
               Do inu = 1, 3
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           worka(ix, iy, iz, iq) = worka(ix, iy, iz, iq) + spincur(ix, iy, iz, imu, inu, iq) * spincur(ix, iy, &
                          & iz, imu, inu, iq)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
! do the -J_q**2 contribution (q=q' part)
         t1a = - (t1*x1+t2*x2) / 8.0_wp
         t1b = - (t2-t1) / 8.0_wp
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  bjt = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2)
                  bjn = worka(ix, iy, iz, 1)
                  bjp = worka(ix, iy, iz, 2)
                  ehfbj2 = ehfbj2 + weight * (t1a*bjt+t1b*(bjn+bjp))
                  edens(ix, iy, iz) = edens(ix, iy, iz) + (t1a*bjt+t1b*(bjn+bjp))
               End Do
            End Do
         End Do
! construct the J_q.J_q' part
         worka = 0.0_wp
         Do iq = 1, 2
            ic = 3 - iq
            Do imu = 1, 3
               Do inu = 1, 3
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           worka(ix, iy, iz, iq) = worka(ix, iy, iz, iq) + spincur(ix, iy, iz, imu, inu, iq) * spincur(ix, iy, &
                          & iz, imu, inu, ic)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  bjt = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2)
                  ehfbj2 = ehfbj2 + weight * t1a * bjt
                  edens(ix, iy, iz) = edens(ix, iy, iz) + t1a * bjt
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        spin-orbit contribution to hf energy
!-----------------------------------------------------------------------
      Do iq = 1, 2
         Call divergence(ncolx, ncoly, ncolz, der1x, der1y, der1z, sodens(1, 1, 1, 1, iq), worka(1, 1, 1, iq))
      End Do
!
      t4a = - t4 / 2.0_wp
      t4b = - t4p / 2.0_wp
      ehfls = 0.0_wp
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               weight = wxyz(ix, iy, iz)
               rhot = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
               rhon = rho(ix, iy, iz, 1)
               rhop = rho(ix, iy, iz, 2)
               dj = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2)
               djn = worka(ix, iy, iz, 1)
               djp = worka(ix, iy, iz, 2)
               ehfls = ehfls + weight * (t4a*rhot*dj+t4b*(rhop*djp+rhon*djn))
               edens(ix, iy, iz) = edens(ix, iy, iz) + t4a * rhot * dj + t4b * (rhop*djp+rhon*djn)
            End Do
         End Do
      End Do
! time-odd term s_q. del x j_q, store the dot product in worka (q=q' part)
      ehfls_odd = 0.0_wp
      If(itimrev == 0 .And. ioddls == 1) Then
         worka = 0.0_wp
         Do iq = 1, 2
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     worka(ix, iy, iz, iq) = worka(ix, iy, iz, iq) + spinden(ix, iy, iz, 1, iq) * (delxj(ix, iy, iz, 2, 3, &
                    & iq)-delxj(ix, iy, iz, 3, 2, iq)) - spinden(ix, iy, iz, 2, iq) * (delxj(ix, iy, iz, 1, 3, iq)-delxj(ix, iy, &
                    & iz, 3, 1, iq)) + spinden(ix, iy, iz, 3, iq) * (delxj(ix, iy, iz, 1, 2, iq)-delxj(ix, iy, iz, 2, 1, iq))
                  End Do
               End Do
            End Do
         End Do
! integrate
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  dj = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2)
                  djn = worka(ix, iy, iz, 1)
                  djp = worka(ix, iy, iz, 2)
                  ehfls_odd = ehfls_odd + weight * (t4a*dj+t4b*(djn+djp))
                  edens(ix, iy, iz) = edens(ix, iy, iz) + t4a * dj + t4b * (djn+djp)
               End Do
            End Do
         End Do
! time-odd term s_q. del x j_q', store the dot product in worka
         worka = 0.0_wp
         Do iq = 1, 2
            ic = 3 - iq
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     worka(ix, iy, iz, iq) = worka(ix, iy, iz, iq) + spinden(ix, iy, iz, 1, iq) * (delxj(ix, iy, iz, 2, 3, &
                    & ic)-delxj(ix, iy, iz, 3, 2, ic)) - spinden(ix, iy, iz, 2, iq) * (delxj(ix, iy, iz, 1, 3, ic)-delxj(ix, iy, &
                    & iz, 3, 1, ic)) + spinden(ix, iy, iz, 3, iq) * (delxj(ix, iy, iz, 1, 2, ic)-delxj(ix, iy, iz, 2, 1, ic))
                  End Do
               End Do
            End Do
         End Do
! integrate
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  dj = worka(ix, iy, iz, 1) + worka(ix, iy, iz, 2)
                  ehfls_odd = ehfls_odd + weight * t4a * dj
                  edens(ix, iy, iz) = edens(ix, iy, iz) + t4a * dj
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
      If(icoul /= 0) Then
         slate = - 3.0_wp / 4.0_wp * (3.0_wp/pi) ** (1.0_wp/3.0_wp) * e2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  weight = wxyz(ix, iy, iz)
                  rhop = rho(ix, iy, iz, 2)
                  wcr = wcoul(ix, iy, iz)
                  ehfc = ehfc + weight * (0.5_wp*rhop*wcr+slate*rhop**(4.0_wp/3.0_wp))
                  ecorc = ecorc + weight * slate / 3.0_wp * rhop ** (4.0_wp/3.0_wp)
                  edens(ix, iy, iz) = edens(ix, iy, iz) + (0.5_wp*rhop*wcr+slate*rhop**(4.0_wp/3.0_wp))
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        Yukawa energy for BKN force
!-----------------------------------------------------------------------
      ehfy = 0.0_wp
      If(ayuk > 1.0e-5_wp) Then
         ehfy = vyuk / 2.0_wp * Sum(wxyz*wyuk(:, :, :)*(rho(:, :, :, 1)+rho(:, :, :, 2)))
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
               edens(ix, iy, iz) = edens(ix, iy, iz) + (h2m(1)*taun+h2m(2)*taup)
            End Do
         End Do
      End Do
!      write(8,'(a,3E15.6)') "rho*tau,j**2,diff = ",ehf22,ehflj2,(ehf22+ehflj2)
!      write(8,'(a,3E15.6)') "S.T,J**2,diff = ",ehf22_odd,ehfbj2,(ehf22_odd+ehfbj2)
!
      Return
End Subroutine energy
