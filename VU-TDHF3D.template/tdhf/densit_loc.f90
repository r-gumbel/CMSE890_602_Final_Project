Subroutine densit_loc(lpsi, npsi, npmin, ncolx, ncoly, ncolz, vocc, der1x, der1y, der1z, psi, psin, pswk1, &
  &                   rhos, taus, drhos, currnts, itimrev)
!-----------------------------------------------------------------------
!        calculates quantities for the localization funtion
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In)   :: lpsi
      Integer, Intent(In)   :: ncolx
      Integer, Intent(In)   :: ncoly
      Integer, Intent(In)   :: ncolz
      Integer, Intent(In)   :: itimrev
      Integer, Intent(In)   :: npsi(2)
      Integer, Intent(In)   :: npmin(2)
      Real(wp), Intent(In)  :: vocc(lpsi, 2)
      Real(wp), Intent(Out) :: rhos(ncolx, ncoly, ncolz, 2, 2, 3)
      Real(wp), Intent(Out) :: taus(ncolx, ncoly, ncolz, 2, 2, 3)
      Real(wp), Intent(Out) :: drhos(ncolx, ncoly, ncolz, 3, 2, 2, 3)
      Real(wp), Intent(Out) :: currnts(ncolx, ncoly, ncolz, 3, 2, 2, 3)
      Real(wp), Intent(In)  :: der1x(ncolx, ncolx)
      Real(wp), Intent(In)  :: der1y(ncoly, ncoly)
      Real(wp), Intent(In)  :: der1z(ncolz, ncolz)
      Complex(wp), Intent(In) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp)             :: psin(ncolx, ncoly, ncolz, 2)
      Complex(wp)             :: pswk1(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: tau(ncolx, ncoly, ncolz, 2)
      Real(wp) :: spinden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: kinvden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp) :: spincur(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp) :: dsqmu(ncolx, ncoly, ncolz, 3, 2, 3)
      Real(wp) :: drho(ncolx, ncoly, ncolz, 3, 2)
      Integer  :: ix, iy, iz, is, iquant, izero, isign, isp, iq, nst
      Real(wp) :: weight
!-----------------------------------------------------------------------
!        sum over all states
!-----------------------------------------------------------------------
      izero = 1.0_wp
      if(itimrev == 1) izero = 0.0_wp
!
      rho = 0.0_wp
      tau = 0.0_wp
      currnt = 0.0_wp
      spinden = 0.0_wp
      kinvden = 0.0_wp
      spincur = 0.0_wp
      dsqmu = 0.0_wp
      drho = 0.0_wp
      Do iq = 1, 2
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(nst,psin,pswk1,weight) REDUCTION(+:rho,tau,currnt,spinden,kinvden,spincur,dsqmu,drho)
!$OMP DO
         Do nst = npmin(iq), npsi(iq)
            weight = vocc(nst, iq)
            psin(:ncolx, :ncoly, :ncolz, :2) = psi(:ncolx, :ncoly, :ncolz, :2, nst, iq)
!-----------------------------------------------------------------------
!        calculate the neutron and proton densities
!-----------------------------------------------------------------------
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     Do is = 1,  2
                        rho(ix, iy, iz, iq) = rho(ix, iy, iz, iq) + weight*psin(ix, iy, iz, is)*conjg(psin(ix, iy, iz, is))
                     End Do
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        calculate the time-odd spin density vector, s
!-----------------------------------------------------------------------
            If(itimrev == 0) Then
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        spinden(ix, iy, iz, 1, iq) = spinden(ix, iy, iz, 1, iq) + 2.0_wp*weight*Real(conjg(psin(ix, iy, iz, 1))* &
                     &                               psin(ix, iy, iz, 2))
                        spinden(ix, iy, iz, 2, iq) = spinden(ix, iy, iz, 2, iq) + 2.0_wp*weight*Dimag(conjg(psin(ix, iy, iz, 1))* &
                     &                               psin(ix, iy, iz, 2))
                        spinden(ix, iy, iz, 3, iq) = spinden(ix, iy, iz, 3, iq) + weight*(conjg(psin(ix, iy, iz, 1))*psin(ix, iy, iz,1)-&
                     &                               conjg(psin(ix, iy, iz, 2))*psin(ix, iy, iz, 2))
                     End Do
                  End Do
               End Do
            End If
!-----------------------------------------------------------------------
!        calculate all densities which require d/dx psi
!-----------------------------------------------------------------------
            Call cmulx(ncolx, ncoly, ncolz, der1x, psin, pswk1, 0)
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     Do is = 1, 2
                        isp = 3 - is
                        isign = 3 - 2*is
                        tau(ix, iy, iz, iq) = tau(ix, iy, iz, iq) + weight*pswk1(ix, iy, iz, is)*conjg(pswk1(ix, iy, iz, is))
                        drho(ix,iy,iz,1,iq) = drho(ix,iy,iz,1,iq) + 2.0_wp*weight*(Real(conjg(psin(ix, iy, iz, is))* &
                    &   pswk1(ix, iy, iz, is)))
                        dsqmu(ix,iy,iz,1,iq,1) = dsqmu(ix,iy,iz,1,iq,1) + 2.0_wp*weight*izero*( &
                        Real(conjg(psin(ix, iy, iz, isp))*pswk1(ix, iy, iz, is)))
                        dsqmu(ix,iy,iz,1,iq,2) = dsqmu(ix,iy,iz,1,iq,2) + 2.0_wp*weight*izero*( &
                    &   isign*Dimag(pswk1(ix, iy, iz, isp)*conjg(psin(ix, iy, iz, is))))
                        dsqmu(ix,iy,iz,1,iq,3) = dsqmu(ix,iy,iz,1,iq,3) + 2.0_wp*weight*izero*( &
                    &   isign*Real(conjg(psin(ix, iy, iz, is))*pswk1(ix, iy, iz, is)))
                     End Do
                  End Do
               End Do
            End Do
! time-odd quantities: current j
            If(itimrev == 0) Then
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        currnt(ix, iy, iz, 1, iq) = currnt(ix, iy, iz, 1, iq) + weight*Dimag(conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy,&
                                        & iz, 1)+conjg(psin(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                     End Do
                  End Do
               End Do
! kinetic energy density (vector part): T
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        kinvden(ix, iy, iz, 1, iq) = kinvden(ix, iy, iz, 1, iq) + 2.0_wp * weight * Real(conjg(pswk1(ix, iy, iz, &
                       & 1))*pswk1(ix, iy, iz, 2))
                        kinvden(ix, iy, iz, 2, iq) = kinvden(ix, iy, iz, 2, iq) + 2.0_wp * weight * Dimag(conjg(pswk1(ix, iy, iz, &
                       & 1))*pswk1(ix, iy, iz, 2))
                        kinvden(ix, iy, iz, 3, iq) = kinvden(ix, iy, iz, 3, iq) + weight * (conjg(pswk1(ix, iy, iz, 1))*pswk1(ix, &
                       & iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                     End Do
                   End Do
               End Do
            End If
! spin-orbit current tensor J_mu_nu
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     spincur(ix, iy, iz, 1, 1, iq) = spincur(ix, iy, iz, 1, 1, iq) - weight * Dimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*psin(ix, iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
                     spincur(ix, iy, iz, 1, 2, iq) = spincur(ix, iy, iz, 1, 2, iq) + weight * Real(conjg(pswk1(ix, iy, iz, 1))*psin(ix,&
                    & iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
                     spincur(ix, iy, iz, 1, 3, iq) = spincur(ix, iy, iz, 1, 3, iq) - weight * Dimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*psin(ix, iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*psin(ix, iy, iz, 2))
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        add all densities which require d/dy psi
!-----------------------------------------------------------------------
            Call cmuly(ncolx, ncoly, ncolz, der1y, psin, pswk1, 0)
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     Do is = 1, 2
                       isp = 3 - is
                       isign = 3 - 2*is
                       tau(ix, iy, iz, iq) = tau(ix, iy, iz, iq) + weight*pswk1(ix, iy, iz, is)*conjg(pswk1(ix, iy, iz, is))
                       drho(ix,iy,iz,2,iq) = drho(ix,iy,iz,2,iq) + 2.0_wp*weight*(Real(conjg(psin(ix, iy, iz, is))* &
                   &   pswk1(ix, iy, iz, is)))
                       dsqmu(ix,iy,iz,2,iq,1) = dsqmu(ix,iy,iz,2,iq,1) + 2.0_wp*weight*izero*( &
                       Real(conjg(psin(ix, iy, iz, isp))*pswk1(ix, iy, iz, is)))
                       dsqmu(ix,iy,iz,2,iq,2) = dsqmu(ix,iy,iz,2,iq,2) + 2.0_wp*weight*izero*( &
                   &   isign*Dimag(pswk1(ix, iy, iz, isp)*conjg(psin(ix, iy, iz, is))))
                       dsqmu(ix,iy,iz,2,iq,3) = dsqmu(ix,iy,iz,2,iq,3) + 2.0_wp*weight*izero*( &
                   &   isign*Real(conjg(psin(ix, iy, iz, is))*pswk1(ix, iy, iz, is)))
                     End Do
                  End Do
               End Do
            End Do
! time-odd quantities: current j
            If(itimrev == 0) Then
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        currnt(ix, iy, iz, 2, iq) = currnt(ix, iy, iz, 2, iq) + weight * Dimag(conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy,&
                       & iz, 1)+conjg(psin(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                     End Do
                  End Do
               End Do
! kinetic energy density T (vector-part)
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        kinvden(ix, iy, iz, 1, iq) = kinvden(ix, iy, iz, 1, iq) + 2.0_wp * weight * Real(conjg(pswk1(ix, iy, iz, &
                      & 1))*pswk1(ix, iy, iz, 2))
                        kinvden(ix, iy, iz, 2, iq) = kinvden(ix, iy, iz, 2, iq) + 2.0_wp * weight * Dimag(conjg(pswk1(ix, iy, iz, &
                      & 1))*pswk1(ix, iy, iz, 2))
                        kinvden(ix, iy, iz, 3, iq) = kinvden(ix, iy, iz, 3, iq) + weight * (conjg(pswk1(ix, iy, iz, 1))*pswk1(ix, &
                      & iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                     End Do
                   End Do
               End Do
            End If
! spin-orbit current tensor J_mu_nu
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     spincur(ix, iy, iz, 2, 1, iq) = spincur(ix, iy, iz, 2, 1, iq) - weight * Dimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*psin(ix, iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
                     spincur(ix, iy, iz, 2, 2, iq) = spincur(ix, iy, iz, 2, 2, iq) + weight * Real(conjg(pswk1(ix, iy, iz, 1))*psin(ix, &
                    & iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
                     spincur(ix, iy, iz, 2, 3, iq) = spincur(ix, iy, iz, 2, 3, iq) - weight * Dimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*psin(ix, iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*psin(ix, iy, iz, 2))
                  End Do
               End Do
            End Do
!-----------------------------------------------------------------------
!        add all densities which require d/dz psi
!-----------------------------------------------------------------------
            Call cmulz(ncolx, ncoly, ncolz, der1z, psin, pswk1, 0)
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     Do is = 1, 2
                       isp = 3 - is
                       isign = 3 - 2*is
                       tau(ix, iy, iz, iq) = tau(ix, iy, iz, iq) + weight*pswk1(ix, iy, iz, is)*conjg(pswk1(ix, iy, iz, is))
                       drho(ix,iy,iz,3,iq) = drho(ix,iy,iz,3,iq) + 2.0_wp*weight*(Real(conjg(psin(ix, iy, iz, is))* &
                   &   pswk1(ix, iy, iz, is)))
                       dsqmu(ix,iy,iz,3,iq,1) = dsqmu(ix,iy,iz,3,iq,1) + 2.0_wp*weight*izero*( &
                       Real(conjg(psin(ix, iy, iz, isp))*pswk1(ix, iy, iz, is)))
                       dsqmu(ix,iy,iz,3,iq,2) = dsqmu(ix,iy,iz,3,iq,2) + 2.0_wp*weight*izero*( &
                   &   isign*Dimag(pswk1(ix, iy, iz, isp)*conjg(psin(ix, iy, iz, is))))
                       dsqmu(ix,iy,iz,3,iq,3) = dsqmu(ix,iy,iz,3,iq,3) + 2.0_wp*weight*izero*( &
                   &   isign*Real(conjg(psin(ix, iy, iz, is))*pswk1(ix, iy, iz, is)))
                     End Do
                  End Do
               End Do
            End Do
! time-odd quantities: current j
            If(itimrev == 0) Then
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        currnt(ix, iy, iz, 3, iq) = currnt(ix, iy, iz, 3, iq) + weight * Dimag(conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, &
                       & iz, 1)+conjg(psin(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                     End Do
                  End Do
               End Do
! kinetic energy density (vector-part)
               Do ix = 1, ncolx
                  Do iy = 1, ncoly
                     Do iz = 1, ncolz
                        kinvden(ix, iy, iz, 1, iq) = kinvden(ix, iy, iz, 1, iq) + 2.0_wp * weight * Real(conjg(pswk1(ix, iy, iz, &
                      & 1))*pswk1(ix, iy, iz, 2))
                        kinvden(ix, iy, iz, 2, iq) = kinvden(ix, iy, iz, 2, iq) + 2.0_wp * weight * Dimag(conjg(pswk1(ix, iy, iz, &
                      & 1))*pswk1(ix, iy, iz, 2))
                        kinvden(ix, iy, iz, 3, iq) = kinvden(ix, iy, iz, 3, iq) + weight * (conjg(pswk1(ix, iy, iz, 1))*pswk1(ix, &
                      & iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                     End Do
                   End Do
               End Do
            End If
! spin-orbit current tensor J_mu_nu
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     spincur(ix, iy, iz, 3, 1, iq) = spincur(ix, iy, iz, 3, 1, iq) - weight * Dimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*psin(ix, iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
                     spincur(ix, iy, iz, 3, 2, iq) = spincur(ix, iy, iz, 3, 2, iq) + weight*Real(conjg(pswk1(ix, iy, iz, 1))*psin(ix, &
                    & iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
                     spincur(ix, iy, iz, 3, 3, iq) = spincur(ix, iy, iz, 3, 3, iq) - weight * Dimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*psin(ix, iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*psin(ix, iy, iz, 2))
                  End Do
               End Do
            End Do
!
         End Do
!$OMP END DO
!$OMP END PARALLEL
      End Do
!-----------------------------------------------------------------------
!        calculate localization densities and currents
!-----------------------------------------------------------------------
        Do iq = 1, 2
           Do ix = 1, ncolx
              Do iy = 1, ncoly
                 Do iz = 1, ncolz
                    Do is = 1,  2
                       isign = 3 - 2*is
                       Do iquant = 1, 3
                          rhos(ix,iy,iz,is,iq,iquant) = 0.5_wp*rho(ix,iy,iz,iq) + isign*spinden(ix,iy,iz,iquant,iq)
                          taus(ix,iy,iz,is,iq,iquant) = 0.5_wp*tau(ix,iy,iz,iq) + isign*kinvden(ix,iy,iz,iquant,iq)
                          drhos(ix,iy,iz,1,is,iq,iquant) = 0.5_wp*drho(ix,iy,iz,1,iq) + isign*dsqmu(ix,iy,iz,1,iq,iquant)!
                          drhos(ix,iy,iz,2,is,iq,iquant) = 0.5_wp*drho(ix,iy,iz,2,iq) + isign*dsqmu(ix,iy,iz,2,iq,iquant)
                          drhos(ix,iy,iz,3,is,iq,iquant) = 0.5_wp*drho(ix,iy,iz,3,iq) + isign*dsqmu(ix,iy,iz,3,iq,iquant)
                          currnts(ix,iy,iz,1,is,iq,iquant) = 0.5_wp*currnt(ix,iy,iz,1,iq) + isign*spincur(ix,iy,iz,iquant,1,iq)
                          currnts(ix,iy,iz,2,is,iq,iquant) = 0.5_wp*currnt(ix,iy,iz,2,iq) + isign*spincur(ix,iy,iz,iquant,2,iq)
                          currnts(ix,iy,iz,3,is,iq,iquant) = 0.5_wp*currnt(ix,iy,iz,3,iq) + isign*spincur(ix,iy,iz,iquant,3,iq)
                       End Do
                    End Do
                 End Do
              End Do
           End Do
        End Do
!
      Return
End Subroutine densit_loc
