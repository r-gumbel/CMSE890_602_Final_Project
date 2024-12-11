Subroutine densit(lpsi, ncolx, ncoly, ncolz, nst, iq, vocc, rho, tau, currnt, sodens, spinden, kinvden, spincur, der1x, der1y, &
& der1z, psin, pswk1, rhos, taus, drho, currnts, itimrev, itheta, iodds)
!
!-----------------------------------------------------------------------
!        calculates quantities for the hamiltonian density
!
!        rho     :   proton and neutron densities
!        tau     :   kinetic energy densities
!        currnt  :   currents in the local part of h
!        sodens  :   spin-orbit currents
!        spinden :   spin-density s
!        kinvden :   kinetic energy density (vector part)
!        spincur :   spin-current tensor
!
!        note that in addition to isospin index currnt and various densities
!        carry the vector-component index in cartesian coordinates
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: nst
      Integer, Intent(In) :: iq
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: itheta
      Integer, Intent(In) :: iodds
      Real(wp), Intent(In) :: vocc(lpsi, 2)
      Real(wp), Intent(Inout) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(Out)   :: rhos(ncolx, ncoly, ncolz, 2, 2)
      Real(wp), Intent(Inout) :: tau(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(Out)   :: taus(ncolx, ncoly, ncolz, 2, 2)
      Real(wp), Intent(Inout) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(Out)   :: drho(ncolx, ncoly, ncolz, 3, 2, 2)
      Real(wp), Intent(Out)   :: currnts(ncolx, ncoly, ncolz, 3, 2, 2)
      Real(wp), Intent(Inout) :: sodens(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(Inout) :: spinden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(Inout) :: kinvden(ncolx, ncoly, ncolz, 3, 2)
      Real(wp), Intent(Inout) :: spincur(ncolx, ncoly, ncolz, 3, 3, 2)
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Complex(wp), Intent(In) :: psin(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk1(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Integer, Parameter :: ioddls = 1 ! should be 1 except special circumstances
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ix, iy, iz, is
      Real(wp) :: weight
!-----------------------------------------------------------------------
!        everything is done for one particular state
!-----------------------------------------------------------------------
      weight = vocc(nst, iq)
!-----------------------------------------------------------------------
!        calculate the neutron and proton densities
!-----------------------------------------------------------------------
      Do ix = 1, ncolx
         Do iy = 1, ncoly
            Do iz = 1, ncolz
               Do is = 1,  2
                  rho(ix, iy, iz, iq) = rho(ix, iy, iz, iq) + weight*psin(ix, iy, iz, is)*conjg(psin(ix, iy, iz, is))
                  rhos(ix,iy,iz,is,iq) = rhos(ix,iy,iz,is,iq) + weight*psin(ix,iy,iz,is)*conjg(psin(ix,iy,iz,is))
               End Do
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        calculate the time-odd spin density vector, s
!-----------------------------------------------------------------------
      If(itimrev == 0 .And. ioddls == 1) Then
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  spinden(ix, iy, iz, 1, iq) = spinden(ix, iy, iz, 1, iq) + 2.0_wp * weight * real(conjg(psin(ix, iy, iz, &
                 & 1))*psin(ix, iy, iz, 2))
                  spinden(ix, iy, iz, 2, iq) = spinden(ix, iy, iz, 2, iq) + 2.0_wp * weight * aimag(conjg(psin(ix, iy, iz, &
                 & 1))*psin(ix, iy, iz, 2))
                  spinden(ix, iy, iz, 3, iq) = spinden(ix, iy, iz, 3, iq) + weight * (conjg(psin(ix, iy, iz, 1))*psin(ix, iy, iz,&
                 & 1)-conjg(psin(ix, iy, iz, 2))*psin(ix, iy, iz, 2))
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
                  tau(ix, iy, iz, iq) = tau(ix, iy, iz, iq) + weight*pswk1(ix, iy, iz, is)*conjg(pswk1(ix, iy, iz, is))
                  taus(ix,iy,iz,is,iq) = taus(ix,iy,iz,is,iq) + weight*pswk1(ix,iy,iz,is)*conjg(pswk1(ix,iy,iz,is))
                  drho(ix,iy,iz,1,is,iq) = drho(ix,iy,iz,1,is,iq) + weight*2.0_wp*conjg(psin(ix,iy,iz,is))*pswk1(ix,iy,iz,is)
               End Do
            End Do
         End Do
      End Do
! time-odd quantities: current j (this is always present traditionally)
      If(itimrev == 0) Then
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  currnt(ix, iy, iz, 1, iq) = currnt(ix, iy, iz, 1, iq) + weight * aimag(conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy,&
                                  & iz, 1)+conjg(psin(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
               End Do
            End Do
         End Do
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  Do is = 1, 2
                     currnts(ix,iy,iz,1,is,iq) = currnts(ix,iy,iz,1,is,iq) + &
                                                & weight*aimag(conjg(psin(ix, iy, iz, is))*pswk1(ix, iy, iz, is))
                  End Do
               End Do
            End Do
         End Do
! kinetic energy density (vector part): T
         If(itheta == 1 .And. iodds == 1) Then
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     kinvden(ix, iy, iz, 1, iq) = kinvden(ix, iy, iz, 1, iq) + 2.0_wp * weight * real(conjg(pswk1(ix, iy, iz, &
                    & 1))*pswk1(ix, iy, iz, 2))
                     kinvden(ix, iy, iz, 2, iq) = kinvden(ix, iy, iz, 2, iq) + 2.0_wp * weight * aimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*pswk1(ix, iy, iz, 2))
                     kinvden(ix, iy, iz, 3, iq) = kinvden(ix, iy, iz, 3, iq) + weight * (conjg(pswk1(ix, iy, iz, 1))*pswk1(ix, &
                    & iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                  End Do
               End Do
            End Do
         End If
      End If
! spin-orbit current tensor J_mu_nu (this we calculate for iodds=0 since we use it
! to calculate spin-orbit density. itheta dependence handled in skyrme.f90).
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               spincur(ix, iy, iz, 1, 1, iq) = spincur(ix, iy, iz, 1, 1, iq) - weight * aimag(conjg(pswk1(ix, iy, iz, &
              & 1))*psin(ix, iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
               spincur(ix, iy, iz, 1, 2, iq) = spincur(ix, iy, iz, 1, 2, iq) + weight * real(conjg(pswk1(ix, iy, iz, 1))*psin(ix,&
              & iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
               spincur(ix, iy, iz, 1, 3, iq) = spincur(ix, iy, iz, 1, 3, iq) - weight * aimag(conjg(pswk1(ix, iy, iz, &
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
                  tau(ix, iy, iz, iq) = tau(ix, iy, iz, iq) + weight*pswk1(ix, iy, iz, is)*conjg(pswk1(ix, iy, iz, is))
                  taus(ix,iy,iz,is,iq) = taus(ix,iy,iz,is,iq) + weight*pswk1(ix,iy,iz,is)*conjg(pswk1(ix,iy,iz,is))
                  drho(ix,iy,iz,2,is,iq) = drho(ix,iy,iz,2,is,iq) + weight*2.0_wp*conjg(psin(ix,iy,iz,is))*pswk1(ix,iy,iz,is)
               End Do
            End Do
         End Do
      End Do
! time-odd quantities: current j (this is always present traditionally)
      If(itimrev == 0) Then
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  currnt(ix, iy, iz, 2, iq) = currnt(ix, iy, iz, 2, iq) + weight * aimag(conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy,&
                 & iz, 1)+conjg(psin(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
               End Do
            End Do
         End Do
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  Do is = 1, 2
                     currnts(ix,iy,iz,2,is,iq) = currnts(ix,iy,iz,2,is,iq) + &
                                                 & weight*aimag(conjg(psin(ix,iy,iz,is))*pswk1(ix,iy,iz,is))
                  End Do
               End Do
            End Do
         End Do
! kinetic energy density T (vector-part)
         If(itheta == 1 .And. iodds == 1) Then
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     kinvden(ix, iy, iz, 1, iq) = kinvden(ix, iy, iz, 1, iq) + 2.0_wp * weight * real(conjg(pswk1(ix, iy, iz, &
                    & 1))*pswk1(ix, iy, iz, 2))
                     kinvden(ix, iy, iz, 2, iq) = kinvden(ix, iy, iz, 2, iq) + 2.0_wp * weight * aimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*pswk1(ix, iy, iz, 2))
                     kinvden(ix, iy, iz, 3, iq) = kinvden(ix, iy, iz, 3, iq) + weight * (conjg(pswk1(ix, iy, iz, 1))*pswk1(ix, &
                    & iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                  End Do
               End Do
            End Do
         End If
      End If
! spin-orbit current tensor J_mu_nu (this we calculate for iodds=0 since we use it
! to calculate spin-orbit density. itheta dependence handled in skyrme.f90).
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               spincur(ix, iy, iz, 2, 1, iq) = spincur(ix, iy, iz, 2, 1, iq) - weight * aimag(conjg(pswk1(ix, iy, iz, &
              & 1))*psin(ix, iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
               spincur(ix, iy, iz, 2, 2, iq) = spincur(ix, iy, iz, 2, 2, iq) + weight * real(conjg(pswk1(ix, iy, iz, 1))*psin(ix, &
              & iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
               spincur(ix, iy, iz, 2, 3, iq) = spincur(ix, iy, iz, 2, 3, iq) - weight * aimag(conjg(pswk1(ix, iy, iz, &
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
                  tau(ix, iy, iz, iq) = tau(ix, iy, iz, iq) + weight*pswk1(ix, iy, iz, is)*conjg(pswk1(ix, iy, iz, is))
                  taus(ix,iy,iz,is,iq) = taus(ix,iy,iz,is,iq) + weight*pswk1(ix,iy,iz,is)*conjg(pswk1(ix,iy,iz,is))
                  drho(ix,iy,iz,3,is,iq) = drho(ix,iy,iz,3,is,iq) + weight*2.0_wp*conjg(psin(ix,iy,iz,is))*pswk1(ix,iy,iz,is)
               End Do
            End Do
         End Do
      End Do
! time-odd quantities: current j (this is always present traditionally)
      If(itimrev == 0) Then
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  currnt(ix, iy, iz, 3, iq) = currnt(ix, iy, iz, 3, iq) + weight * aimag(conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, &
                 & iz, 1)+conjg(psin(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
               End Do
            End Do
         End Do
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  Do is = 1, 2
                     currnts(ix,iy,iz,3,is,iq) = currnts(ix,iy,iz,3,is,iq) + &
                                                 & weight*aimag(conjg(psin(ix,iy,iz,is))*pswk1(ix,iy,iz,is))
                  End Do
               End Do
            End Do
         End Do
! kinetic energy density (vector-part)
         If(itheta == 1 .And. iodds == 1) Then
            Do ix = 1, ncolx
               Do iy = 1, ncoly
                  Do iz = 1, ncolz
                     kinvden(ix, iy, iz, 1, iq) = kinvden(ix, iy, iz, 1, iq) + 2.0_wp * weight * real(conjg(pswk1(ix, iy, iz, &
                    & 1))*pswk1(ix, iy, iz, 2))
                     kinvden(ix, iy, iz, 2, iq) = kinvden(ix, iy, iz, 2, iq) + 2.0_wp * weight * aimag(conjg(pswk1(ix, iy, iz, &
                    & 1))*pswk1(ix, iy, iz, 2))
                     kinvden(ix, iy, iz, 3, iq) = kinvden(ix, iy, iz, 3, iq) + weight * (conjg(pswk1(ix, iy, iz, 1))*pswk1(ix, &
                    & iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*pswk1(ix, iy, iz, 2))
                  End Do
               End Do
            End Do
         End If
      End If
! spin-orbit current tensor J_mu_nu (this we calculate for iodds=0 since we use it
! to calculate spin-orbit density. itheta dependence handled in skyrme.f90).
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               spincur(ix, iy, iz, 3, 1, iq) = spincur(ix, iy, iz, 3, 1, iq) - weight * aimag(conjg(pswk1(ix, iy, iz, &
              & 1))*psin(ix, iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
               spincur(ix, iy, iz, 3, 2, iq) = spincur(ix, iy, iz, 3, 2, iq) + weight * real(conjg(pswk1(ix, iy, iz, 1))*psin(ix, &
              & iy, iz, 2)-conjg(psin(ix, iy, iz, 1))*pswk1(ix, iy, iz, 2))
               spincur(ix, iy, iz, 3, 3, iq) = spincur(ix, iy, iz, 3, 3, iq) - weight * aimag(conjg(pswk1(ix, iy, iz, &
              & 1))*psin(ix, iy, iz, 1)-conjg(pswk1(ix, iy, iz, 2))*psin(ix, iy, iz, 2))
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        calculate the old spin-orbit current vector J
!-----------------------------------------------------------------------
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               sodens(ix, iy, iz, 1, iq) = spincur(ix, iy, iz, 2, 3, iq) - spincur(ix, iy, iz, 3, 2, iq)
               sodens(ix, iy, iz, 2, iq) = spincur(ix, iy, iz, 3, 1, iq) - spincur(ix, iy, iz, 1, 3, iq)
               sodens(ix, iy, iz, 3, iq) = spincur(ix, iy, iz, 1, 2, iq) - spincur(ix, iy, iz, 2, 1, iq)
            End Do
         End Do
      End Do
!
      Return
End Subroutine densit
