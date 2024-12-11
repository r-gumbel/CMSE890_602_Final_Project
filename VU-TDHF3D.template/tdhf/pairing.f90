Subroutine pairing(ipairn, ipairp, iter, irest, nx, ny, nz, lpsi, npmin, npsi, nneut, nprot, epairn, epairp,&
                   gapn, gapp, fermin, fermip, spenrg, vocc, v2, wxyz, rho, v0prot, v0neut, rho0pr, psi,&
                   ireord, itimrev)
!-----------------------------------------------------------------------
!      Driving subroutine to calculate pairing
!      Adopted from Sky3D pairing routine but time-reversal is imposed
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In)     :: nx, ny, nz
      Integer, Intent(In)     :: iter
      Integer, Intent(In)     :: irest
      Integer, Intent(In)     :: ipairn
      Integer, Intent(In)     :: ipairp
      Integer, Intent(In)     :: lpsi
      Integer, Intent(In)     :: nneut
      Integer, Intent(In)     :: nprot
      Integer, Intent(In)     :: npmin(2)
      Integer, Intent(In)     :: npsi(2)
      Integer, Intent(In)     :: itimrev
      Integer, Intent(In)     :: ireord
      Real(wp)                :: epairn(2)
      Real(wp)                :: epairp(2)
      Real(wp)                :: gapn
      Real(wp)                :: gapp
      Real(wp)                :: fermin(2)
      Real(wp)                :: fermip(2)
      Real(wp), Intent(In)    :: spenrg(lpsi,2)
      Real(wp), Intent(Out)   :: vocc(lpsi,2)
      Real(wp), Intent(Out)   :: v2(lpsi,2)
      Real(wp), Intent(In)    :: v0prot, v0neut, rho0pr
      Real(wp), Intent(In)    :: rho(nx, ny, nz, 2)
      Real(wp), Intent(In)    :: wxyz(nx, ny, nz)
      Complex(wp), Intent(In) :: psi(nx, ny, nz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp)                :: deltaf(lpsi, 2)
      Real(wp)                :: fk(lpsi,2)
      Real(wp)                :: eferm(2), epair(2), avdelt(2), avg(2)
      Real(wp)                :: particle_number, mass_number

      Integer, Dimension(:), Allocatable :: iorder
      Integer                 :: iq, ipair, n1, n2
      Save                       eferm
!
      If(irest /= 0) Then
         eferm(1) = fermin(1)
         eferm(2) = fermip(1)
      End If
      If(itimrev == 1) Then
         mass_number = (nneut + nprot)/2
      Else
         mass_number = (nneut + nprot)
      End If
      ipair = max(ipairn, ipairp)
! calculate gaps
      If(ipair /= 0) Then
         Call pairgap(iter, irest, nx, ny, nz, spenrg, deltaf, eferm, mass_number, npmin, npsi, v2, lpsi, psi,&
                     &ipair, itimrev, wxyz, rho, v0prot, v0neut, rho0pr, fk)
      End If
! solve pairing problem for each isospin iq
      Do iq = 1, 2
         If(iq == 1 .and. ipairn /= 0) Then
            If(itimrev == 1) Then
               particle_number = nneut/2
            Else
               particle_number = nneut
            End If
            Call pairdn(particle_number, v2, eferm, spenrg, deltaf, avdelt, epair, avg, npmin, npsi, lpsi, iq)
         Else If(iq == 2 .and. ipairp /= 0) Then
            If(itimrev == 1) Then
               particle_number = nprot/2
            Else
               particle_number = nprot
            End If
            Call pairdn(particle_number, v2, eferm, spenrg, deltaf, avdelt, epair, avg, npmin, npsi, lpsi, iq)
         End If
      End Do
! T-R
      If(itimrev == 1) Then
         If(ipairn /= 0) Then
            gapn = avg(1)/2.0_wp
            epairn(2) = 2*epair(1)
            fermin(1) = eferm(1)
            fermin(2) = avdelt(1)
            vocc(npmin(1):npsi(1),1) = 2*v2(npmin(1):npsi(1),1)
         Else
            gapn = 0.0_wp
            epairn(2) = 0.0_wp
            fermin(1) = spenrg(npsi(1),1)
            n1 = npsi(1) - npmin(1) + 1
            n2 = nneut/2
            If(ireord == 1) Then
               Allocate(iorder(n1))
               Call eorder(spenrg(npmin(1):npsi(1),1), iorder, n1)
               Call reorder(v2(npmin(1):npsi(1),1), iorder, n1, n2)
               vocc(npmin(1):npsi(1),1) = 2*v2(npmin(1):npsi(1),1)
               Deallocate(iorder)
            End If
         End If
         If(ipairp /= 0) Then
            gapp = avg(2)/2.0_wp
            epairp(2) = 2*epair(2)
            fermip(1) = eferm(2)
            fermip(2) = avdelt(2)
            vocc(npmin(2):npsi(2),2) = 2*v2(npmin(2):npsi(2),2)
         Else
            gapp = 0.0_wp
            epairp(2) = 0.0_wp
            fermip(1) = spenrg(npsi(2),2)
            n1 = npsi(2) - npmin(2) + 1
            n2 = nprot/2
            If(ireord == 1) Then
               Allocate(iorder(n1))
               Call eorder(spenrg(npmin(2):npsi(2),2), iorder, n1)
               Call reorder(v2(npmin(2):npsi(2),2), iorder, n1, n2)
               vocc(npmin(2):npsi(2),2) = 2*v2(npmin(2):npsi(2),2)
               Deallocate(iorder)
            End If
          End If
      Else
         If(ipairn /= 0) Then
            gapn = avg(1)
            epairn(2) = epair(1)
            fermin(1) = eferm(1)
            fermin(2) = avdelt(1)
            vocc(npmin(1):npsi(1),1) = v2(npmin(1):npsi(1),1)
         Else
            gapn = 0.0_wp
            epairn(2) = 0.0_wp
            fermin(1) = spenrg(npsi(1),1)
            n1 = npsi(1) - npmin(1) + 1
            n2 = nneut
            If(ireord == 1) Then
               Allocate(iorder(n1))
               Call eorder(spenrg(npmin(1):npsi(1),1), iorder, n1)
               Call reorder(v2(npmin(1):npsi(1),1), iorder, n1, n2)
               vocc(npmin(1):npsi(1),1) = v2(npmin(1):npsi(1),1)
               Deallocate(iorder)
            End If
         End If
         If(ipairp /= 0) Then
            gapp = avg(2)
            epairp(2) = epair(2)
            fermip(1) = eferm(2)
            fermip(2) = avdelt(2)
            vocc(npmin(2):npsi(2),2) = v2(npmin(2):npsi(2),2)
         Else
            gapp = 0.0_wp
            epairp(2) = 0.0_wp
            fermip(1) = spenrg(npsi(2),2)
            n1 = npsi(2) - npmin(2) + 1
            n2 = nprot
            If(ireord == 1) Then
               Allocate(iorder(n1))
               Call eorder(spenrg(npmin(2):npsi(2),2), iorder, n1)
               Call reorder(v2(npmin(2):npsi(2),2), iorder, n1, n2)
               vocc(npmin(2):npsi(2),2) = v2(npmin(2):npsi(2),2)
               Deallocate(iorder)
            End If
         End If
      End If
!      DO iq = 1, 2
!         WRITE(*,'(a,i2,a,4(1pg12.4))') 'iq=',iq,': ', eferm(iq), epair(iq), avdelt(iq), avg(iq)
!      ENDDO
End Subroutine pairing
Subroutine pairgap(iter, irest, nx, ny, nz, spenrg, deltaf, eferm, mass_number, npmin, npsi, vocc, lpsi, psi,&
                  &ipair, itimrev, wxyz, rho, v0prot, v0neut, rho0pr, fk)
      Implicit None
      Integer, Parameter   :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In)     :: nx, ny, nz
      Integer, Intent(In)     :: lpsi
      Integer, Intent(In)     :: iter
      Integer, Intent(In)     :: irest
      Integer, Intent(In)     :: ipair
      Integer, Intent(In)     :: itimrev
      Integer, Intent(In)     :: npmin(2), npsi(2)
      Real(wp), Intent(In)    :: v0prot, v0neut, rho0pr
      Real(wp), Intent(In)    :: rho(nx, ny, nz, 2)
      Real(wp), Intent(In)    :: wxyz(nx, ny, nz)
      Real(wp), Intent(In)    :: mass_number
      Real(wp), Intent(In)    :: spenrg(lpsi,2)
      Real(wp), Intent(Out)   :: deltaf(lpsi, 2)
      Real(wp), Intent(In)    :: eferm(2)
      Real(wp), Intent(In)    :: vocc(lpsi, 2)
      Real(wp), Intent(InOut) :: fk(lpsi,2)
      Complex(wp), Intent(In) :: psi(nx, ny, nz, 2, lpsi, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp), Parameter :: xsmall = 1.e-20_wp
      Real(wp), Parameter :: smallp = 0.000001_wp
      Integer, Parameter  :: itrsin = 10
      Integer             :: iq, nst, is
      Real(wp)            :: v0act, fact
      Real(wp)            :: paircrf(nx, ny, nz,2)
      Real(wp)            :: pairptf(nx, ny, nz,2)
      Real(wp)            :: deq, uq, xcut
!
      xcut = 0.0_wp
      If(ipair == 4) xcut = 1.0_wp
      deq = 5.0_wp
      uq = deq/10.0_wp
! constant gap in early stages of iteration
      If(iter <= itrsin .And. irest == 0) Then
         deltaf = 11.2_wp / Sqrt(mass_number)
         Return
      End If
      Do iq = 1, 2
         Do nst = npmin(iq), npsi(iq)
            fk(nst,iq) = (1.0_wp/(1.0_wp + exp((spenrg(nst,iq)-eferm(iq)-deq)/uq)))
            fk(nst,iq) = fk(nst,iq)*(1.0_wp-xcut+xcut/(1.0_wp + exp(-(spenrg(nst,iq)-eferm(iq)+deq)/uq)))
            If(spenrg(nst,iq) > 0.0_wp) fk(nst,iq) = fk(nst,iq)*1.0e-12_wp
         End Do
         If(ipair /= 3 .Or. ipair /= 4) fk = 1.0_wp
         If(itimrev == 1) then
            fact = 1.0_wp
         Else
            fact = 0.5_wp
         End If
! now the detailed gaps
         paircrf(:,:,:,iq) = 0.0_wp
! accumulate new pair-density
         Do nst = npmin(iq), npsi(iq)
            Do is = 1, 2
               paircrf(:,:,:,iq) = paircrf(:,:,:,iq)+fk(nst,iq)*Sqrt(Max(vocc(nst, iq)-vocc(nst, iq)**2, smallp))*fact* &
               (REAL(psi(:, :, :, is, nst, iq))**2 + AIMAG(psi(:, :, :, is, nst, iq))**2)
            End Do
         End Do
! determine pairing strength
         If(iq == 1) Then
            v0act = v0neut
         Else
            v0act = v0prot
         End If
! now multiply with strength to obtain local pair-potential
         If(ipair == 2 .Or. ipair == 3 .Or. ipair == 4) Then
            pairptf(:,:,:,iq) = v0act*(1.0_wp - (rho(:, :, :, 1) + rho(:, :, :, 2))/rho0pr)*paircrf(:,:,:,iq)
         Else
            pairptf(:,:,:,iq) = v0act*paircrf(:,:,:,iq)
         End If
!        write(*,*) "pairing energy", iq, SUM(wxyz*paircrf(:,:,:,iq)*pairptf(:,:,:,iq))
! finally compute the actual gaps as s.p. expectation values with the pair potential
         Do nst = npmin(iq), npsi(iq)
            deltaf(nst, iq) = 0.0_wp
            Do is = 1, 2
               deltaf(nst,iq) = deltaf(nst,iq) + SUM(wxyz*pairptf(:,:,:,iq)*conjg(psi(:,:,:,is,nst,iq))*psi(:,:,:,is,nst,iq))
            End Do
         End Do
      End Do
End Subroutine pairgap
Subroutine pairdn(particle_number, vocc, eferm, sp_energy, deltaf, avdelt, epair, avg, npmin, npsi, lpsi, iq)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Real(wp), Intent(In)  :: particle_number
      Real(wp), Intent(In)  :: sp_energy(lpsi,2), deltaf(lpsi,2)
      Real(wp), Intent(Out) :: vocc(lpsi,2), avdelt(2), epair(2), avg(2), eferm(2)
      Integer, Intent(In)   :: lpsi
      Integer, Intent(In)   :: iq
      Integer, Intent(In)   :: npmin(2), npsi(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp), Parameter   :: xsmall = 1.0e-20_wp
      Integer               :: it, nst
      Real(wp)              :: sumuv, sumduv, edif, delta2, equasi, v2, vol, res
!
! start with non-pairing value of Fermi energy
      it = npmin(iq) + Nint(particle_number) - 1
      eferm(iq) = 0.5_wp * (sp_energy(it, iq) + sp_energy(it+1, iq))
      vocc(npmin(iq):it, iq) = 1.0_wp
      vocc(it+1:npsi(iq), iq) = 0.0_wp
! determine true fermi energy for given 'delta' and 'e'
      call rbrent(particle_number, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq, res)
      eferm(iq) = res
! finally compute average gap and pairing energy
      sumuv = 0.0_wp
      sumduv = 0.0_wp
      Do nst = npmin(iq), npsi(iq)
         delta2 = deltaf(nst, iq)*deltaf(nst, iq)
         edif = sp_energy(nst, iq) - eferm(iq)
         equasi = Sqrt(edif*edif + delta2)
         v2 = 0.5_wp - 0.5_wp*edif / equasi
         vol = 0.5_wp * Sqrt(Max(v2-v2*v2, xsmall))
         sumuv = vol + sumuv
         sumduv = vol*deltaf(nst, iq) + sumduv
      End Do
      sumuv = Max(sumuv, xsmall)
      avdelt(iq) = sumduv/sumuv
      epair(iq) = sumduv
      avg(iq) = epair(iq)/sumuv**2
End Subroutine pairdn
Subroutine rbrent(particle_number, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq, res)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Real(wp), Intent(In) :: particle_number
      Real(wp), Intent(In) :: sp_energy(lpsi, 2), deltaf(lpsi, 2)
      Real(wp), Intent(Out) :: vocc(lpsi, 2)
      Integer, Intent(In)  :: npmin(2), npsi(2), lpsi, iq
      Real(wp), Intent(Out) :: res
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      data for wijngaarden-dekker-brent method for root
      Real(wp), Parameter :: ra0 = - 100.0_wp, rb0 = 100.0_wp, reps = 1.e-14_wp, rtol = 1.e-16_wp
      Integer, Parameter  :: iitmax = 100
      Real(wp)            :: ra, rb, rfa, rfb, rfc, rc, rd, re, rtol1, rxm, rp, rq, rr, rs, rmin1, rmin2
      Integer             :: ii
!
      ra = ra0
      rb = rb0
      Call bcs_occupation(ra, rfa, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq)
      rfa = particle_number - rfa
      Call bcs_occupation(rb, rfb, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq)
      rfb = particle_number - rfb
      rfc = rfb
      rc = rb
      rd = rb - ra
      re = rd
      If((rfa > 0.0_wp) .And. (rfb > 0.0_wp)) Then
         If(rfa > rfb) Then
            rb = 1.0_wp
         Else
            rb = 0.0_wp
         End If
      Else If((rfa < 0.0_wp) .And. (rfb < 0.0_wp)) Then
         If(rfa > rfb) Then
            rb = 0.0_wp
         Else
            rb = 1.0_wp
         End If
      Else
         iteration: Do ii = 1, iitmax
            If(((rfb > 0.0_wp) .And. (rfc > 0.0_wp)) .Or. ((rfb < 0.0_wp) .And. (rfc < 0.0_wp))) Then
               rc = ra
               rfc = rfa
               rd = rb - ra
               re = rd
            End If
            If(Abs(rfc) < Abs(rfb)) Then
               ra = rb
               rb = rc
               rc = ra
               rfa = rfb
               rfb = rfc
               rfc = rfa
            End If
!  convergence check
            rtol1 = 2.0_wp * reps * Abs(rb) + 0.5_wp * rtol
            rxm = 0.5_wp * (rc-rb)
            If((Abs(rxm) <= rtol1) .Or. ((Abs(rfb) == 0.0_wp))) Then
               res = rb
               Call bcs_occupation(res, rfa, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq)
               Return
            End If
            If((Abs(re) >= rtol1) .Or. (Abs(rfa) > Abs(rfb))) Then
               rs = rfb / rfa
               If(ra == rc) Then
                  rp = 2.0_wp * rxm * rs
                  rq = 1.0_wp - rs
               Else
                  rq = rfa / rfc
                  rr = rfb / rfc
                  rp = rs * (2.0_wp*rxm*rq*(rq-rr)-(rb-ra)*(rr-1.0_wp))
                  rq = (rq-1.0_wp) * (rr-1.0_wp) * (rs-1.0_wp)
               End If
               If(rp > 0.0_wp) Then
                  rq = - rq
               End If
               rp = Abs(rp)
               rmin1 = 3.0_wp * rxm * rq - Abs(rtol1*rq)
               rmin2 = Abs(rq*re)
               If(2.0_wp*rp < Min(rmin1, rmin2)) Then
                  re = rd
                  rd = rp / rq
               Else
                  rd = rxm
                  re = rd
               End If
            Else
               rd = rxm
               re = rd
            End If
            ra = rb
            rfa = rfb
            If(Abs(rd) > rtol1) Then
               rb = rb + rd
            Else
               rb = rb + SIGN(rtol1, rxm)
            End If
            Call bcs_occupation(rb, rfb, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq)
            rfb = particle_number - rfb
         End Do iteration
         Stop 'No solution found in pairing iterations'
      End If
End Subroutine rbrent
Subroutine bcs_occupation(efermi, bcs_partnum, npmin, npsi, sp_energy, vocc, deltaf, lpsi, iq)
      Implicit None
      Integer, Parameter    :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Real(wp), Parameter   :: smal = 1.0e-10_wp
      Real(wp), Intent(In)  :: efermi
      Real(wp), Intent(In)  :: sp_energy(lpsi,2), deltaf(lpsi,2)
      Real(wp), Intent(Out) :: bcs_partnum, vocc(lpsi,2)
      Integer, Intent(In)   :: npmin(2), npsi(2), lpsi, iq
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer               :: nst
      Real(wp)              :: delta2, equasi, edif
!
      bcs_partnum = 0.0_wp
      Do nst = npmin(iq), npsi(iq)
         delta2 = deltaf(nst, iq)*deltaf(nst, iq)
         edif = sp_energy(nst, iq) - efermi
         equasi = Sqrt(edif*edif+delta2)
         vocc(nst, iq) = 0.5_wp * (1.0_wp - edif/equasi)
         vocc(nst, iq) = Min(Max(vocc(nst, iq), smal), 1.0_wp - smal)
         bcs_partnum = bcs_partnum + vocc(nst, iq)
      End Do
End Subroutine bcs_occupation
