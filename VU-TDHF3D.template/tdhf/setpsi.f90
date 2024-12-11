Subroutine setpsi(lpsi, ncolx, ncoly, ncolz, npsitot, npmin, nshell, radinx, radiny, radinz, wx, wy, wz, xclx, xcly, xclz, psi, &
& tpsi, xcm, ycm, zcm, atheta, btheta, gtheta, nneut, nprot, vocc, listis, itimrev, inuc, &
& nsplx, nsply, nsplz, nnotx, nnoty, nnotz, nordx, nordy, nordz, xkvx, xkvy, xkvz, iunit, isrestart)
!---------------------------------------------------------------------------
!     Initializes the wave functions with gaussian times polynomial for
!     a single nucleus. The width of the gaussians for each direction is
!     given by radinx(3). The polynomial is simply x**nx*y**ny*z**nz where
!     nx,ny,nz are Cartesian oscillator quanta.
!
!     In case of time-reversal invariance we initialize spinors with upper
!     component and zero lower component. If we Do not impose time-reversal
!     invariance Then the quantum degenerate states to the above set are
!     created with zero upper component and the same lower component. Due
!     to this choice the quantum degenerate states are orthogonal to each
!     other.
!
!     We also have the possibility of rotating the initial deformed nucleus
!     with respect to the unrotated solution (initial orientation determines
!     final orientation). The unrotated orientation is determined by the way
!     oscillator quanta are chosen. The rotation is done by evaluating the
!     initital quess on a rotated axes system determined by the Euler angles.
!     This is followed by a rotation in spin space (although it is may not be
!     neccessary here because of the way spinors are created here).
!
!     Quantization axis:
!     iquant = 1 (x-axis), =2 (y-axis), =3 (z-axis)
!     Default iquant=3 (for DC-TDHF iquant=1 is more suitable)
!     Subroutine sinfo chooses this axis from the principal axis.
!----------------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Logical, Intent(In) :: isrestart
      Integer :: ipers = 3
      Integer, Intent(In) :: lpsi
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: itimrev
      Integer, Intent(In) :: inuc
      Integer, Intent(In) :: npsitot(2)
      Integer, Intent(In) :: npmin(2)
      Integer, Intent(In) :: iunit
      Integer, Intent(In) :: nsplx
      Integer, Intent(In) :: nsply
      Integer, Intent(In) :: nsplz
      Integer, Intent(In) :: nnotx
      Integer, Intent(In) :: nnoty
      Integer, Intent(In) :: nnotz
      Integer, Intent(In) :: nordx
      Integer, Intent(In) :: nordy
      Integer, Intent(In) :: nordz
      Real(wp), Intent(In) :: xkvx(nnotx)
      Real(wp), Intent(In) :: xkvy(nnoty)
      Real(wp), Intent(In) :: xkvz(nnotz)
      Real(wp), Intent(In) :: radinx
      Real(wp), Intent(In) :: radiny
      Real(wp), Intent(In) :: radinz
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: xcm(2), ycm(2), zcm(2)
      Real(wp), Intent(In) :: atheta, btheta, gtheta
      Real(wp), Intent(Out) :: vocc(lpsi, 2)
      Real(wp)              :: bwrkx(nordx)
      Real(wp)              :: bwrky(nordy)
      Real(wp)              :: bwrkz(nordz)
      Real(wp)              :: tmpx(nordx,nordx+1)
      Real(wp)              :: tmpy(nordy,nordy+1)
      Real(wp)              :: tmpz(nordz,nordz+1)
      Integer, Intent(Out) :: nshell(3, lpsi, 2)
      Integer, Intent(Out) :: listis(lpsi, 2)
      Complex(wp), Intent(Inout) :: psi(ncolx, ncoly, ncolz, 2, lpsi, 2)
      Complex(wp), Intent(Inout) :: tpsi(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ix, iy, iz, iq, nneut, nprot, ncolxr, ncolyr, ncolzr, maxpart
      Integer :: nper, nps, nbode, nst, nbod, ka, shiftx, shifty, shiftz
      Integer :: i, j, k, kat, ismax, is, mmn, lm, l1, lx, il
      Integer :: jy, jx, jz, iiz, iiy, iix, nordxr, nordyr, nordzr
      Integer :: indexx, indexy, indexz, shiftxr, shiftyr, shiftzr
      Real(wp) :: anorm, x2, y2, z2, z, y, x, gshift(3), gshiftr(3), ethetar(3), xcmr(2), ycmr(2), zcmr(2)
      Real(wp) :: xp, yp, zp, xx, yy, zz, etheta1, etheta2, etheta3, voccr(lpsi,2)
      Real(wp), Dimension(:), Allocatable :: xclxr, xclyr, xclzr
      Complex(wp), Dimension(:,:,:,:,:,:), Allocatable :: bspl_coef
      Complex(wp), Dimension(:,:,:,:,:,:), Allocatable :: bspl_coefr
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Integer, Parameter  :: iquant = 3
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: small = 1.0e-25_wp
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
!-----------------------------------------------------------------------
!        calculate the shell occupation numbers and store the spins
!        of each single particle state
!
!                   nshell    =  nshell(direction,state,isospin)
!                                defines the node structure of the
!                                initial wave functions;
!                                for each state with given isospin
!                                one needs to specify the number of
!                                nodes in the directions
!                                1 = x, 2 = y, 3 = z .
!                                essentially the cartesian principle
!                                quantum numbers
!
!                   vocc      =  vocc(state,isospin)
!                                defines the s.p. occupations
!-----------------------------------------------------------------------
      nper = itimrev + 1
!
      iso: Do iq = 1, 2
         nps = npsitot(iq) - npmin(iq) + 1
         nbode = nneut
         If(iq == 2) nbode = nprot
!
         nst = 0
         nbod = 0
         Do ka = 0, nps
            Do k = 0, nps
               Do j = 0, nps
                  Do i = 0, nps
                     kat = i + j + k
                     If(kat /= ka) Cycle
                     ismax = 1
! if no time-reversal generate two sets of the same quantum numbers
                     If(itimrev == 0) ismax = 2
                     Do is = 1, ismax
                        nbod = nbod + nper
                        nst = nst + 1
                        If(nst > nps) Cycle iso
                        If(inuc == 1) Then
                           nshell(1, nst, iq) = i
                           nshell(2, nst, iq) = j
                           nshell(3, nst, iq) = k
                        Else
                           nshell(1, nst, iq) = i
                           nshell(2, nst, iq) = j
                           nshell(3, nst, iq) = k
                        End If
                        vocc(nst+npmin(iq)-1, iq) = nper
                        If(nbod > nbode) vocc(nst+npmin(iq)-1, iq) = 0
                        listis(nst, iq) = is
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do iso
!
      Do iq = 1, 2
         Write(8, '(/,A,I3,/)') ' ordering and nodes of the states with isospin ', iq
         mmn = npsitot(iq) - npmin(iq) + 1
         lm = 1
         Do While(mmn >  0)
            l1 = Min(mmn, 20)
            lx = lm + l1 - 1
            Write(8, '(/,A,20(I5))') ' state    = ', (il, il=lm, lx)
            Write(8, '(A,20(I5))') ' x-node     ', (nshell(1, nst, iq), nst=lm, lx)
            Write(8, '(A,20(I5))') ' y-node     ', (nshell(2, nst, iq), nst=lm, lx)
            Write(8, '(A,20(I5))') ' z-node     ', (nshell(3, nst, iq), nst=lm, lx)
            Write(8, '(A,20(F5.2))') ' occupation ', (vocc(nst+npmin(iq)-1, iq), nst=lm, lx)
            mmn = mmn - l1
            lm = lm + 20
         End Do
      End Do
      If (isrestart) Then
         Rewind(iunit)
         gshift(1) = (xclx(1)+xclx(ncolx))/2.0_wp
         gshift(2) = (xcly(1)+xcly(ncoly))/2.0_wp
         gshift(3) = (xclz(1)+xclz(ncolz))/2.0_wp
         maxpart = max(npsitot(1) - npmin(1) + 1,npsitot(2) - npmin(2) + 1)
         Read(iunit) nordxr, nordyr, nordzr
         If(nordxr /= nordx .OR. nordyr /= nordy .OR. nordzr /= nordz) Then
            Write(8, *) "Spline orders in restart of nucleus ",iunit-18," do not match! Terminating program."
            Stop
         End If
         Read(iunit) ncolxr, ncolyr, ncolzr
!
         Allocate(xclxr(ncolxr), xclyr(ncolyr), xclzr(ncolzr))
         Read(iunit) xclxr, xclyr, xclzr
         gshiftr(1) = (xclxr(1)+xclxr(ncolxr))/2.0_wp
         gshiftr(2) = (xclyr(1)+xclyr(ncolyr))/2.0_wp
         gshiftr(3) = (xclzr(1)+xclzr(ncolzr))/2.0_wp
!
         Allocate(bspl_coefr(ncolxr,ncolyr,ncolzr,2,maxpart,2))
         Allocate(bspl_coef(ncolx,ncoly,ncolz,2,lpsi,2))
         bspl_coef = (0.0_wp, 0.0_wp)
!
         Read(iunit) ethetar
         Read(iunit) xcmr, ycmr, zcmr
         Read(iunit) voccr(1:npsitot(1)-npmin(1)+1,1), voccr(1:npsitot(2)-npmin(2)+1,2)
         vocc(npmin(1):npsitot(1),1) = voccr(1:npsitot(1)-npmin(1)+1,1)
         vocc(npmin(2):npsitot(2),2) = voccr(1:npsitot(2)-npmin(2)+1,2)
         shiftxr = 0
         shiftyr = 0
         shiftzr = 0
         shiftx = 0
         shifty = 0
         shiftz = 0
         If(ncolxr > ncolx) shiftxr = (ncolxr - ncolx)/2
         If(ncolyr > ncoly) shiftyr = (ncolyr - ncoly)/2
         If(ncolzr > ncolz) shiftzr = (ncolzr - ncolz)/2
         If(ncolxr < ncolx) shiftx = (ncolx - ncolxr)/2
         If(ncolyr < ncoly) shifty = (ncoly - ncolyr)/2
         If(ncolzr < ncolz) shiftz = (ncolz - ncolzr)/2
         Do iq = 1, 2
            Read(iunit) bspl_coefr(:,:,:,:,1:npsitot(iq)-npmin(iq)+1,iq)
            bspl_coef(1+shiftx:ncolx-shiftx,1+shifty:ncoly-shifty,1+shiftz:ncolz-shiftz,:,npmin(iq):npsitot(iq),iq) =&
            & bspl_coefr(1+shiftxr:ncolxr-shiftxr,1+shiftyr:ncolyr-shiftyr,1+shiftzr:ncolzr-shiftzr,:,1:npsitot(iq)-npmin(iq)+1,iq)
            psi(:,:,:,:,npmin(iq):npsitot(iq),iq) = (small, small)
         End Do
         Do iq = 1, 2
            etheta1 = atheta/180.0_wp*pi !+ ethetar(1)
            etheta2 = btheta/180.0_wp*pi !- ethetar(2)
            etheta3 = gtheta/180.0_wp*pi
!$OMP PARALLEL DO DEFAULT(SHARED)&
!$OMP REDUCTION(+:psi)&
!$OMP PRIVATE(x, y, z, xp, yp, zp, iix, iiy, iiz, indexx, indexy, indexz, tmpx, tmpy, tmpz, bwrkx, bwrky, bwrkz)
             Do iz = 1, ncolz
                Do iy = 1, ncoly
                   Do ix = 1, ncolx
                      x = xclx(ix) - xcm(iq)
                      y = xcly(iy) - ycm(iq)
                      z = xclz(iz) - zcm(iq)
                      xp =  x*(cos(etheta3)*cos(etheta2)*cos(etheta1) - sin(etheta3)*sin(etheta1)) +&
                            y*(cos(etheta3)*cos(etheta2)*sin(etheta1) + sin(etheta3)*cos(etheta1)) -&
                            z*cos(etheta3)*sin(etheta2) + xcmr(iq) + (gshift(1)-gshiftr(1))
                      yp =  x*(-sin(etheta3)*cos(etheta2)*cos(etheta1) - cos(etheta3)*sin(etheta1)) +&
                            y*(-sin(etheta3)*cos(etheta2)*sin(etheta1) + cos(etheta3)*cos(etheta1)) +&
                            z*sin(etheta3)*sin(etheta2) + ycmr(iq) + (gshift(2)-gshiftr(2))
                      zp =  x*sin(etheta2)*cos(etheta1) +&
                            y*sin(etheta2)*sin(etheta1) +&
                            z*cos(etheta2) + zcmr(iq) + (gshift(3)-gshiftr(3))
                      If(xp > xclx(ncolx)) Then
                         If (Abs(xp - xclx(ncolx)) < 1.0e-8_wp ) Then
                            xp = xclx(ncolx)
                         Else
                           Cycle
                         End If
                      End If
                      If(yp > xcly(ncoly)) Then
                         If (Abs(yp - xcly(ncoly)) < 1.0e-8_wp ) Then
                             yp = xcly(ncoly)
                           Else
                             Cycle
                         End If
                      End If
                      If(zp > xclz(ncolz)) Then
                         If (Abs(zp - xclz(ncolz)) < 1.0e-8_wp ) Then
                             zp = xclz(ncolz)
                         Else
                            Cycle
                         End If
                      End If
                      If(xp < xclx(1)) Then
                         If (Abs(xp - xclx(1)) < 1.0e-8_wp ) Then
                             xp = xclx(1)
                         Else
                             Cycle
                         End If
                      End If
                      If(yp < xcly(1)) Then
                         If (Abs(yp - xcly(1)) < 1.0e-8_wp ) Then
                             yp = xcly(1)
                         Else
                             Cycle
                         End If
                      End If
                      If(zp < xclz(1)) Then
                         If (Abs(zp - xclz(1)) < 1.0e-8_wp ) Then
                             zp = xclz(1)
                         Else
                             Cycle
                         End If
                      End If
                      Call splgens(0, xkvz, zp, bwrkz, tmpz, indexz, nordz, nordz + 1, nnotz, nsplz)
                      Call splgens(0, xkvy, yp, bwrky, tmpy, indexy, nordy, nordy + 1, nnoty, nsply)
                      Call splgens(0, xkvx, xp, bwrkx, tmpx, indexx, nordx, nordx + 1, nnotx, nsplx)
                      Do nst = npmin(iq), npsitot(iq)
                         Do is = 1, 2
                            Do jz = 1, nordz
                               iiz = indexz - jz + 1
                               If (iiz > ncolz .and. ipers == 3) iiz = iiz - ncolz
                               Do jy = 1, nordy
                                  iiy = indexy - jy + 1
                                  If (iiy > ncoly .and. ipers == 3) iiy = iiy - ncoly
                                  Do jx = 1, nordx
                                     iix = indexx - jx + 1
                                     If (iix > ncolx .and. ipers == 3) iix = iix - ncolx
                                     psi(ix,iy,iz,is,nst,iq) = psi(ix,iy,iz,is,nst,iq)&
                                                             + bspl_coef(iix,iiy,iiz,is,nst,iq)*bwrkx(jx)&
                                                             * bwrky(jy)*bwrkz(jz)
                                  End Do
                               End Do
                            End Do
                         End Do
                      End Do
                   End Do
                End Do
             End Do
!$OMP END PARALLEL DO
         End Do
! Rotation in spin space
         If(itimrev == 0) Then
            Do iq = 1, 2
               Do nst = npmin(iq), npsitot(iq)
                  Do iz = 1, ncolz
                     Do iy = 1, ncoly
                        Do ix = 1, ncolx
                           tpsi(ix,iy,iz,1) = psi(ix,iy,iz,1,nst,iq)
                           psi(ix,iy,iz,1,nst,iq) = &
                           exp(-eye*etheta1/2.0_wp)*cos(etheta2/2.0_wp)*exp(-eye*etheta3/2.0_wp)*psi(ix,iy,iz,1,nst,iq) -&
                           exp(-eye*etheta1/2.0_wp)*sin(etheta2/2.0_wp)*exp(+eye*etheta3/2.0_wp)*psi(ix,iy,iz,2,nst,iq)
                           psi(ix,iy,iz,2,nst,iq) = &
                           exp(+eye*etheta1/2.0_wp)*sin(etheta2/2.0_wp)*exp(-eye*etheta3/2.0_wp)*tpsi(ix,iy,iz,1) +&
                           exp(+eye*etheta1/2.0_wp)*cos(etheta2/2.0_wp)*exp(+eye*etheta3/2.0_wp)*psi(ix,iy,iz,2,nst,iq)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End If
! If no restart for this nucleus, setup the psi array in the standard way
         Deallocate(bspl_coef, bspl_coefr)
      Else
!-----------------------------------------------------------------------
!        The lowest state is a plain gaussian. Note that unequal radinx,
!        radiny, and radinz will produce a deformed initial guess.
!        etheta? are the Euler angles of rotation.
!-----------------------------------------------------------------------
         anorm = 0.0_wp
         etheta1 = atheta / 180.0_wp * pi
         etheta2 = btheta / 180.0_wp * pi
         etheta3 = gtheta / 180.0_wp * pi
         Do ix = 1, ncolx
            Do iy = 1, ncoly
               Do iz = 1, ncolz
                  x = xclx(ix) - xcm(1)
                  y = xcly(iy) - ycm(1)
                  z = xclz(iz) - zcm(1)
                  xp = x * (Cos(etheta3)*Cos(etheta2)*Cos(etheta1)-Sin(etheta3)*Sin(etheta1)) + y * &
                 & (Cos(etheta3)*Cos(etheta2)*Sin(etheta1)+Sin(etheta3)*Cos(etheta1)) - z * Cos(etheta3) * Sin(etheta2)
                  yp = x * (-Sin(etheta3)*Cos(etheta2)*Cos(etheta1)-Cos(etheta3)*Sin(etheta1)) + y * &
                 & (-Sin(etheta3)*Cos(etheta2)*Sin(etheta1)+Cos(etheta3)*Cos(etheta1)) + z * Sin(etheta3) * Sin(etheta2)
                  zp = x * Sin(etheta2) * Cos(etheta1) + y * Sin(etheta2) * Sin(etheta1) + z * Cos(etheta2)
                  x2 = (xp/radinx) ** 2
                  y2 = (yp/radiny) ** 2
                  z2 = (zp/radinz) ** 2
                  tpsi(ix, iy, iz, 1) = Exp(-x2) * Exp(-y2) * Exp(-z2)
                  anorm = anorm + wx(ix) * wy(iy) * wz(iz) * conjg(tpsi(ix, iy, iz, 1)) * tpsi(ix, iy, iz, 1)
               End Do
            End Do
         End Do
!
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  tpsi(ix, iy, iz, 1) = tpsi(ix, iy, iz, 1) / (Sqrt(anorm)+small)
               End Do
            End Do
         End Do
!-----------------------------------------------------------------------
!        now the rest of the states as polynomials times the
!        gaussian. the order of the polynomial is transferred
!        via the shell-ordering fields nshell.
!-----------------------------------------------------------------------
         Do iq = 1, 2
            Do nst = 1, (npsitot(iq)-npmin(iq)+1)
               nps = nst + npmin(iq) - 1
! with time-reversal we generate spinors with non-zero upper component
               is = 1
! if no time-reversal generate the degenerate spinor with non-zero lower component
               If(itimrev == 0) is = listis(nst, iq)
!
               Do iz = 1, ncolz
                  Do iy = 1, ncoly
                     Do ix = 1, ncolx
                        x = xclx(ix) - xcm(1)
                        y = xcly(iy) - ycm(1)
                        z = xclz(iz) - zcm(1)
                        xp = x * (Cos(etheta3)*Cos(etheta2)*Cos(etheta1)-Sin(etheta3)*Sin(etheta1)) + y * &
                       & (Cos(etheta3)*Cos(etheta2)*Sin(etheta1)+Sin(etheta3)*Cos(etheta1)) - z * Cos(etheta3) * Sin(etheta2)
                        yp = x * (-Sin(etheta3)*Cos(etheta2)*Cos(etheta1)-Cos(etheta3)*Sin(etheta1)) + y * &
                       & (-Sin(etheta3)*Cos(etheta2)*Sin(etheta1)+Cos(etheta3)*Cos(etheta1)) + z * Sin(etheta3) * Sin(etheta2)
                        zp = x * Sin(etheta2) * Cos(etheta1) + y * Sin(etheta2) * Sin(etheta1) + z * Cos(etheta2)
                        xx = xp ** nshell(1, nst, iq)
                        yy = yp ** nshell(2, nst, iq)
                        zz = zp ** nshell(3, nst, iq)
                        psi(ix, iy, iz, is, nps, iq) = tpsi(ix, iy, iz, 1) * xx * yy * zz
! these lines are ok but mess up detailed comparison with other codes
                        If(iquant == 1 .and. is == 1) Then
                           psi(ix,iy,iz,1,nps,iq) = tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                           psi(ix,iy,iz,2,nps,iq) = tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                        ElseIf(iquant == 1 .and. is == 2) Then
                           psi(ix,iy,iz,1,nps,iq) =  tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                           psi(ix,iy,iz,2,nps,iq) = -tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                        End If
                        If(iquant == 2 .and. is == 1) Then
                           psi(ix,iy,iz,1,nps,iq) = tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                           psi(ix,iy,iz,2,nps,iq) = eye*tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                        ElseIf(iquant == 2 .and. is == 2) Then
                           psi(ix,iy,iz,1,nps,iq) =  tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                           psi(ix,iy,iz,2,nps,iq) = -eye*tpsi(ix,iy,iz,1)*xx*yy*zz/sqrt(2.0_wp)
                        End If
                        If(iquant == 3 .and. is == 1) Then
                           psi(ix,iy,iz,1,nps,iq) = tpsi(ix,iy,iz,1)*xx*yy*zz
                           psi(ix,iy,iz,2,nps,iq) = (0.0_wp,0.0_wp)
                        ElseIf(iquant == 3 .and. is == 2) Then
                           psi(ix,iy,iz,1,nps,iq) = (0.0_wp,0.0_wp)
                           psi(ix,iy,iz,2,nps,iq) = tpsi(ix,iy,iz,1)*xx*yy*zz
                        End If
                    End Do
                  End Do
               End Do
!
            End Do
         End Do
!
!       Do iq = 1, 2
!          Do nst = 1, (npsitot(iq) - npmin(iq) + 1)
!             nps = nst + npmin(iq) - 1
!
!             Do iz = 1, ncolz
!                Do iy = 1, ncoly
!                   Do ix = 1, ncolx
!                      tpsi(ix,iy,iz,1) = psi(ix,iy,iz,1,nps,iq)
!                      psi(ix,iy,iz,1,nps,iq) = &
!                      exp(-eye*etheta1/2.0_wp)*cos(etheta2/2.0_wp)*exp(-eye*etheta3/2.0_wp)*psi(ix,iy,iz,1,nps,iq) -&
!                      exp(-eye*etheta1/2.0_wp)*sin(etheta2/2.0_wp)*exp(+eye*etheta3/2.0_wp)*psi(ix,iy,iz,2,nps,iq)
!                      psi(ix,iy,iz,2,nps,iq) = &
!                      exp(+eye*etheta1/2.0_wp)*sin(etheta2/2.0_wp)*exp(-eye*etheta3/2.0_wp)*tpsi(ix,iy,iz,1) +&
!                      exp(+eye*etheta1/2.0_wp)*cos(etheta2/2.0_wp)*exp(+eye*etheta3/2.0_wp)*psi(ix,iy,iz,2,nps,iq)
!                   End Do
!                End Do
!             End Do
!
!          End Do
!       End Do
!
      End If
!
      Return
End Subroutine setpsi
