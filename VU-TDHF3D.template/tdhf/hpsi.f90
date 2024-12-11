Subroutine hpsi(ncolx, ncoly, ncolz, upotq, bmassq, xkinpq, hsigmaq, cqq, dcqq, bmunuq, dbmuq, dbmasq, dxknpq, der1x, der1y, &
& der1z, der2x, der2y, der2z, eshift, pinn, pout, pswk, itimrev)
!-----------------------------------------------------------------------
!        hpsi:  program to form the product of h on the
!        state vector pinn and return the result in pout
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: itimrev
      Real(wp), Intent(In) :: eshift
      Real(wp), Intent(In) :: upotq(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: bmassq(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: xkinpq(ncolx, ncoly, ncolz, 3)
      Real(wp), Intent(In) :: hsigmaq(ncolx, ncoly, ncolz, 3)
      Real(wp), Intent(In) :: cqq(ncolx, ncoly, ncolz, 3)
      Real(wp), Intent(In) :: dcqq(ncolx, ncoly, ncolz, 3, 3)
      Real(wp), Intent(In) :: bmunuq(ncolx, ncoly, ncolz, 3, 3)
      Real(wp), Intent(In) :: dbmuq(ncolx, ncoly, ncolz, 3)
      Real(wp), Intent(In) :: dbmasq(ncolx, ncoly, ncolz, 3)
      Real(wp), Intent(In) :: dxknpq(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: der1x(ncolx, ncolx)
      Real(wp), Intent(In) :: der1y(ncoly, ncoly)
      Real(wp), Intent(In) :: der1z(ncolz, ncolz)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Complex(wp), Intent(Inout) :: pout(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pinn(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pswk(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Complex(wp), Parameter :: eye = (0.0_wp, 1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: is, iz, iy, ix, ic
      Real(wp) :: sigis, p, q, r
!-----------------------------------------------------------------------
!        diagonal part of h times psi (U_q part)
!-----------------------------------------------------------------------
      Do is = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = (upotq(ix, iy, iz)-eshift) * pinn(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! diagonal part of I_q term ( -i/2[Del.I_q]psi ) (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye / 2.0_wp * dxknpq(ix, iy, iz) * pinn(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
! this is the big sigma_q term (time-odd)
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     p = hsigmaq(ix, iy, iz, 1)
                     q = hsigmaq(ix, iy, iz, 2)
                     r = hsigmaq(ix, iy, iz, 3)
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) + (p*pinn(ix, iy, iz, ic)-eye*q*sigis*pinn(ix, iy, iz, &
                    & ic)+r*sigis*pinn(ix, iy, iz, is))
                  End Do
               End Do
            End Do
         End Do
      End If
! this is the -i/2[Del.B_q]psi (tensor) part of the h*psi
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  p = dbmuq(ix, iy, iz, 1)
                  q = dbmuq(ix, iy, iz, 2)
                  r = dbmuq(ix, iy, iz, 3)
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye / 2.0_wp * (p*pinn(ix, iy, iz, ic)-eye*q*sigis*pinn(ix, iy, &
                 & iz, ic)+r*sigis*pinn(ix, iy, iz, is))
               End Do
            End Do
         End Do
      End Do
!-----------------------------------------------------------------------
!        terms which require d/dx psi
!-----------------------------------------------------------------------
      Call cmulx(ncolx, ncoly, ncolz, der1x, pinn, pswk, 0)
! part of the effective mass term
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - dbmasq(ix, iy, iz, 1) * pswk(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! the other part of the I_q term ( -eye*I_q*Dpsi term ) (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye * xkinpq(ix, iy, iz, 1) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
! this is the -eye*B_q.Dsigma (tensor) part of h*psi
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  p = bmunuq(ix, iy, iz, 1, 1)
                  q = bmunuq(ix, iy, iz, 1, 2)
                  r = bmunuq(ix, iy, iz, 1, 3)
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye * (p*pswk(ix, iy, iz, ic)-eye*q*sigis*pswk(ix, iy, iz, &
                 & ic)+r*sigis*pswk(ix, iy, iz, is))
               End Do
            End Do
         End Do
      End Do
! this is the [Del.sigma.C_q]Del psi part of the C_q term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - (dcqq(ix, iy, iz, 1, 1)-sigis*eye*dcqq(ix, iy, iz, 1, 2)) * &
                    & pswk(ix, iy, iz, ic) - sigis * dcqq(ix, iy, iz, 1, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms which require d/dy psi
!-----------------------------------------------------------------------
      Call cmuly(ncolx, ncoly, ncolz, der1y, pinn, pswk, 0)
! part of the effective mass term
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - dbmasq(ix, iy, iz, 2) * pswk(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! part of -eye*I_q*Dpsi term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye * xkinpq(ix, iy, iz, 2) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
! this is the -eye*B_q.Dsigma (tensor) part of h*psi
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  p = bmunuq(ix, iy, iz, 2, 1)
                  q = bmunuq(ix, iy, iz, 2, 2)
                  r = bmunuq(ix, iy, iz, 2, 3)
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye * (p*pswk(ix, iy, iz, ic)-eye*q*sigis*pswk(ix, iy, iz, &
                 & ic)+r*sigis*pswk(ix, iy, iz, is))
               End Do
            End Do
         End Do
      End Do
! this is the [Del.sigma.C_q]Del psi part of the C_q term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - (dcqq(ix, iy, iz, 2, 1)-sigis*eye*dcqq(ix, iy, iz, 2, 2)) * &
                    & pswk(ix, iy, iz, ic) - sigis * dcqq(ix, iy, iz, 2, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms which require d/dz psi
!-----------------------------------------------------------------------
      Call cmulz(ncolx, ncoly, ncolz, der1z, pinn, pswk, 0)
! part of the effective mass term
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - dbmasq(ix, iy, iz, 3) * pswk(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! part of the -eye*I_q*Dpsi term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye * xkinpq(ix, iy, iz, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
! this is the -eye*B_q.Dsigma (tensor) part of h*psi
      Do is = 1, 2
         ic = 3 - is
         sigis = 3 - 2 * is
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  p = bmunuq(ix, iy, iz, 3, 1)
                  q = bmunuq(ix, iy, iz, 3, 2)
                  r = bmunuq(ix, iy, iz, 3, 3)
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - eye * (p*pswk(ix, iy, iz, ic)-eye*q*sigis*pswk(ix, iy, iz, &
                 & ic)+r*sigis*pswk(ix, iy, iz, is))
               End Do
            End Do
         End Do
      End Do
! this is the [Del.sigma.C_q]Del psi part of the C_q term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - (dcqq(ix, iy, iz, 3, 1)-sigis*eye*dcqq(ix, iy, iz, 3, 2)) * &
                    & pswk(ix, iy, iz, ic) - sigis * dcqq(ix, iy, iz, 3, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms which require (d/dx)**2 psi
!-----------------------------------------------------------------------
      Call cmulx(ncolx, ncoly, ncolz, der2x, pinn, pswk, 0)
! Effective mass, kinetic energy term
      Do is = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - bmassq(ix, iy, iz) * pswk(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! This is the second part of C_q term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - (cqq(ix, iy, iz, 1)-sigis*eye*cqq(ix, iy, iz, 2)) * pswk(ix, &
                    & iy, iz, ic) - sigis * cqq(ix, iy, iz, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms which require (d/dy)**2 psi
!-----------------------------------------------------------------------
      Call cmuly(ncolx, ncoly, ncolz, der2y, pinn, pswk, 0)
! Effective mass, kinetic energy term
      Do is = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - bmassq(ix, iy, iz) * pswk(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! This is the second part of C_q term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - (cqq(ix, iy, iz, 1)-sigis*eye*cqq(ix, iy, iz, 2)) * pswk(ix, &
                    & iy, iz, ic) - sigis * cqq(ix, iy, iz, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!-----------------------------------------------------------------------
!        terms which require (d/dz)**2 psi
!-----------------------------------------------------------------------
      Call cmulz(ncolx, ncoly, ncolz, der2z, pinn, pswk, 0)
! Effective mass, kinetic energy term
      Do is = 1, 2
         Do iz = 1, ncolz
            Do iy = 1, ncoly
               Do ix = 1, ncolx
                  pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - bmassq(ix, iy, iz) * pswk(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
! This is the second part of C_q term (time-odd)
      If(itimrev == 0) Then
         Do is = 1, 2
            ic = 3 - is
            sigis = 3 - 2 * is
            Do iz = 1, ncolz
               Do iy = 1, ncoly
                  Do ix = 1, ncolx
                     pout(ix, iy, iz, is) = pout(ix, iy, iz, is) - (cqq(ix, iy, iz, 1)-sigis*eye*cqq(ix, iy, iz, 2)) * pswk(ix, &
                    & iy, iz, ic) - sigis * cqq(ix, iy, iz, 3) * pswk(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!
      Return
End Subroutine hpsi
