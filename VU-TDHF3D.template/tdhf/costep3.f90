Subroutine costep3(ncolx, ncoly, ncolz, xclx, xcly, xclz, xlam, damp, psin, psout, idp)
!-----------------------------------------------------------------------
!     extra constraint steps
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: idp
      Real(wp), Intent(In) :: xlam(7)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: damp(ncolx, ncoly, ncolz)
      Complex(wp), Intent(In) :: psin(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(Inout) :: psout(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: x, y, z, xqx
      Integer :: is, iz, iy, ix
!-----------------------------------------------
!
      If(idp == 1) Then !x-axis symmetry axis
         Do is = 1, 2
            Do iz = 1, ncolz
               z = xclz(iz)
               Do iy = 1, ncoly
                  y = xcly(iy)
                  Do ix = 1, ncolx
                     x = xclx(ix)
                     xqx = Sqrt(3.0_wp/(4.0_wp*pi))*x*xlam(1)! q10
                     xqx = xqx - Sqrt(3.0_wp/(8.0_wp*pi))*y*xlam(2)! Re(q11)
                     xqx = xqx + Sqrt(3.0_wp/(8.0_wp*pi))*z*xlam(3)! Im(q11)
                     xqx = xqx - Sqrt(15.0_wp/(8.0_wp*pi))*x*y*xlam(4)! Re(q21)
                     xqx = xqx + Sqrt(15.0_wp/(8.0_wp*pi))*x*z*xlam(5)! Im(q21)
                     xqx = xqx + Sqrt(15.0_wp/(8.0_wp*pi))*y*z*xlam(6)! extra
                     xqx = xqx*damp(ix, iy, iz)
                     psout(ix, iy, iz, is) = psout(ix, iy, iz, is) + xqx*psin(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      Else If(idp == 2) Then
         Do is = 1, 2
            Do iz = 1, ncolz
               z = xclz(iz)
               Do iy = 1, ncoly
                  y = xcly(iy)
                  Do ix = 1, ncolx
                     x = xclx(ix)
                     xqx = Sqrt(3.0_wp/(4.0_wp*pi))*y*xlam(1)
                     xqx = xqx - Sqrt(3.0_wp/(8.0_wp*pi))*z*xlam(2)
                     xqx = xqx + Sqrt(3.0_wp/(8.0_wp*pi))*x*xlam(3)
                     xqx = xqx - Sqrt(15.0_wp/(8.0_wp*pi))*y*z*xlam(4)
                     xqx = xqx + Sqrt(15.0_wp/(8.0_wp*pi))*x*y*xlam(5)
                     xqx = xqx + Sqrt(15.0_wp/(8.0_wp*pi))*x*z*xlam(6)
                     xqx = xqx*damp(ix, iy, iz)
                     psout(ix, iy, iz, is) = psout(ix, iy, iz, is) + xqx*psin(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      Else If(idp == 3) Then
         Do is = 1, 2
            Do iz = 1, ncolz
               z = xclz(iz)
               Do iy = 1, ncoly
                  y = xcly(iy)
                  Do ix = 1, ncolx
                     x = xclx(ix)
                     xqx = Sqrt(3.0_wp/(4.0_wp*pi))*z*xlam(1)
                     xqx = xqx - Sqrt(3.0_wp/(8.0_wp*pi))*x*xlam(2)
                     xqx = xqx + Sqrt(3.0_wp/(8.0_wp*pi))*y*xlam(3)
                     xqx = xqx - Sqrt(15.0_wp/(8.0_wp*pi))*z*x*xlam(4)
                     xqx = xqx + Sqrt(15.0_wp/(8.0_wp*pi))*z*y*xlam(5)
                     xqx = xqx + Sqrt(15.0_wp/(8.0_wp*pi))*x*y*xlam(6)
                     xqx = xqx*damp(ix, iy, iz)
                     psout(ix, iy, iz, is) = psout(ix, iy, iz, is) + xqx*psin(ix, iy, iz, is)
                  End Do
               End Do
            End Do
         End Do
      End If
!
      Return
End Subroutine costep3
