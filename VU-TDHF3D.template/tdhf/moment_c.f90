Subroutine moment_c(time, xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, idp)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: idp
      Real(wp), Intent(In) :: time
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: x, y, z, vol
      Real(wp) :: q21R, q21I, q22R, q22I, q31R, q31I, q32R, q32I, q33R, q33I
      Integer :: iq, ix, iy, iz
!-----------------------------------------------------------------------
!        calculate the cartesian octupole moments
!-----------------------------------------------------------------------
      q21R = 0.0_wp
      q21I = 0.0_wp
      q22R = 0.0_wp
      q22I = 0.0_wp
      q31R = 0.0_wp
      q31I = 0.0_wp
      q32R = 0.0_wp
      q32I = 0.0_wp
      q33R = 0.0_wp
      q33I = 0.0_wp
!
      If(idp == 1) Then
         Do iq = 1, 2
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                     q21R = q21R - Sqrt(15.0_wp/(8.0_wp*pi)) * vol * x * y
                     q21I = q21I + Sqrt(15.0_wp/(8.0_wp*pi)) * vol * x * z
                     q22R = q22R + Sqrt(15.0_wp/(32.0_wp*pi)) * vol * (y**2-z**2)
                     q22I = q22I - Sqrt(15.0_wp/(8.0_wp*pi)) * vol * y * z
                     q31R = q31R + Sqrt(21.0_wp/(64.0_wp*pi)) * vol * y * (y**2+z**2-4.0_wp*x**2)
                     q31I = q31I - Sqrt(21.0_wp/(64.0_wp*pi)) * vol * z * (y**2+z**2-4.0_wp*x**2)
                     q32R = q32R + Sqrt(105.0_wp/(32.0_wp*pi)) * vol * x * (y**2-z**2)
                     q32I = q32I - Sqrt(105.0_wp/(8.0_wp*pi)) * vol * x * y * z
                     q33R = q33R - Sqrt(35.0_wp/(64.0_wp*pi)) * vol * (y**3-3.0_wp*y*z**2)
                     q33I = q33I - Sqrt(35.0_wp/(64.0_wp*pi)) * vol * (z**3-3.0_wp*y**2*z)
                  End Do
               End Do
            End Do
         End Do
      End If
      If(idp == 2) Then
         Do iq = 1, 2
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                     q21R = q21R - Sqrt(15.0_wp/(8.0_wp*pi)) * vol * y * z
                     q21I = q21I + Sqrt(15.0_wp/(8.0_wp*pi)) * vol * x * y
                     q22R = q22R + Sqrt(15.0_wp/(32.0_wp*pi)) * vol * (z**2-x**2)
                     q22I = q22I - Sqrt(15.0_wp/(8.0_wp*pi)) * vol * z * x
                     q31R = q31R + Sqrt(21.0_wp/(64.0_wp*pi)) * vol * z * (z**2+x**2-4.0_wp*y**2)
                     q31I = q31I - Sqrt(21.0_wp/(64.0_wp*pi)) * vol * x * (z**2+x**2-4.0_wp*y**2)
                     q32R = q32R + Sqrt(105.0_wp/(32.0_wp*pi)) * vol * y * (z**2-x**2)
                     q32I = q32I - Sqrt(105.0_wp/(8.0_wp*pi)) * vol * x * y * z
                     q33R = q33R - Sqrt(35.0_wp/(64.0_wp*pi)) * vol * (z**3-3.0_wp*z*x**2)
                     q33I = q33I - Sqrt(35.0_wp/(64.0_wp*pi)) * vol * (x**3-3.0_wp*z**2*x)
                  End Do
               End Do
            End Do
         End Do
      End If
      If(idp == 3) Then
         Do iq = 1, 2
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                     q21R = q21R - Sqrt(15.0_wp/(8.0_wp*pi)) * vol * z * x
                     q21I = q21I + Sqrt(15.0_wp/(8.0_wp*pi)) * vol * z * y
                     q22R = q22R + Sqrt(15.0_wp/(32.0_wp*pi)) * vol * (x**2-y**2)
                     q22I = q22I - Sqrt(15.0_wp/(8.0_wp*pi)) * vol * x * y
                     q31R = q31R + Sqrt(21.0_wp/(64.0_wp*pi)) * vol * x * (x**2+y**2-4.0_wp*z**2)
                     q31I = q31I - Sqrt(21.0_wp/(64.0_wp*pi)) * vol * y * (x**2+y**2-4.0_wp*z**2)
                     q32R = q32R + Sqrt(105.0_wp/(32.0_wp*pi)) * vol * z * (x**2-y**2)
                     q32I = q32I - Sqrt(105.0_wp/(8.0_wp*pi)) * vol * x * y * z
                     q33R = q33R - Sqrt(35.0_wp/(64.0_wp*pi)) * vol * (x**3-3.0_wp*x*y**2)
                     q33I = q33I - Sqrt(35.0_wp/(64.0_wp*pi)) * vol * (y**3-3.0_wp*x**2*y)
                  End Do
               End Do
            End Do
         End Do
      End If
      Write(80, '(F10.4,10E12.4)') time, q21R, q21I, q22R, q22I, q31R, q31I, q32R, q32I, q33R, q33I
!
End Subroutine moment_c
