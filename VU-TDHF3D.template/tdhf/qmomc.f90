Subroutine qmomc(idp, ncolx, ncoly, ncolz, wx, wy, wz, xclx, xcly, xclz, xcm, ycm, zcm, rho, q20c, q202c)
!-----------------------------------------------------------------------
!        calculate moments q20 and square of q20 for constraint step
!        use the principal axis for calculation
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: idp
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(Out) :: q20c
      Real(wp), Intent(Out) :: q202c
      Real(wp), Parameter :: pi = acos(-1.0_wp)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iq, iz, iy, ix
      Real(wp) :: z, zsq, y, ysq, x, xsq, vol, cylm
!-----------------------------------------------
      q20c = 0.0_wp
      q202c = 0.0_wp
!
      Do iq = 1, 2
         If(idp == 1) Then
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               zsq = z*z
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  ysq = y*y
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     xsq = x*x
                     vol = wx(ix)*wy(iy)*wz(iz)*rho(ix, iy, iz, iq)
                     q20c = q20c + vol*(2.0_wp*xsq-zsq-ysq)
                     q202c = q202c + vol*(2.0_wp*xsq-zsq-ysq)**2
                  End Do
               End Do
            End Do
         Else If(idp == 2) Then
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               zsq = z*z
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  ysq = y*y
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     xsq = x*x
                     vol = wx(ix)*wy(iy)*wz(iz)*rho(ix, iy, iz, iq)
                     q20c = q20c + vol*(2.0_wp*ysq-xsq-zsq)
                     q202c = q202c + vol*(2.0_wp*ysq-xsq-zsq)**2
                  End Do
               End Do
            End Do
         Else If(idp == 3) Then
            Do iz = 1, ncolz
               z = xclz(iz) - zcm(iq)
               zsq = z*z
               Do iy = 1, ncoly
                  y = xcly(iy) - ycm(iq)
                  ysq = y*y
                  Do ix = 1, ncolx
                     x = xclx(ix) - xcm(iq)
                     xsq = x*x
                     vol = wx(ix)*wy(iy)*wz(iz)*rho(ix, iy, iz, iq)
                     q20c = q20c + vol*(2.0_wp*zsq-xsq-ysq)
                     q202c = q202c + vol*(2.0_wp*zsq-xsq-ysq)**2
                  End Do
               End Do
            End Do
         End If
      End Do
!
      cylm = Sqrt(5.0_wp/(16.0_wp*pi))
      q20c = cylm*q20c
      q202c = cylm*cylm*q202c
!
      Return
End Subroutine qmomc
