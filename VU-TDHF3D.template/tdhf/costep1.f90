Subroutine costep1(iq, ncolx, ncoly, ncolz, xclx, xcly, xclz, xcm, ycm, zcm, xlam, damp, psin, psout, idp)
!-----------------------------------------------------------------------
!     quadrupole constraint steps
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: iq
      Integer, Intent(In) :: idp
      Real(wp), Intent(In) :: xlam
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(In) :: damp(ncolx, ncoly, ncolz)
      Complex(wp), Intent(In) :: psin(ncolx, ncoly, ncolz, 2)
      Complex(wp), Intent(Inout) :: psout(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Integer :: is, iz, iy, ix
      Real(wp) :: z, z2, y, y2, x, x2, xq20, cylm
!-----------------------------------------------
      cylm = Sqrt(5.0_wp/(16.0_wp*pi))
      Do is = 1, 2
         Do iz = 1, ncolz
            z = xclz(iz) - zcm(iq)
            z2 = z * z
            Do iy = 1, ncoly
               y = xcly(iy) - ycm(iq)
               y2 = y * y
               Do ix = 1, ncolx
                  x = xclx(ix) - xcm(iq)
                  x2 = x * x
                  If(idp == 1) Then
                     xq20 = cylm * (x2+x2-y2-z2) * xlam
                  Else If(idp == 2) Then
                     xq20 = cylm * (y2+y2-x2-z2) * xlam
                  Else If(idp == 3) Then
                     xq20 = cylm * (z2+z2-x2-y2) * xlam
                  End If
                  psout(ix, iy, iz, is) = psout(ix, iy, iz, is) + damp(ix,iy,iz)*xq20*psin(ix, iy, iz, is)
               End Do
            End Do
         End Do
      End Do
!
      Return
End Subroutine costep1
