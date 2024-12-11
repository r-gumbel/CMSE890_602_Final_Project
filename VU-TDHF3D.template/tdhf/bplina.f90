Function bplina(nx, nz, xar, zar, f, xcu, zcu)
!-----------------------------------------------------------------------
!       to do a bilinear inter. of f(nx,nz), on xar(nx) and zar(nz)
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: nx
      Integer, Intent(In) :: nz
      Real(wp), Intent(In) :: xcu
      Real(wp), Intent(In) :: zcu
      Real(wp), Intent(In) :: xar(nx)
      Real(wp), Intent(In) :: zar(nz)
      Real(wp), Intent(In) :: f(nx, nz)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: nxm, nzm, i, icu, jcu
      Real(wp) :: bplina, dxf, dzf, dxc, dzc
!-----------------------------------------------
      bplina = 0.0_wp
      nxm = nx - 1
      nzm = nz - 1
      Do i = 1, nxm
         If(xar(i) <= xcu .And. xcu <= xar(i+1)) Go To 20
      End Do
      Return
20    Continue
      icu = i
      Do i = 1, nzm
         If(zar(i) <= zcu .And. zcu <= zar(i+1)) Go To 40
      End Do
      Return
40    Continue
      jcu = i
      dxf = (xcu-xar(icu)) / (xar(icu+1)-xar(icu))
      dzf = (zcu-zar(jcu)) / (zar(jcu+1)-zar(jcu))
      dxc = 1.0_wp - dxf
      dzc = 1.0_wp - dzf
      bplina = dxf * (f(icu+1, jcu+1)*dzf+f(icu+1, jcu)*dzc) + dxc * (f(icu, jcu+1)*dzf+f(icu, jcu)*dzc)
!
      Return
End Function bplina
