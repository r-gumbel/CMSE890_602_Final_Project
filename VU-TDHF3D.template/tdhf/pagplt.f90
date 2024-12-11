Subroutine pagplt(f, xar, nx, zar, nz, f1, xperi, xco, zco)
!-----------------------------------------------------------------------
!       to do a 2-d bilinear interpolation of f(x,z).
!       x -- is the horizontal axis
!       z -- is the vertical axis
!       f1   =  value of f which gets printed as "5".
!       xperi= no of x units per horizontal inch of page plot
!       xc,zc= work arrays
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: nx
      Integer, Intent(In) :: nz
      Real(wp), Intent(In) :: f1
      Real(wp), Intent(In) :: xperi
      Real(wp), Intent(In) :: f(nx, nz)
      Real(wp), Intent(In) :: xar(nx)
      Real(wp), Intent(In) :: zar(nz)
      Real(wp), Intent(Inout) :: xco(nx+1)
      Real(wp), Intent(Inout) :: zco(nz+1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: ixsc, izsc, i, nhor, nver, ntkx, ntkz, j, nverm, nhorm, jcar, ntz
      Real(wp) :: dimx, dimz, dxp, dzp, sm, zcu, xcu, bplina
      Character :: ipl, ibr, ids
      Character, Dimension(25) :: icar
      Character, Dimension(121) :: ifun, ibor
!-----------------------------------------------
      Data icar / ' ', '1', ' ', '3', ' ', '5', ' ', '7', ' ', '9', ' ', 'b', ' ', 'd', ' ', 'f', ' ', 'h', ' ', 'j', ' ', 'l', '&
     & ', 'n', '*' /
      Data ipl, ibr, ids, ixsc, izsc / '+', 'i', '-', 10, 6 /
!
      Do i = 1, 121
         ibor(i) = ids
         If(Mod(i-1, 10) /= 0) Cycle
         ibor(i) = ipl
      End Do
!
      dimx = xar(nx) - xar(1)
      dimz = zar(nz) - zar(1)
      dxp = xperi / ixsc
      dzp = xperi / izsc
      nhor = dimx / dxp + 1
      nver = dimz / dzp + 1
      nver = min0(nver, 121)
      ntkx = (nhor+ixsc-1) / ixsc
      ntkz = (nver+izsc-1) / izsc
!
      Write(8, '(A,/,A,F12.4,A)') '1', ' contour 5=', f1, ' n/fm**3'
!
      Do i = 1, ntkx
         xco(i) = xar(1) + (i-1) * ixsc * dxp
      End Do
      Do i = 1, ntkz
         zco(i) = zar(nz) - (i-1) * izsc * dzp
      End Do
!
      Write(8, '(A)') '  z '
      Write(8, '(1X,F6.2,1X,121A)') zco(1), (ibor(j), j=1, nhor)
!
      nverm = nver - 1
      nhorm = nhor - 1
      sm = 5.0_wp / f1
      Do i = 2, nverm
         zcu = zar(nz) - (i-1) * dzp
         Do j = 2, nhorm
            xcu = xar(1) + (j-1) * dxp
            jcar = 1 + Int(0.5_wp+bplina(nx, nz, xar, zar, f, xcu, zcu)*sm)
            jcar = min0(jcar, 25)
            ifun(j) = icar(jcar)
         End Do
         If(Mod(i-1, izsc) == 0) Then
            ntz = (i+izsc-1) / izsc
            ifun(nhor) = ipl
            Write(8, '(1X,F6.2,1X,A,120A)') zco(ntz), '+', (ifun(j), j=2, nhor)
         Else
            ifun(nhor) = ibr
            Write(8, '(8X,A,120A)') 'i', (ifun(j), j=2, nhor)
         End If
      End Do
      If(Mod(nver-1, izsc) == 0) Then
         Write(8, '(1X,F6.2,1X,121A)') zco(ntkz), (ibor(i), i=1, nhor)
      Else
         Write(8, '(8X,A,120A)') '-', (ibor(i), i=2, nhor)
      End If
      Write(8, '(A,12(F6.2,4X),F6.2)') '  x= ', (xco(i), i=1, ntkx)
!
      Return
End Subroutine pagplt
