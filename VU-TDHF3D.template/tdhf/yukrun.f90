Subroutine yukrun(rho, wyuk, pwk1, pwk2, pwk3, pwk4, der2x, der2y, der2z, wxyz, ayuk, ncolx, ncoly, ncolz)
!-----------------------------------------------------------------------
!       solve the Helmholtz equation to obtain the Yukawa potential at
!       collocation points.
!       Details in Umar et al, J. Comp. Phys. 93, 426 (1991).
!
!           output:
!                  wyuk : Yukawa potential at collocation points.
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: ayuk
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Real(wp) :: wyuk(ncolx, ncoly, ncolz)
      Real(wp) :: pwk1(ncolx, ncoly, ncolz)
      Real(wp) :: pwk2(ncolx, ncoly, ncolz)
      Real(wp) :: pwk3(ncolx, ncoly, ncolz)
      Real(wp) :: pwk4(ncolx, ncoly, ncolz)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp), Parameter :: tolmin = 1.0e-10_wp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: iter
      Real(wp) :: zznorm, ztnorm, akp1, zznor1, wdifn, wnorm, ckp1, toler
!----------------------------------------------------------------------
!        setup the right hand vector of the linear system, store in pwk1
!----------------------------------------------------------------------
      pwk1 = - 4.0D0 * pi * ayuk * (rho(:, :, :, 1)+rho(:, :, :, 2))
!----------------------------------------------------------------------
!        main iteration loop begins here
!----------------------------------------------------------------------
      Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, wyuk, pwk2)
!
      pwk1 = pwk1 - pwk2 + wyuk / ayuk ** 2
      pwk2 = pwk1
!
      iter = 0
9999  Continue
      iter = iter + 1
!
      Call laplacian(ncolx, ncoly, ncolz, der2x, der2y, der2z, pwk2, pwk3)
!
      pwk3 = pwk3 - 1.0_wp / ayuk ** 2 * pwk2
!
      zznorm = Sum(wxyz*pwk1*pwk1)
      ztnorm = Sum(wxyz*pwk2*pwk3)
      akp1 = zznorm / ztnorm
!
      pwk4 = wyuk
      wyuk = wyuk + akp1 * pwk2
      pwk4 = pwk4 - wyuk
      pwk1 = pwk1 - akp1 * pwk3
!
      zznor1 = Sum(wxyz*pwk1*pwk1)
      wdifn = Sum(wxyz*pwk4*pwk4)
      wnorm = Sum(wxyz*wyuk*wyuk)
      ckp1 = zznor1 / zznorm
!
      pwk2 = pwk1 + ckp1 * pwk2
!
      toler = wdifn / wnorm
      If(dabs(toler) <= tolmin) Go To 5000
      If(iter > 150) Then
         Write(6, '(2/,A,2/)') ' yukrun: yukrun did not converge'
         Stop
      End If
      Go To 9999
!
!
5000  Continue
      Return
End Subroutine yukrun
