Subroutine reorder(voc, iorder, n, ns)
!-------------------------------------------------
! subroutine to order energies based index iorder
!-------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: n
      Integer, Intent(In) :: ns
      Integer, Intent(In) :: iorder(n)
      Real(wp), Intent(Inout) :: voc(n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: i, k
!
      k = 0
      voc = 0.0_wp
      Do i = 1, ns
         k = k + 1
         voc(iorder(k)) = 1.0_wp
      End Do
!      write(*,*) "===================="
!      write(*,*) n, ns
!      write(*,*) iorder
!
      Return
End Subroutine reorder
