Subroutine write_density(iseq, time, ncolx, ncoly, ncolz, xclx, xcly, xclz, rho, currnt, clocal, directname)
!-------------------------------------------------------------------------------
!       write data onto a binary file that can then be used for further analysis
!       and plotting. The file name is determined by "iseq", which is passed as
!       either the static iteration or the timestep number.
!-------------------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: iseq
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: time
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: clocal(ncolx, ncoly, ncolz, 2, 2, 3)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: currnt(ncolx, ncoly, ncolz, 3, 2)
      Character(len=80)     :: directname
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Character(10) :: filename
! convert iseq into a string for the filename
      Write(filename, '(I6.6,A4)') iseq, '.tdd'
! create file and write general information
      Open(Unit=99, File='./'//TRIM(directname)//'/tdd'//'/'//filename, Form='UNFORMATTED', Status='REPLACE')
      Write(99) iseq, time, ncolx, ncoly, ncolz
      Write(99) xclx, xcly, xclz
! write density
      Write(99) rho
      Write(99) currnt
      Write(99) clocal!(:,:,:,1,:)
!
End Subroutine write_density
