Subroutine getin(h2m, t0, t1, t2, t3, t4, t4p, x0, x1, x2, x3, alpha, ayuk, vyuk, v0neut, v0prot, rho0pr, &
& radinf, itrbx, mtrbx, serr, derr, x0dmp, e0dmp, iprint, mprint, mplots, mplott, icoul, irest, mrests, mrestt, &
& imode, ifixcm, itimrev, itimrevs, ipairf, c0, d0, iconstr, q2in, mconstr, itheta, nt, dt, mnof, nof, centf, boostf, &
& euler_alpha, euler_beta, euler_gamma, fmass, fcharg, ecm, rsep, xb, mxp, terr, ifixb, irwgs, iodds, ihdiag, &
& nexadd, nexiter, nextra_n, nextra_p, directname, isrestart)
!-----------------------------------------------------------------------
!       getin  - subroutine to read in parameters
!
!       the parameters defining the calculation are as follows:
!
!       BSPL.INP (read in from tdhf3d.f90 using grid routine)
!       Read in mesh/b-spline parameters
!
!       PAIR.INP (contains VDI AND/OR DDDI pairing strengts etc.)
!
!       SKYRME.INP
!       parameters for the Skyrme potential:
!                   t0,..,t4  -  the standard Skyrme parameters
!                   x0,x3     -  the exchange parameters
!                   alpha     -  the power of the density dependence
!                                in the t3,x3 term
!                                zero-range t3 term: x3=alpha=1.0
!                   ayuk,vyuk -  range and strength of the original BKN force
!                   itheta    -  include J**2 terms
!
!       TDHF3D.INP
!                   FORCE     -  name of the skyrme interaction used (see skyrme.inp for names)
!                   nof       -  number of nuclei (maximum mnof in tdhf3d.f90)
!                   imode     -  static/dynamic control (0=static, 1=dynamic)
!                  centf(3,i) -  location of the i th nucleus (if ifixb=0 or nof/=2)
!                                1=x , 2=y, 3=z
!                 boostf(3,i) -  boost for the i th nucleus (dynamic run only)
!              euler_alpha(i) -  euler alpha for each nuclei (i=1,nof)
!               euler_beta(i) -  euler beta for each nuclei (i=1,nof)
!              euler_gamma(i) -  euler gamma for each nuclei (i=1,nof)
!                                Note:  Euler angles are relative to original orientation
!                    fmass(i) -  mass of the i th nucleus
!                   fcharg(i) -  charge of the i th nucleus
!       parameters for wave function initialization for static runs:
!                 radinf(3,i) -  spreads in x,y,z directions of the Gaussian for the
!                                initial guess for the ith nucleus.
!                                this is usually in the range of 2.5-3.5. Small
!                                numbers may result in NaN's for large boxes.
!
!       parameters tuning the iteration:
!                   itrbx     -  maximum number of static iterations
!                   mtrbx     -  maximum number of iterations for density constraint
!                   serr      -  error in generating the bnd wfns
!                   derr      -  error in generating the bnd wfns for density constraint
!                   x0dmp     -  overall factor in front of the
!                                damping operator.
!                   e0dmp     -  energy by which the kinetic energy
!                                is divided by. this regulates the
!                                strength of the suppression of high
!                                Fourier components in the state.
!                   ifixcm    -  include the (1-1/A) cm fixing term in the k.e.
!                                1=yes. only for static calculations.
!                   icoul     -  turn Coulomb on and off
!                                icoul=0 - Coulomb and exchange off
!                                icoul=1 - Coulomb and exchange on
!                   iconstr   -  constraint flag (1=quadrupole, 2=denstiy constraint)
!                                (3=quadrupole and q10, q21, 4=quadrupole, q10, q21, q22 axial)
!                   q2in      -  the quadrupole moment to be constrained to this value
!                   mconstr   -  how often to do density constraint (w.r.t. time steps)
!                   c0,d0     -  density constraint parameters
!       the following have same meaning for static/dynamic:
!                   iprint    -  print level. standard output is 0.
!                                higher values produce more output.
!                   mprint    -  print modulus
!                   mplots,t  -  plot modulus for static and time-dependent runs
!                   irest     -  flag for the restart option
!                                =0 run is not a restart
!                                =n restart from nth iteration
!                   mrests    -  how often to store static restart file
!                   mrestt    -  how often to store dynamic restart file
!       time reversal control
!                   itimrev   - =0 Time Reversal Symmetry not
!                                  explicitly incorporated. Otherwise
!                                  # of states cut in 2 via TR sym.
!                                Always =0 for dynamic calculations.
!                   itimrevs  -  TR symmetry for static calculations carried
!                                out in dynamic mode (density constraint with even-even
!                                system, itimrevs=1, odd-A system itimrevs=0 always)
!                   iodds     -  turn on/off all the time-odd spin terms, in which case
!                                itheta also determines which terms are included:
!                                if itheta=1 and iodds=1 all terms are present
!                                if itheta=0 and iodds=1 (s.T-J**2) terms omitted
!                                if itheta=1 and iodds=0 only J's in B_q included
!                                NOTE: subrotines densit, energy, and skyrme contain an
!                                      internal parameter ioddls=1. This leaves the odd
!                                      part of l-s force even when iodds=0. For completely
!                                      turning off odd terms one should have ioddls=0.
!                   ihdiag    -  H-diagonalization flag (0=off, 1= on)
!       pairing calculation
!              ipairf =       - pairing flags for each nucleus
!                               =0 no pairing
!                               =1 VDI pairing (when available)
!                               =2 DDDI pairing (when available)
!                               =3 DDDI pairing with Fermi cutoff from below
!                               =4 DDDI pairing with Fermi cutoff from above and below
!       Extra states for static runs
!               nexadd          =0 or 1 for adding extra states to either nucleus
!               nexiter         =# number of iterations after which reordering of levels stop
!               nextra_n(i)     = extra neutron levels for the ith nucleus (i=1,nof)
!               nextra_p(i)     = extra proton levels for the ith nucleus (i=1,nof)
!
!       dynamic input parameters:
!                   nt        -  number of time steps
!                   dt        -  time step in fm/c
!                   ifixb     -  flag for dynamic determination of
!                                boosts for 2 fragments, using ecm
!                                and rsep. (=1 yes)
!                   irwgs     -  read g.s. wavefunctions for time-evolution from restart
!                                e.g. running same impact parameter at different Ecm's
!                   ecm       -  center of mass energy for 2 fragments
!                   rsep      -  relative separation for 2 fragments
!                   xb        -  initial impact parameter (Fermi)
!                   mxp       -  number of exponential terms
!                   terr      -  error in the norm of each s.p. state
!                                for the expansion of propagator
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer :: itrbx
      Integer :: mtrbx
      Integer :: iprint
      Integer :: mprint
      Integer :: mplots
      Integer :: mplott
      Integer :: icoul
      Integer :: irest
      Integer :: irwgs
      Integer :: mrests
      Integer :: mrestt
      Integer :: imode
      Integer :: itimrev
      Integer :: iodds
      Integer :: ihdiag
      Integer :: itimrevs
      Integer :: iconstr
      Integer :: mconstr
      Integer :: itheta
      Integer :: nt
      Integer :: mnof
      Integer :: nof
      Integer :: mxp
      Integer :: ifixb
      Integer :: nexadd
      Integer :: nexiter
      Real(wp) :: h2m(2)
      Real(wp) :: t0
      Real(wp) :: t1
      Real(wp) :: t2
      Real(wp) :: t3
      Real(wp) :: t4
      Real(wp) :: t4p
      Real(wp) :: x0
      Real(wp) :: x1
      Real(wp) :: x2
      Real(wp) :: x3
      Real(wp) :: alpha
      Real(wp) :: ayuk
      Real(wp) :: vyuk
      Real(wp) :: serr
      Real(wp) :: derr
      Real(wp) :: x0dmp
      Real(wp) :: e0dmp
      Real(wp) :: c0
      Real(wp) :: d0
      Real(wp) :: q2in
      Real(wp) :: ecm, rsep, xb, terr, dt
      Real(wp), Dimension(3, mnof) :: centf, boostf, radinf
      Real(wp), Dimension(mnof)    :: fmass, fcharg, euler_alpha, euler_beta, euler_gamma
      Integer, Dimension(2, mnof)  :: ipairf
      Integer, Dimension(mnof)     :: nextra_n, nextra_p
      Logical, Dimension(mnof)     :: isrestart
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: i, ifixcm, mass, ibeta
      Real(wp) :: v0prot_vdi, v0neut_vdi, rho0pr_vdi, v0prot_dddi, v0neut_dddi, rho0pr_dddi
      Real(wp) :: v0prot, v0neut, rho0pr
      Logical :: exists
      Character(8) :: force_type, tipo
      Character(len=80) :: directname, a1_c, a2_c, z1_c, z2_c, ecm_c, xb_c, filename, it_c, beta_c
      Character(len=80) :: par11_c, par12_c, par21_c, par22_c
      Character(35), Parameter :: version = "         VU-TDHF3D@ Version 2020.05"
!
      Read(5,*) force_type
      Read(5,*) nof
      If(nof > mnof) Then
         Write(*, '(2/,A,2/)') ' getin: nof > mnof '
         Stop
      End If
      Read(5,*) imode
      Read(5,*) (centf(1, i), i=1, nof)
      Read(5,*) (centf(2, i), i=1, nof)
      Read(5,*) (centf(3, i), i=1, nof)
!      if(imode == 0) centf = 0.0_wp
      Read(5,*) (boostf(1, i), i=1, nof)
      Read(5,*) (boostf(2, i), i=1, nof)
      Read(5,*) (boostf(3, i), i=1, nof)
      Read(5,*) (euler_alpha(i), i=1, nof)
      Read(5,*) (euler_beta(i), i=1, nof)
      Read(5,*) (euler_gamma(i), i=1, nof)
      Read(5,*) (fmass(i), i=1, nof)
      Read(5,*) (fcharg(i), i=1, nof)
      Read(5,*) (radinf(1, i), i=1, nof)
      Read(5,*) (radinf(2, i), i=1, nof)
      Read(5,*) (radinf(3, i), i=1, nof)
      Read(5,*) x0dmp, e0dmp
      Read(5,*) itrbx, mtrbx, serr, derr
      Read(5,*) iprint, mprint, mplots, mplott
      Read(5,*) irest, irwgs, mrests, mrestt
      Read(5,*) ifixcm, icoul
      Read(5,*) iconstr, q2in, mconstr, c0, d0
      Read(5,*) itimrev, itimrevs, iodds, ihdiag
      Read(5,*) (ipairf(1, i), i=1, nof)
      Read(5,*) (ipairf(2, i), i=1, nof)
      Read(5,*) nexadd, nexiter
      Read(5,*) (nextra_n(i), i=1, nof)
      Read(5,*) (nextra_p(i), i=1, nof)
! read force parameters from skyrme.inp using force name from tdhf3d.inp
      Do
         Read(Unit=11, Advance="no", Fmt='(A)', Eor=20) tipo
         If(TRIM(tipo) == TRIM(force_type)) Then
            Read(11,*) t0, t1, t2, t3, x0, x1, x2, x3, t4, t4p, alpha, ayuk, vyuk, h2m(1), h2m(2), itheta
            Go To 30
         End If
         Read(11,*)
         Read(Unit=11, Advance="no", Fmt='(A)', end=20)
      End Do
20    Continue
      Stop 'Unknown type of Skyrme force in file skyrme.inp'
30    Continue
! read pairing parameters from pair.inp using force name from tdhf3d.inp
      Do
         Read(Unit=10, Advance="no", Fmt='(A)', Eor=40) tipo
         If(TRIM(tipo) == TRIM(force_type)) Then
            Read(10,*) v0prot_vdi, v0neut_vdi, rho0pr_vdi, v0prot_dddi, v0neut_dddi, rho0pr_dddi
            Go To 50
         End If
         Read(10,*)
         Read(Unit=10, Advance="no", Fmt='(A)', end=40)
      End Do
40    Continue
      Stop 'Unknown type of pairing force in file pair.inp'
50    Continue
!
      If(imode == 1) Then
         Read(5,*) ifixb
         If(nof /= 2 .And. ifixb == 1) Then
            Write(8, '(2/,A,2/)') ' getin: boost fix for nof = 2 only'
            Stop
         End If
         Read(5,*) ecm, rsep, xb
!        impact parameter is always positive (right ion up from axis)
         xb = Abs(xb)
!        avoid exactly zero b for numerical reasons for coulomb trajectory
         If(xb == 0.0_wp .And. nof == 2) Then
            xb = 0.000010_wp
         End If
         Read(5,*) mxp, terr
         Read(5,*) nt, dt
      End If
! setup the directory name for output files
      Write(unit=a1_c,fmt='(i3)') NINT(fmass(1))
      Write(unit=z1_c,fmt='(i3)') NINT(fcharg(1))
      Write(unit=a2_c,fmt='(i3)') NINT(fmass(2))
      Write(unit=z2_c,fmt='(i3)') NINT(fcharg(2))
      Write(unit=it_c,fmt='(i1)') itimrevs
      If(imode == 0) Write(unit=it_c,fmt='(i1)') itimrev
      Write(unit=ecm_c,fmt='(f10.3)') ecm
      Write(unit=xb_c,fmt='(f10.3)') xb
      ibeta = euler_beta(2)
      Write(unit=beta_c,fmt='(i3)') ibeta
      Write(unit=par11_c,fmt='(i1)') ipairf(1,1)
      Write(unit=par12_c,fmt='(i1)') ipairf(2,1)
      Write(unit=par21_c,fmt='(i1)') ipairf(1,2)
      Write(unit=par22_c,fmt='(i1)') ipairf(2,2)
      a1_c = Adjustl(a1_c)
      a2_c = Adjustl(a2_c)
      z1_c = Adjustl(z1_c)
      z2_c = Adjustl(z2_c)
      ecm_c = Adjustl(ecm_c)
      beta_c = Adjustl(beta_c)
      xb_c = Adjustl(xb_c)
      it_c = Adjustl(it_c)
      par11_c = Adjustl(par11_c)
      par12_c = Adjustl(par12_c)
      par21_c = Adjustl(par21_c)
      par22_c = Adjustl(par22_c)
      directname = "AL_"//TRIM(a1_c)//"_ZL_"//TRIM(z1_c)//"_AR_"//TRIM(a2_c)//"_ZR_"//TRIM(z2_c)//"_"//&
      &TRIM(force_type)//"_Ecm"//TRIM(ecm_c)//"_beta"//TRIM(beta_c)//"_b"//TRIM(xb_c)
      call createdir(TRIM(directname))
      call createdir(TRIM(directname)//"/tdd")
!-----------------------------------------------------------------------
!        open files for output
!-----------------------------------------------------------------------
      Open(unit=8, File='./'//TRIM(directname)//'/'//'tdhf3d.out',Status='unknown', Form='formatted', Position='asis')
      Write(8, '(/,a,/)') version
      Open(Unit=12, File='./'//TRIM(directname)//'/'//'shf3d.restart', Status='unknown', Form='unformatted', Position='asis')
      Open(Unit=13, File='./'//TRIM(directname)//'/'//'tdhf3d.plot', Status='unknown', Form='formatted', Position='asis')
      Open(Unit=14, File='./'//TRIM(directname)//'/'//'tdhf3d.restart', Status='unknown', Form='unformatted', Position='asis')
      Open(Unit=15, File='./'//TRIM(directname)//'/'//'t0gs.restart', Status='unknown', Form='unformatted', Position='asis')
      Open(Unit=78, File='./'//TRIM(directname)//'/'//'moments_isoscalar.out', Status='unknown', Form='formatted', Position='asis')
      Open(Unit=79, File='./'//TRIM(directname)//'/'//'moments_isovector.out', Status='unknown', Form='formatted', Position='asis')
      Open(Unit=80, File='./'//TRIM(directname)//'/'//'moments_complex.out', Status='unknown', Form='formatted', Position='asis')
      isrestart = .FALSE.
      If (imode == 0) Then
         filename="A_"//TRIM(a1_c)//"_Z_"//TRIM(z1_c)//"_pair_"//TRIM(par11_c)//"_"//TRIM(par12_c)//"_itim_"//TRIM(it_c)//"_"//TRIM(force_type)//".static"
         Open(Unit=19, File=TRIM(filename), Status='unknown', Form='unformatted', Position='asis')
      Else If (NINT(fmass(1))==NINT(fmass(2)) .AND. NINT(fcharg(1))==NINT(fcharg(2)) .AND. ipairf(1,1)==ipairf(2,1) .AND. ipairf(1,2)==ipairf(2,2)) Then
         filename="A_"//TRIM(a1_c)//"_Z_"//TRIM(z1_c)//"_pair_"//TRIM(par11_c)//"_"//TRIM(par12_c)//"_itim_"//TRIM(it_c)//"_"//TRIM(force_type)//".static"
         Inquire(FILE=filename,EXIST=exists)
         If (exists .AND. irwgs == 0) Then
            nextra_n(1) = 0
            nextra_p(1) = 0
            nextra_n(2) = 0
            nextra_p(2) = 0
            Write(8,'(/A/)') "Restarting from file for fragment 1 and 2."
            isrestart(1) = .TRUE.
            isrestart(2) = .TRUE.
            Open(Unit=19, File=TRIM(filename), Status='unknown', Form='unformatted', Position='asis')
         End If
      Else
         filename="A_"//TRIM(a1_c)//"_Z_"//TRIM(z1_c)//"_pair_"//TRIM(par11_c)//"_"//TRIM(par12_c)//"_itim_"//TRIM(it_c)//"_"//TRIM(force_type)//".static"
         Inquire(FILE=filename,EXIST=exists)
         If (exists .AND. irwgs == 0) Then
            nextra_n(1) = 0
            nextra_p(1) = 0
            Write(8,'(/A/)') "Restarting from file for fragment 1."
            isrestart(1) = .TRUE.
            Open(Unit=19, File=TRIM(filename), Status='unknown', Form='unformatted', Position='asis')
         End If
         filename="A_"//TRIM(a2_c)//"_Z_"//TRIM(z2_c)//"_pair_"//TRIM(par21_c)//"_"//TRIM(par22_c)//"_itim_"//TRIM(it_c)//"_"//TRIM(force_type)//".static"
         Inquire(FILE=filename,EXIST=exists)
         If (exists .AND. irwgs == 0) Then
            nextra_n(2) = 0
            nextra_p(2) = 0
            Write(8,'(/A/)') "Restarting from file for fragment 2."
            isrestart(2) = .TRUE.
            Open(Unit=20, File=TRIM(filename), Status='unknown', Form='unformatted', Position='asis')
         End If
      End If
!
      If(force_type /= "BKN  " .And. ayuk >= 1.0e-5_wp) Then
         Write(8, '(/A/)') ' Force parameter ayuk should be zero if force not BKN'
         Stop
      End If
!
      If(imode == 0 .And. iconstr == 2) Then
         Write(8, '(/A/)') ' No density contraint for a static calculation is implemented'
         Write(8, '(/A/)') ' Only density consraint during TDHF is allowed'
         Stop
      End If
      If(imode == 1 .And. ifixcm == 1) Then
         Write(8, '(/A/)') ' Incompatible c.m. correction for tdhf'
         Stop
      End If
      If(imode == 1 .And. itimrev == 1) Then
         Write(8, '(/A/)') ' Incompatible TDHF - time-reversal input'
         Stop
      End If
      Do i = 1, nof
         mass = Int(fmass(i)) - Int(fcharg(i))
         If(Mod(mass, 2) /= 0 .And. itimrev == 1 .And. imode == 0) Then
            Write(8, '(/A/)') ' Odd n nucleus not TR invariant (itimrev)'
            Stop
         End If
         If(Mod(mass, 2) /= 0 .And. itimrevs == 1 .And. imode == 1) Then
            Write(8, '(/A/)') ' Odd n nucleus not TR invariant (itimrevs)'
            Stop
         End If
         mass = Int(fcharg(i))
         If(Mod(mass, 2) /= 0 .And. itimrev == 1 .And. imode == 0) Then
            Write(8, '(/A/)') ' Odd charge nucleus not TR invariant (itimrev)'
            Stop
         End If
         If(Mod(mass, 2) /= 0 .And. itimrevs == 1 .And. imode == 1) Then
            Write(8, '(/A/)') ' Odd charge nucleus not TR invariant (itimrevs)'
            Stop
         End If
      End Do
      If(itheta == 1 .And. imode == 1 .And. iodds == 0) Then
         Write(8, '(/A/)') ' (s.T-J**2) Galilean invariance violated'
         Stop
      End If
!
      Write(8, '(2/,A,A8,A)') ' skyrme ', force_type, ' parameters'
      Write(8, '(A,F15.5)') ' t0           =  ', t0
      Write(8, '(A,F15.5)') ' t1           =  ', t1
      Write(8, '(A,F15.5)') ' t2           =  ', t2
      Write(8, '(A,F15.5)') ' t3           =  ', t3
      Write(8, '(A,F15.5)') ' t4           =  ', t4
      Write(8, '(A,F15.5)') ' t4p          =  ', t4p
      Write(8, '(A,F15.5)') ' x0           =  ', x0
      Write(8, '(A,F15.5)') ' x1           =  ', x1
      Write(8, '(A,F15.5)') ' x2           =  ', x2
      Write(8, '(A,F15.5)') ' x3           =  ', x3
      Write(8, '(A,F15.5)') ' alpha        =  ', alpha
      Write(8, '(A,F15.5)') ' ayuk(BKN)    =  ', ayuk
      Write(8, '(A,F15.5)') ' vyuk(BKN)    =  ', vyuk
      Write(8, '(A,I15  )') ' theta-ls     =  ', itheta
      Write(8, '(A,I15  )') ' iodds        =  ', iodds
!
      If(iconstr == 1) Write(8, '(2/,A,F15.5,/)') ' quadrupole moment constrained to: ', q2in
      If(iconstr == 3) Write(8, '(2/,A,F15.5,/)') ' quadrupole moment constrained to: ', q2in
      If(iconstr == 4) Write(8, '(2/,A,F15.5,/)') ' quadrupole moment constrained to: ', q2in
      If(iconstr == 3) Write(8, '(A,/)') ' also constrain q10 and q21 '
      If(iconstr == 4) Write(8, '(A,/)') ' also constrain q10, q21, and q22 (axial) '
      If(iconstr == 2) Write(8, '(2/,A,I4,A,/)') ' perform density constraint every ', mconstr, ' time-steps'
!
      Write(8, '(/,A,I12,/)') ' number of nuclei =  ', nof
      Do i = 1, nof
         Write(8, '(/,A,I2,/)') ' initial data for nucleus ', i
         Write(8, '(A,F12.4)') ' mass               =  ', fmass(i)
         Write(8, '(A,F12.4)') ' charge             =  ', fcharg(i)
         Write(8, '(A,F12.4)') ' euler angle alpha  =  ', euler_alpha(i)
         Write(8, '(A,F12.4)') ' euler angle beta   =  ', euler_beta(i)
         Write(8, '(A,F12.4)') ' euler angle gamma  =  ', euler_gamma(i)
         Write(8, '(A,1P,3E12.4)') ' origin (x,y,z)     =  ', centf(1, i), centf(2, i), centf(3, i)
         Write(8, '(A,1P,3E12.4)') ' boost  (x,y,z)     =  ', boostf(1, i), boostf(2, i), boostf(3, i)
         Write(8, '(A,3F12.4)') ' osc. widths        =  ', radinf(1, i), radinf(2, i), radinf(3, i)
!
         If(ipairf(1, i) /= 0 .Or. ipairf(2, i) /= 0) Then
            Select Case(ipairf(1, i))
            Case(1)
               Write(8, '(/,A,/)') ' Neutrons use VDI pairing'
               v0neut = v0neut_vdi
               v0prot = v0prot_vdi
               rho0pr = rho0pr_vdi
               Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
               Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
               Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
               If(v0neut < 0.0_wp) Then
                  Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                  flush(8)
                  STOP 'Pairing of this type not available for this force'
               End If
            Case(2)
               Write(8, '(/,A,/)') ' Neutrons use DDDI pairing'
               v0neut = v0neut_dddi
               v0prot = v0prot_dddi
               rho0pr = rho0pr_dddi
               Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
               Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
               Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
               If(v0neut < 0.0_wp) Then
                  Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                  flush(8)
                  STOP 'Pairing of this type not available for this force'
               End If
             Case(3)
                Write(8, '(/,A,/)') ' Neutrons use DDDI pairing with Fermi cutoff from above'
                v0neut = v0neut_dddi
                v0prot = v0prot_dddi
                rho0pr = rho0pr_dddi
                Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
                Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
                Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
                If(v0neut < 0.0_wp) Then
                   Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                   flush(8)
                   STOP 'Pairing of this type not available for this force'
                End If
             Case(4)
                Write(8, '(/,A,/)') ' Neutrons use DDDI pairing with Fermi cutoff from above and below'
                v0neut = v0neut_dddi
                v0prot = v0prot_dddi
                rho0pr = rho0pr_dddi
                Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
                Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
                Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
                If(v0neut < 0.0_wp) Then
                   Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                   flush(8)
                   STOP 'Pairing of this type not available for this force'
                End If
            Case Default
               Write(8, '(/,A)') ' No pairing for neutrons'
            End Select
            Select Case(ipairf(2, i))
            Case(1)
               Write(8, '(/,A,/)') ' Protons use VDI pairing'
               v0neut = v0neut_vdi
               v0prot = v0prot_vdi
               rho0pr = rho0pr_vdi
               Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
               Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
               Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
               If(v0prot < 0.0_wp) Then
                  Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                  flush(8)
                  STOP 'Pairing of this type not available for this force'
               End If
            Case(2)
               Write(8, '(/,A,/)') ' Protons use DDDI pairing'
               v0neut = v0neut_dddi
               v0prot = v0prot_dddi
               rho0pr = rho0pr_dddi
               Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
               Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
               Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
               If(v0prot < 0.0_wp) Then
                  Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                  flush(8)
                  STOP 'Pairing of this type not available for this force'
               End If
             Case(3)
                Write(8, '(/,A,/)') ' Protons use DDDI pairing with Fermi cutoff from above'
                v0neut = v0neut_dddi
                v0prot = v0prot_dddi
                rho0pr = rho0pr_dddi
                Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
                Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
                Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
                If(v0neut < 0.0_wp) Then
                   Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                   flush(8)
                   STOP 'Pairing of this type not available for this force'
                End If
             Case(4)
                Write(8, '(/,A,/)') ' Protons use DDDI pairing with Fermi cutoff from above and below'
                v0neut = v0neut_dddi
                v0prot = v0prot_dddi
                rho0pr = rho0pr_dddi
                Write(8, '(A,F15.5)') ' v0neut           =  ', v0neut
                Write(8, '(A,F15.5)') ' v0prot           =  ', v0prot
                Write(8, '(A,F15.5)') ' rho0pr           =  ', rho0pr
                If(v0neut < 0.0_wp) Then
                   Write(8, '(/,A,/)') 'Pairing of this type not available for this force'
                   flush(8)
                   STOP 'Pairing of this type not available for this force'
                End If
            Case Default
               Write(8, '(/,A)') ' No pairing for protons'
            End Select
         Else If(ipairf(1, i) == 0 .And. ipairf(2, i) == 0) Then
            Write(8, '(2/,A)') ' pairing not used'
         Else
            Write(8, '(2/,A)') ' Wrong input value for pairing flag'
            Stop
         End If
         If(((ipairf(1, i) /= 0) .or. (ipairf(2, i) /= 0)) .and. (nexadd /= 0)) Then
             Write(8,'(/,a,/)') "Pairing and extra added states conflict"
             Stop
         End If
         If(itimrev == 0 .And. (ipairf(1, i) /= 0 .Or. ipairf(2, i) /= 0) .And. imode /= 0) Then
            Write(8, '(/,A)') ' Using fixed partial occupations in TDHF'
         End If
         If(itimrevs /= 1 .And. (ipairf(1, i) /= 0 .Or. ipairf(2, i) /= 0) .And. imode == 1) Then
            Write(8, '(//,A)') ' Incompatible pairing - time-reversal input for static calculations'
!            Stop
         End If
         If(itimrev /= 1 .And. (ipairf(1, i) /= 0 .Or. ipairf(2, i) /= 0) .And. imode == 0) Then
            Write(8, '(//,A)') ' Incompatible pairing - time-reversal input for static calculations'
!            Stop
         End If
      End Do
!
      If(imode == 1 .And. iconstr == 0) Then
         Write(8, '(2/,A)') ' damping parameters and flags:'
         Write(8, '(A,1P,E12.4)') ' damping coefficient         =  ', x0dmp
         Write(8, '(A,1P,E12.4)') ' damping energy scale        =  ', e0dmp
         Write(8, '(A,1P,E12.4)') ' energy fluctuation error    =  ', serr
         Write(8, '(A,I12)')      ' no. of static iterations    =  ', itrbx
         Write(8, '(A,I12)')      ' no. of constraint iterations=  ', mtrbx
      End If
      If(imode == 0 .And. iconstr == 0) Then
         Write(8, '(2/,A)') ' damping parameters and flags:'
         Write(8, '(A,1P,E12.4)') ' damping coefficient         =  ', x0dmp
         Write(8, '(A,1P,E12.4)') ' damping energy scale        =  ', e0dmp
         Write(8, '(A,1P,E12.4)') ' energy fluctuation error    =  ', serr
         Write(8, '(A,I12)')      ' no. of static iterations    =  ', itrbx
         Write(8, '(A,I12)')      ' no. of constraint iterations=  ', mtrbx
      End If
      If(imode == 1 .And. iconstr /= 0) Then
         Write(8, '(2/,A)') ' damping parameters and flags:'
         Write(8, '(A,1P,E12.4)') ' damping coefficient         =  ', x0dmp
         Write(8, '(A,1P,E12.4)') ' damping energy scale        =  ', e0dmp
         Write(8, '(A,1P,E12.4)') ' energy fluctuation error    =  ', serr
         Write(8, '(A,1P,E12.4)') ' energy fluctuation error(dc)=  ', derr
         Write(8, '(A,I12)')      ' no. of static iterations    =  ', itrbx
         Write(8, '(A,I12)')      ' no. of constraint iterations=  ', mtrbx
         Write(8, '(A,1P,E12.4)') ' constraint amplitude (c0)   =  ', c0
         Write(8, '(A,1P,E12.4)') ' constraint exchange  (d0)   =  ', d0
      End If
      If(imode == 0 .And. iconstr /= 0) Then
         Write(8, '(2/,A)') ' damping parameters and flags:'
         Write(8, '(A,1P,E12.4)') ' damping coefficient         =  ', x0dmp
         Write(8, '(A,1P,E12.4)') ' damping energy scale        =  ', e0dmp
         Write(8, '(A,1P,E12.4)') ' energy fluctuation error    =  ', serr
         Write(8, '(A,1P,E12.4)') ' energy fluctuation error(dc)=  ', derr
         Write(8, '(A,I12)')      ' no. of static iterations    =  ', itrbx
         Write(8, '(A,I12)')      ' no. of constraint iterations=  ', mtrbx
         Write(8, '(A,1P,E12.4)') ' constraint amplitude (c0)   =  ', c0
         Write(8, '(A,1P,E12.4)') ' constraint exchange  (d0)   =  ', d0
      End If
      Write(8, '(2/,A,I12)') ' print flag                  =  ', iprint
      Write(8, '(A,I12)') ' print modulus               =  ', mprint
      Write(8, '(A,I12)') ' static plot modulus         =  ', mplots
      Write(8, '(A,I12)') ' dynamic plot modulus        =  ', mplott
      Write(8, '(A,I12)') ' restart flag                =  ', irest
      Write(8, '(A,I12)') ' static restart store freq.  =  ', mrests
      Write(8, '(A,I12)') ' dynamic restart store freq. =  ', mrestt
      Write(8, '(A,I12)') ' c.m. correction flag        =  ', ifixcm
      Write(8, '(A,I12)') ' coulomb flag                =  ', icoul
!
      If(imode == 1) Then
         Write(8, '(2/,A)')    ' dynamic parameters and flags:'
         Write(8, '(A,I12)')   ' number of time steps =  ', nt
         Write(8, '(A,F12.5)') ' time step (fm/c)     =  ', dt
         Write(8, '(A,I12)')   ' dynamic boost flag   =  ', ifixb
         Write(8, '(A,F12.4)') ' ecm  (dynamic)  (mev)=  ', ecm
         Write(8, '(A,F12.4)') ' rsep (dynamic)   (fm)=  ', rsep
         Write(8, '(A,F12.4)') ' xb                   =  ', xb
         Write(8, '(A,I12)')   ' mxp                  =  ', mxp
         Write(8, '(A,1P,E12.4)') ' terr (norm error)    =  ', terr
!
      End If
!
      Return
!
End Subroutine getin
