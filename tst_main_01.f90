!
!  Compile with gfortran:
!
!     gfortran -O3 -c mod_mlf_garrappa.f90
!     gfortran -O3 tst_main_01.f90 mod_mlf_garrappa.o -o tst_main_01.exe.gfc
!
!  Compile with ifort:
!
!     ifort -O3 -c mod_mlf_garrappa.f90
!     ifort -O3 tst_main_01.f90 mod_mlf_garrappa.o -o tst_main_01.exe.ifc
!     
!  To run, e.g. the test case c40, 
!
!     ./tst_main_01.exe.gfc cas=c40 eep=6 
!     ./tst_main_01.exe.gfc cas=c40 eep=8 
!     ./tst_main_01.exe.gfc cas=c40 eep=10
!     ./tst_main_01.exe.gfc cas=c40 eep=15
!     ./tst_main_01.exe.gfc cas=c40
!
!
      program tst_main_01
!    
      use mod_mlf_garrappa, only: mlf_garrappa, mld_garrappa, & 
                                  mlf_set_epsilon
!      
      implicit none 
! 
!     Dependences: a simple RNG to generate data for the extra test
!
      real(8) :: auxifun_uniran
!
      integer, parameter :: io=123
!
      complex(8),dimension(:),allocatable :: zj, ej, ee
      real(8) :: afa,bta,argz, dt, rho  
      complex(8) :: z, tmp
      complex(4) :: zsp
      complex(8) :: zfi, zdi, ufi, udi
      real(8) :: rr, ri, rz(4), ffnom, dfnom, ufaer, udaer
      integer  :: fi 
      integer  ::  i, j, k 
      complex(8),parameter :: ic = cmplx(0,1,kind=8)
      logical  :: lg 
      character(len=1) :: ch
      character(len=3) :: cas
      integer :: nr, n0
      character(len=200) :: charg
      character(len=:),allocatable :: chnam, chval 
!
      complex(8),dimension(:,:),allocatable   :: z2d, e2d
      complex(8),dimension(:,:,:),allocatable :: z3d, e3d
      integer :: nd1, nd2, nd3, eep  
!
!
      cas = 'c01'
      eep = 15 
!
      do j = 1, command_argument_count()
         charg = '' 
         call get_command_argument( j, charg )
         call pair_of_name_value ( charg, chnam, chval )
         select case ( trim(adjustl(chnam)) )
         case('cas')
            cas = adjustl(trim(chval))
         case('eep')
            read(chval,*) eep 

         case default 
            cycle 
         end select 
      enddo
!

      rho = 10.0d0**(-eep)

      call mlf_set_epsilon( rho ) 

!      
!     Case cas-th: output from calcualting M-L directly using 
!     the multi-precision package FM.
!     
      inquire ( file='tcases/tt_mlfm_'//cas//'.txt', exist=lg )

      if ( .not. lg ) goto 506

      open(unit=io,file='tcases/tt_mlfm_'//cas//'.txt',action='read')
!
!
      read(io,*) ch, nr, afa, bta, argz

      write(*,100)

      write(*,200) cas, cas 

      write(*,*) 'eep  =', eep 
      write(*,*) 'cas  =', cas 
      write(*,*) 'eps  =', rho 


      write(*,*) ch, nr, afa, bta, argz
      write(*,*) 'afa  =', afa
      write(*,*) 'bta  =', bta
      write(*,*) 'argz =', argz 


      write(*,101)

      write(*,302) 'Robert Garrappa method.'

!
!
      write(*,302) &
      'Calc. Mittag-Leffler function: treat input as 0D (scalar)' 

      call evaltime ( 0, dt )

      ffnom = 0
      ufaer = 0

      do i = 1, nr

         read(io,300,end=505) rr,ri,rz(1),rz(2),rz(3),rz(4)

         zfi = cmplx(rz(1),rz(2),8)

         ffnom = ffnom + abs(zfi)**2

         z = cmplx( rr, ri, kind=8 )

         ufi = mlf_garrappa ( afa, bta, z ) 

         ufaer = ufaer + abs(ufi-zfi)**2

      enddo

505   call evaltime ( 1, dt )

      write(*,*) 'time: ', dt 

      if ( ffnom .eq. 0.0d0 ) ffnom = epsilon(1.0d0)

      if ( i .gt. nr ) then 
         write(*,90)  &
         sqrt(ufaer/nr), sqrt(ufaer/ffnom)
      else 
         if ( i .gt. 1 ) then 
            write(*,90)  &
            sqrt(ufaer/(i-1)), sqrt(ufaer/ffnom)
         else 
            write(*,*) "ERROR: Nothing to be read!"
         endif 
      endif 

!
!
!     Again,
!
      rewind(unit=io)
!
      read(io,*) ch, nr, afa, bta, argz
!
      write(*,302) &
      'Calc. Mittag-Leffler function: treat input as 1D (array)' 

      allocate( zj(1:nr), ej(1:nr), ee(1:nr) )


      call evaltime ( 0, dt )

      ffnom = 0
      ufaer = 0

      do i = 1, nr
         read(io,300,end=605) rr,ri,rz(1),rz(2),rz(3),rz(4)
         zj(i) = cmplx( rr,   ri,   8)
         ee(i) = cmplx( rz(1),rz(2),8)
      enddo

605   continue 

      n0 = nr 
      if ( i .lt. nr )  n0 = i-1
      
      ej = mlf_garrappa (afa,bta, zj)

      do i = 1, n0
         ffnom = ffnom + abs( ee(i) )**2
         ufaer = ufaer + abs( ej(i) - ee(i) )**2
      enddo

      call evaltime ( 1, dt )

      write(*,*) 'time: ', dt 


      if ( ffnom .eq. 0.0d0 ) ffnom = epsilon(1.0d0)
      if ( dfnom .eq. 0.0d0 ) dfnom = epsilon(1.0d0)

      if ( i .gt. nr ) then 
         write(*,90)  &
         sqrt(ufaer/nr), sqrt(ufaer/ffnom)
      else 
         if ( i .gt. 1 ) then 
            write(*,90)  &
            sqrt(ufaer/(i-1)), sqrt(ufaer/ffnom)
         else 
            write(*,*) "ERROR: Nothing to be read!"
         endif 
      endif 


!  
!     Derivative of ML function:
!
!     Again,
!
      rewind(unit=io)

      read(io,*) ch, nr, afa, bta, argz
!
      write(*,302) &
      'Calc. its derivative: treat input as 0D (scalar)' 


      call evaltime ( 0, dt )

      ffnom = 0
      ufaer = 0

      do i = 1, nr

         read(io,300,end=705) rr,ri,rz(1),rz(2),rz(3),rz(4)

         zfi = cmplx(rz(3),rz(4),8)

         ffnom = ffnom + abs(zfi)**2

         z = cmplx( rr, ri, kind=8 )
!
         ufi = mld_garrappa ( afa, bta, z ) 

         ufaer = ufaer + abs(ufi-zfi)**2

      enddo

705   call evaltime ( 1, dt )

      write(*,*) 'time: ', dt 

      if ( ffnom .eq. 0.0d0 ) ffnom = epsilon(1.0d0)

      if ( i .gt. nr ) then 
         write(*,90)  &
         sqrt(ufaer/nr), sqrt(ufaer/ffnom)
      else 
         if ( i .gt. 1 ) then 
            write(*,90)  &
            sqrt(ufaer/(i-1)), sqrt(ufaer/ffnom)
         else 
            write(*,*) "ERROR: Nothing to be read!"
         endif 
      endif 



!
!
!     Again,
!
      rewind(unit=io)
!
      read(io,*) ch, nr, afa, bta, argz
!
      write(*,302) &
      'Calc. its derivative: treat input as 1D (array)' 


      call evaltime ( 0, dt )

      ffnom = 0
      ufaer = 0

      do i = 1, nr
         read(io,300,end=805) rr,ri,rz(1),rz(2),rz(3),rz(4)
         zj(i) = cmplx( rr,   ri,   8)
         ee(i) = cmplx( rz(3),rz(4),8)
      enddo

805   continue 

      n0 = nr 
      if ( i .lt. nr )  n0 = i-1
      
      ej = mld_garrappa (afa,bta, zj)

      do i = 1, n0
         ffnom = ffnom + abs( ee(i) )**2
         ufaer = ufaer + abs( ej(i) - ee(i) )**2
      enddo

      call evaltime ( 1, dt )

      write(*,*) 'time: ', dt 


      if ( ffnom .eq. 0.0d0 ) ffnom = epsilon(1.0d0)
      if ( dfnom .eq. 0.0d0 ) dfnom = epsilon(1.0d0)

      if ( i .gt. nr ) then 
         write(*,90)  &
         sqrt(ufaer/nr), sqrt(ufaer/ffnom)
      else 
         if ( i .gt. 1 ) then 
            write(*,90)  &
            sqrt(ufaer/(i-1)), sqrt(ufaer/ffnom)
         else 
            write(*,*) "ERROR: Nothing to be read!"
         endif 
      endif 












!
!
!
      write(*,302) &
      'For EXTRA-TESTS, treating input as 2D and 3D, read fort.13'


      write(13,302) &
         'Suppose that you now are convinced be the latter performance.'
      write(13,302) &
         'We shall check the consistency of numerical results '
      write(13,302) &
         'for input given randomly as 0D, 1D, 2D and 3D.'


      nd1 = 3
      nd2 = 7
      nd3 = 9 

      allocate ( z2d(nd1,nd2), e2d(nd1,nd2) ) 
      allocate ( z3d(nd1,nd2,nd3), e3d(nd1,nd2,nd3) ) 


      do j=1,nd2
      do i=1,nd1
         z2d(i,j) = cmplx( i*auxifun_uniran(),&
                           j*auxifun_uniran(), kind=8)
      enddo
      enddo


      do k=1,nd3 
      do j=1,nd2
      do i=1,nd1
         z3d(i,j,k) = cmplx( i*auxifun_uniran()*sin(dble(k)),&
                             j*auxifun_uniran()*cos(dble(k)), kind=8)
      enddo
      enddo
      enddo

      

      write(13,304) '2D: Col1 = z(i,j)', 'Col2 = | e(i,j) - "Eexact" |'

      e2d = mlf_garrappa ( afa, bta, z2d )  
      e3d = mlf_garrappa ( afa, bta, z3d )  


      do j=1,nd2
         do i=1,nd1
            ufi = mlf_garrappa ( afa, bta, z2d(i,j) ) 
            write(13,303) z2d(i,j), abs(e2d(i,j)-ufi)
         enddo
         write(13,*)
      enddo


      write(13,304) &
         '3D: Col1 = z(i,j,k)', 'Col2 = | e(i,j,k) - "Eexact" |'


      do k=1,nd3
      do j=1,nd2
         do i=1,nd1
            ufi = mlf_garrappa ( afa, bta, z3d(i,j,k) )
            write(13,303) z3d(i,j,k), abs(e3d(i,j,k)-ufi)
         enddo
         write(13,*)
      enddo
      enddo





!
!
      close(io)
!
      deallocate( zj, ej, ee )
!
      stop 
!
506   write(*,'(a)') 'Not found file: tcases/tt_mlfm_'//cas//'.txt'
      stop 
!
90    format( '-->> absolute error E =', 1pe10.3, 5x, &
                   'relative error E =', 1pe10.3 )
!
100   format( /,72('*') )
101   format( 72('-') )
103   format( 9(a14) )

203   format( 20(e14.6) )

200   format('CASE ',a3,', data file tcases/tt_mlfm_',a3,'.txt')

300   format( 6(1pe23.15e3,1x) )
301   format(35(1pe23.15e3,1x) )
302   format(/,a)
303   format( '(',1pe23.15e3,',',1pe23.15e3,')', 3x, 1x1pe23.15e3 )
304   format( /, a, t55,  a, / )
!
      contains 
!======================================================================
!
!     call evaltime ( 0, dt )
!     ...
!     call evaltime ( 1, dt )
!
!=====
      subroutine evaltime ( istat, dt )
!
      implicit none 
!
      integer,intent(in)  :: istat 
      real(8),intent(out) :: dt 
!
      integer :: clock_rate, clock_max, t  
!
      if ( istat .eq. 0 ) then 
         call system_clock ( t, clock_rate, clock_max )
         dt = dble(t) 
      else 
         call system_clock ( t, clock_rate, clock_max )
         dt =  ( real(t,8) - dt )/ real(clock_rate,8)
      endif 
!
      return 
      end subroutine 
!======================================================================
!     Use some of Fortran 2008 standard.
!
!     ch_input = 'a=1.234'   input 
!     ch_name  = 'a'         output 
!     ch_value = '1.234'     --
!
!=====
      subroutine pair_of_name_value ( ch_input, ch_name, ch_value )
      implicit none 
      character(len=*),intent(in)              :: ch_input 
      character(len=:),allocatable,intent(out) :: ch_name, ch_value 
!
      integer :: k
! 
      do k = 1,len_trim(ch_input) 
         if ( ch_input(k:k) .eq. '=' ) exit 
      enddo
!
      ch_name  = ch_input(:k-1)
      ch_value = ch_input(k+1:)
!
      return 
      end subroutine 
!======================================================================
      end program tst_main_01
