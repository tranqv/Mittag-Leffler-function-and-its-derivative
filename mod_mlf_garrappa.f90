!=======================================================================
!
!  1) The following Fortran code is developed mainly from the Matlab 
!     code of Prof. Robert Garrappa with a minor modification. Feel free
!     to use or modify it. In doing so, please cite the following 
!     references:
!
!     [1] https://www.mathworks.com/matlabcentral/fileexchange/
!                                      48154-the-mittag-leffler-function 
!         (the Matlab code)
!
!     [2] R. Garrappa, Numerical evaluation of two and three parameter
!         Mittag-Leffler functions, SIAM Journal of Numerical Analysis, 
!         2015, 53(3), 1350-1369
!         (Theoretical background)
!
!     [3] Gorenflo, Rudolf, Joulia Loutchko, and Yuri Luchko. 
!         "Computation of the Mittag-Leffler function E alpha, beta (z) 
!         and its derivative." Fract. Calc. Appl. Anal. 2002.
!         (For calculation of derivative of the Mittag-Leffler function)
!
!     [4] "An Overview of Software Development for Special Functions", 
!         W. J. Cody, Lecture Notes in Mathematics, 506, Numerical 
!         Analysis Dundee, 1975, G. A. Watson (ed.), Springer Verlag, 
!         Berlin, 1976.
!         (For a portable gamma function with under/overflow trapping.)
!
!  2) Please report on this Fortran code to to tranquocviet@tdt.edu.vn 
!     or viet204@gmail.com. 
!
!=======================================================================
!     Version 02, dated Tue May 1 2018 (current version)
!     
!     +  passed more than 70 test cases fantastically. The test cases 
!        are released together with this code, however, with a rough
!        explanation. I'm going to update more testcases.
!
!     +  compiler for the test cases: gfortran version 4.9.2, ifort 
!        version 18.0.0.
!
!     +  several comments of the Matlab code are still remained in this
!        version. Most of them are marked by !@.
!
!     +  there are lots of irrelevant comments inside this version. I 
!        don't have time to clean. So just ignore or clean them by
!        yourself :D 
!
!     +  I have packed everything into only one module to not conflict 
!        with other package. Now the code is really portable and safe.
!
!     +  what's else? ...
!
!=======================================================================
!     gfortran -O3 -c mod_mlf_garrappa.f90
!
      module mod_mlf_garrappa
!
      implicit none 
!      
      private 
!
!-----------------------------------------------------------------------
!
!     group 1: fixed  
!
      integer,parameter :: r4 = kind(1.0e0), r8 = kind(1.0d0)
      integer,parameter :: i4 = selected_int_kind(9)
      integer,parameter :: i8 = selected_int_kind(18)
!      
!     group 2: choose precision
!
      integer,parameter :: rk = r8  ! double precision
!     integer,parameter :: rk = r4  ! single precision (not available)
      integer,parameter :: ik = i4  ! integers with 4 bytes, common use
!
!     group 3: machine dependent constants
!
      integer(ik),parameter :: iinf = huge(iinf)
      real(rk),parameter    :: rinf = huge(rinf)
      real(rk),parameter    :: rtin = tiny(rtin)
      real(rk),parameter    :: reps = epsilon(reps)
      real(rk),parameter    :: lnep = log(reps)
      real(rk),parameter    :: xbig = 171.624_rk      ! if rk=r8 
!     real(rk),parameter    :: xbig =  35.040_rk      ! if rk=r4 
!
      real(rk),parameter ::                                      &!
         picons= 3.141592653589793238462643383279502884197e+0_rk,&!pi
         piinve= 3.183098861837906715377675267450287240689e-1_rk,&!1/pi
         pidiv2= 1.570796326794896619231321691639751442099e+0_rk,&!pi/2
         pipow2= 9.869604401089358618834490999876151135314e+0_rk,&!pi^2
         pimul2= 6.283185307179586476925286766559005768394e+0_rk,&!pi*2
         sqrtpi= 1.772453850905516027298167483341145182798e+0_rk,&!pi^0.5
         loge10= 2.302585092994045684017991454684364207601e+0_rk  !ln10
!         
      complex(rk),parameter :: iz = cmplx(0.0_rk,1.0_rk,rk)
!
      real(rk),parameter :: default_epsilon = 1.0e-15_rk  
      real(rk),save      :: present_epsilon = default_epsilon
!
!-----------------------------------------------------------------------
!
!     To reset the prepenst_epsilon above on demand. If we do not call 
!     this, present_epsilon is set to default_epsilon = 10^(-15).
!
      public :: mlf_set_epsilon
!----------------------------------------------------------------------
!
!     General Mittag Leffler function for various kinds of input:
!
      public :: genmlf
      interface genmlf
         module procedure genmlf_garrappa_01 
         module procedure genmlf_garrappa_02 
         module procedure genmlf_garrappa_03 
         module procedure genmlf_garrappa_04 
      end interface 
!
!     Usage: 
!
!        E = genmlf ( afa, bta, gma, z )
!
!     for z and E are defined as scalar (0D) or arrays (1D, 2D, or 3D).
!
!----------------------------------------------------------------------
!
!     Mittag Leffler function for various shapes of input:
!
      public :: mlf_garrappa
      interface mlf_garrappa
         module procedure mlf_garrappa_01 
         module procedure mlf_garrappa_02 
         module procedure mlf_garrappa_03 
         module procedure mlf_garrappa_04 
      end interface 
!
!     Usage: 
!
!        E = mlf_garrappa ( afa, bta, z )
!
!     for z and E are defined as scalar (0D) or arrays (1D, 2D, or 3D).
!
!----------------------------------------------------------------------
!
!     Derivative of Mittag Leffler function for shapes of input:
!
      public :: mld_garrappa
      interface mld_garrappa
         module procedure mld_garrappa_01
         module procedure mld_garrappa_02
         module procedure mld_garrappa_03
         module procedure mld_garrappa_04
      end interface 
!      
!----------------------------------------------------------------------
      contains 
!======================================================================
      function mlf_garrappa_01 ( afa, bta, z ) result(E)
!
      implicit none 
!
      complex(rk),intent(in)  :: z
      real(rk),intent(in)     :: afa, bta
      complex(rk)             :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
!     Local variables:
!
      complex(rk) :: zloc(2)
!      
      zloc(1) = z
!      
      call sub_genmlf_multishot (                      &
           afa, bta, 1.0_rk, 1, 1, zloc(1), 1, zloc(2) )
!           
      e = zloc(2) 
!     
      return 
      end function
!=====      
!
!     1D input
!
      function mlf_garrappa_02 ( afa, bta, z ) result(E)
!
      implicit none 
!
      real(rk),intent(in)                 :: afa, bta
      complex(rk),dimension(:),intent(in) :: z
      complex(rk),dimension(size(z))      :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
      call sub_genmlf_multishot (afa,bta,1.0_rk,size(z),1,z,1,E)
!
      return 
      end function
!=====
!
!     2D input
!
      function mlf_garrappa_03 ( afa, bta, z ) result(E)
!
      implicit none 
!
      real(rk),intent(in)                        :: afa, bta
      complex(rk),dimension(:,:),intent(in)      :: z
      complex(rk),dimension(size(z,1),size(z,2)) :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
      call sub_genmlf_multishot (                    &
           afa, bta, 1.0_rk,                         &
           size(z,1)*size(z,2), 1, z(:,1), 1, E(:,1) )
!
      return 
      end function
!=====     
!     3D input
!
      function mlf_garrappa_04 ( afa, bta, z ) result(E)
!
      implicit none 
!
      real(rk),intent(in)                     :: afa, bta
      complex(rk),dimension(:,:,:),intent(in) :: z
      complex(rk),dimension(size(z,1),size(z,2),size(z,3)) :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
      call sub_genmlf_multishot (                                  &
           afa,  bta, 1.0_rk,                                      &
           size(z,1)*size(z,2)*size(z,3), 1, z(:,1,1), 1, E(:,1,1) )
!
      return 
      end function
!======================================================================
!
!     General ML function:
!
      function genmlf_garrappa_01 ( afa, bta, gma, z ) result(E)
!
      implicit none 
!
      complex(rk),intent(in)  :: z
      real(rk),intent(in)     :: afa, bta, gma
      complex(rk)             :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
!     Local variables:
!
      complex(rk) :: zloc(2)
!      
      zloc(1) = z
!      
      call sub_genmlf_multishot (                   &
           afa, bta, gma, 1, 1, zloc(1), 1, zloc(2) )
!           
      e = zloc(2) 
!     
      return 
      end function
!=====     
!
!     1D input
!
      function genmlf_garrappa_02 ( afa, bta, gma, z ) result(E)
!
      implicit none 
!
      real(rk),intent(in)                 :: afa, bta, gma
      complex(rk),dimension(:),intent(in) :: z
      complex(rk),dimension(size(z))      :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
      call sub_genmlf_multishot (afa,bta,gma,size(z),1,z,1,E)
!
      return 
      end function
!=====
!
!     2D input
!
      function genmlf_garrappa_03 ( afa, bta, gma, z ) result(E)
!
      implicit none 
!
      real(rk),intent(in)                        :: afa, bta, gma
      complex(rk),dimension(:,:),intent(in)      :: z
      complex(rk),dimension(size(z,1),size(z,2)) :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
      call sub_genmlf_multishot (                    &
           afa, bta, gma,                            &
           size(z,1)*size(z,2), 1, z(:,1), 1, E(:,1) )
!
      return 
      end function
!=====     
!     3D input
!
      function genmlf_garrappa_04 ( afa, bta, gma, z ) result(E)
!
      implicit none 
!
      real(rk),intent(in)                     :: afa, bta, gma
      complex(rk),dimension(:,:,:),intent(in) :: z
      complex(rk),dimension(size(z,1),size(z,2),size(z,3)) :: E
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot
!
      call sub_genmlf_multishot (                                  &
           afa,  bta, gma,                                         &
           size(z,1)*size(z,2)*size(z,3), 1, z(:,1,1), 1, E(:,1,1) )
!
      return 
      end function
!=======================================================================
!
!     Derivative of Mittag-Leffer function with two parameters. 
!
!     Applying THEOREM 4.1 of the Gorenflo's paper
!
!=====
      function mld_garrappa_01 ( afa, bta, z )  result(f)
!
      implicit none 
!
      real(rk),intent(in)    :: afa, bta
      complex(rk),intent(in) :: z
      complex(rk)            :: f  
!
!     Dependencies:
!
!     external :: sub_dermlf_multishot
!
!     Local variables:
!
      complex(rk) :: zloc(2)
!
!     Check if the arguments afa>0. Otherwise, do nothing and report. 
!
      if ( afa .le. 0.0_rk ) then 
         write(*,'(a,1pe10.3)') &
            'ERROR: Wrong input to mld_garrappa_01, afa<=0,',afa
         return 
      endif 
!
      zloc(1) = z
!
      call sub_dermlf_multishot (              &
           afa, bta, 1, 1, zloc(1), 1, zloc(2) )
!
      f = zloc(2)
!
      return 
      end function
!=====
      function mld_garrappa_02 ( afa, bta, z )  result(f)
!
      implicit none 
!
      real(rk),intent(in)                 :: afa, bta
      complex(rk),dimension(:),intent(in) :: z
      complex(rk),dimension(size(z))      :: f 
!
!     Dependencies:
!
!     external :: sub_dermlf_multishot
!
!     Check if the arguments afa>0. Otherwise, do nothing and report. 
!
      if ( afa .le. 0.0_rk ) then 
         write(*,'(a,1pe10.3)') &
            'ERROR: Wrong input to mld_garrappa_02, afa<=0,',afa
         return 
      endif 
!
      call sub_dermlf_multishot ( afa, bta, size(z), 1, z, 1, f )
!
      return 
      end function
!=====
      function mld_garrappa_03 ( afa, bta, z )  result(f)
!
      implicit none 
!
      real(rk),intent(in)                   :: afa, bta
      complex(rk),dimension(:,:),intent(in) :: z
      complex(rk),dimension(size(z,1),size(z,2)) :: f 
!
!     Dependencies:
!
!     external :: sub_dermlf_multishot
!
!     Check if the arguments afa>0. Otherwise, do nothing and report. 
!
      if ( afa .le. 0.0_rk ) then 
         write(*,'(a,1pe10.3)') &
            'ERROR: Wrong input to mld_garrappa_03, afa<=0,',afa
         return 
      endif 
!
      call sub_dermlf_multishot (                              &
           afa, bta, size(z,1)*size(z,2), 1, z(:,1), 1, f(:,1) )
!
      return 
      end function
!=====
      function mld_garrappa_04 ( afa, bta, z )  result(f)
!
      implicit none 
!
      real(rk),intent(in)                     :: afa, bta
      complex(rk),dimension(:,:,:),intent(in) :: z
      complex(rk),dimension(size(z,1),size(z,2),size(z,3)) :: f 
!
!     Dependencies:
!
!     external :: sub_dermlf_multishot
!
!     Check if the arguments afa>0. Otherwise, do nothing and report. 
!
      if ( afa .le. 0.0_rk ) then 
         write(*,'(a,1pe10.3)') &
            'ERROR: Wrong input to mld_garrappa_03, afa<=0,',afa
         return 
      endif 
!
      call sub_dermlf_multishot (                                  &
           afa, bta,                                               &
           size(z,1)*size(z,2)*size(z,3), 1, z(:,1,1), 1, f(:,1,1) )
!
      return 
      end function

!======================================================================
!
!     Set the value of epsilon in the range: 4.44e-16 < eps < 1.0e-1 
!
!=====
      subroutine mlf_set_epsilon (eps)
      real(rk),intent(in) :: eps 
!     
      if ( 2*reps .lt. eps .and. eps .lt. 1.0e-1_rk ) then 
         present_epsilon = eps 
      endif 
!      
      return 
      end subroutine 
!======================================================================
!
!  Description coppied from the Matlab code of Prof. R. Garrappa:
!
! Evaluation of the Mittag-Leffler (ML) function with 1, 2 or 3 parameters
! by means of the OPC algorithm [1]. The routine evaluates an approximation
! Et of the ML function E such that |E-Et|/(1+|E|) approx 1.0e-15   
!     
!
! E = ML(z,afa) evaluates the ML function with one parameter afa for
! the corresponding elements of z; afa must be a real and positive
! scalar. The one parameter ML function is defined as
!
! E = sum_{k=0}^{infty} z^k/Gamma(afa*k+1)
!
! with Gamma the Euler's gamma function.
!
!
! E = ML(z,afa,bta) evaluates the ML function with two parameters afa
! and bta for the corresponding elements of z; afa must be a real and
! positive scalar and bta a real scalar. The two parameters ML function is
! defined as
!
! E = sum_{k=0}^{infty} z^k/Gamma(afa*k+bta)
!
!
! E = ML(z,afa,bta,gma) evaluates the ML function with three parameters
! afa, bta and gma for the corresponding elements of z; afa must be a
! real scalar such that 0<afa<1, bta any real scalar and gma a real and
! positive scalar; the arguments z must satisfy |Arg(z)| > afa*pi. The
! three parameters ML function is defined as 
!
! E = sum_{k=0}^{infty} Gamma(gma+k)*z^k/Gamma(gma)/k!/Gamma(afa*k+bta)
!
!
! NOTE: 
! This routine implements the optimal parabolic contour (OPC) algorithm
! described in [1] and based on the inversion of the Laplace transform on a
! parabolic contour suitably choosen in one of the regions of analyticity
! of the Laplace transform.
!
!
! REFERENCES
!
!   [1] R. Garrappa, Numerical evaluation of two and three parameter
!   Mittag-Leffler functions, SIAM Journal of Numerical Analysis, 2015,
!   53(3), 1350-1369
!
!
!   Please, report any problem or comment to : 
!        roberto dot garrappa at uniba dot it
!
!   Copyright (c) 2015, Roberto Garrappa, University of Bari, Italy
!   roberto dot garrappa at uniba dot it
!   Homepage: http://www.dm.uniba.it/Members/garrappa
!   Revision: 1.4 - Date: October 8 2015
!
!======================================================================
!
!     Here we are! This is the steersman.
!
!=====
!
      subroutine sub_genmlf_multishot (             &
                 afa, bta, gma, n, incz, z, ince, e )
!
!     use mod_mlf_garrappa
!
      implicit none 
!
!     Input & Output: 
!
      real(rk),intent(in)     :: afa, bta, gma 
      integer(ik),intent(in)  :: n, incz, ince 
      complex(rk),dimension(incz,n),intent(in)  :: z
      complex(rk),dimension(ince,n),intent(out) :: e
!
!     Dependencies:
!
!     external :: worksub_ltinv_multishot 
!     real(rk) :: specfun_gamma
!     external :: auxisub_initrandom
!
!     work-space: 
!
      complex(rk),dimension(:),allocatable :: zws  ! length = 2*m
      real(rk),dimension(:),allocatable    :: rws  ! length = 6*m+1
      integer(ik),dimension(:),allocatable :: iws  ! length = 2*m
!
      integer(ik) :: m 
!     
!     Local variables:
!
      complex(rk) :: zinp, eout 
      real(rk)    :: eps, log_eps, epsmul10  
      real(rk)    :: t, absz, argz, aagz
      integer(ik) :: j
      logical     :: gma_is_not_one
!
!     HERE WE GO 
!
      if ( afa .le. 0.0_rk .or. gma .le. 0.0_rk ) then 
         write(*,907) 
         return 
      endif 
!
      gma_is_not_one = abs(gma-1) .gt. reps
!
      if ( gma_is_not_one .and. (afa > 1.0_rk) ) then 
         write(*,908) 
         return 
      endif 
!
!     afa=bta=gma=1
!
      if ( (.not. gma_is_not_one) .and. &
            abs(afa-1) .lt. reps  .and. &
            abs(bta-1) .lt. reps        ) then 
         do j = 1,n      
            e(1,j) = exp(z(1,j))
         enddo 
         return 
      endif
!
!     Initializing RNG for the quick sort:
!
      call auxisub_initrandom
!
!     User defines epsilon or just use default values, i.e. 1.0e-15.
!     This epsilon_value variable is maintained in the module above.
!     Its value can be reset by, e.g., call mlf_set_epsilon( 1.0d-10 )
!
      eps      = present_epsilon
!  
!     eps      = 1.0e-15_rk   ! the most accurate setting.
!     eps      = reps * 2     ! value smaller than reps gives errors
!     eps      = 1.0e-07_rk   ! with a larger values it runs faster.
!
      log_eps  = log(eps)  
      epsmul10 = eps * 10 
!
!     Other values of t=1 may not working. Hence, let us fix t=1.
!
      t = 1.0_rk 
!
!     Preparing for work-space:
!
      m = int(afa) + 2  
!
      allocate( zws(1:2*m), rws(1:6*m+1), iws(1:2*m) ) 
!
!     Let's fly 
!
      do j = 1,n      
!
!        Input z(*):
!
         zinp = z(1,j)
!         
         absz = abs( zinp )
         argz = atan2( aimag(zinp), real(zinp) ) 
         aagz = abs( argz )
!
!        Check parameters and arguments for the three parameter case
!
         if ( gma_is_not_one .and. ( aagz <= afa*picons ) ) then 
            write(*,909) 
            return 
         endif 
!
!        Inversion of the LT for each element of z
!
         if ( absz .lt. eps ) then 
!
            eout = 1.0_rk / specfun_gamma(bta) 
!
         else
!
            zws(:) = cmplx(0.0_rk,0.0_rk,rk)
            rws(:) = 0.0_rk
            iws(:) = 0 
!
            call worksub_ltinv_multishot (            &
                 t, zinp, afa, bta, gma,              &
                 absz, argz,                          &
                 eps, log_eps, epsmul10, m,           &
                 zws(1), zws(m+1), rws(1), rws(m+1),  &
                 rws(2*m+2), rws(3*m+2),              &
                 rws(4*m+2), rws(5*m+2), iws(1),      &
                 iws(m+1), eout                       )
!
         endif 
!
!        Output:
!
         e(1,j) = eout 
!
      enddo
!
!     Dismiss:
!
      deallocate( zws, rws, iws ) 
!
      return 
!
907   format('ERROR: ml requires afa>0 and gma>0' ) 
908   format('ERROR: |gma-1|>eps, ml requires 0 < ALPHA < 1' ) 
909   format('ERROR: |gma-1|>eps, ml requires |Arg(z)| > afa*pi.')
!
      end subroutine 
!
!======================================================================
!
      subroutine worksub_ltinv_multishot (          &
                 t, z, afa, bta, gma,               &
                 absz, argz,                        &
                 eps, log_eps, epsmul10, m,         &
                 ztmp, s_star, rtmp,  phi_s_star,   &
                 p, q,  array_mu, array_h, array_n, &
                 admissible_regions, E              )
!
!     use mod_mlf_garrappa
!      
      implicit none 
!
!     Input & Output:
!
      complex(rk),intent(out) :: E
      real(rk),intent(in)     :: t
      complex(rk),intent(in)  :: z
      real(rk),intent(in)     :: afa, bta, gma
      real(rk),intent(in)     :: absz, argz      
      real(rk),intent(in)     :: eps, log_eps, epsmul10 
!
!     Workspace: 
!
      integer,intent(in)       :: m
      complex(rk),dimension(m) :: ztmp
      complex(rk),dimension(m) :: s_star 
      real(rk),dimension(m)    :: rtmp 
      real(rk),dimension(m+1)  :: phi_s_star 
      real(rk),dimension(m)    :: p, q 
      real(rk),dimension(m)    :: array_mu
      real(rk),dimension(m)    :: array_h
      integer(ik),dimension(m) :: array_n
      integer(ik),dimension(m) :: admissible_regions 
!
!     Dependencies:
!
!     external :: primsub_optimalparam_rb
!
!     Local variables:
!
      complex(rk) :: e_integral, e_residues 
      complex(rk) :: zk, zd, fk
      real(rk)    :: rtm1, uk, muj, hj, mu, h, local_logeps
      integer(ik) :: jj1, j1, nj, n, idmin, klp1  
      integer(ik) :: kmin, kmax, k_vett 
      integer(ik) :: j, k, klen, n_admissible_regions 
      integer(ik) :: len_s_star, len_s_star_p1
      logical     :: not_found_region 
!
!     GO
!
      kmin = ceiling(-afa/2 - argz/pimul2 )
      kmax = floor  ( afa/2 - argz/pimul2 )
!      
      klen = kmax-kmin+1
      klp1 = klen+1
!
!     NOTE on the workspace:
!
!     Most of actual workspaces require the length klen+1, except 
!     the array phi_s_star, which needs the length klen+2. 
!
!     Most of the provided workspaces have length m, except the 
!     array of the phi_s_star that is coming with length m+1. 
!
!     What is m ? 
!
!           m = int(afa) + 2
!
!     PROOF of m+1 >= klen+2 
!
!     For any x in R, we have 
!
!           kmin = ceiling(-afa/2 - x ) >= -afa/2 - x
!           kmax = floor  ( afa/2 - x ) <=  afa/2 - x
!
!     Hence 
!
!           kmax-kmin <=  afa/2 - x - (-afa/2 - x) = afa < int(afa) + 1
!
!     for any afa in R. Note that the inequality afa < int(afa)+1 
!     does not allow the equality "=" to be happened. Hence, we have
!
!           kmax-kmin < int(afa) + 1
!
!     Looking to both the sides, they all are integers. Therefore we 
!     conclude that
!
!           kmax-kmin <= int(afa)   holds for any afa in R.
!
!     Wow ... so far, we conclude that 
!
!           klen = kmax-kmin+1 <= int(afa)+1
!     and
!           klp1 = klen+1      <= int(afa)+2 = m  (As defined above) 
!
!     Or, m+1 >= klen+2 (DONE). 
!
!     In this case, the longest workspace phi_s_star should have length
!     of m+1, while others have that of m, where m = int(afa) + 2. For 
!     each input afa, we have the workspaces to works with an array of 
!     input z(:).
!
!
!@Evaluation of phi(s_star) for each pole
!
      rtm1 = absz**(1.0_rk/afa)
!
      do k = 1, klen
         k_vett = kmin - 1 + k      
         s_star(k+1) = rtm1 * exp(iz*( argz + pimul2*k_vett )/afa)
         phi_s_star(k+1) = (real(s_star(k+1),rk) + abs(s_star(k+1)))/2
      enddo
!
!@Sorting of the poles according to the value of phi(s_star)
!
!
      if ( klen .gt. 1 ) then 
!
         call primsub_qsort_r8_idx ( klen, phi_s_star(2), array_n(1) )
!
!        Then, phi_s_star is sorted only by its indices. To get the task
!        to be done, we perform
!
         do k = 1,klen
            rtmp(k) = phi_s_star(array_n(k)+1) 
         enddo
!         
         do k = 1,klen
            phi_s_star(k+1) = rtmp(k)  
         enddo
!
!     Then we arrange the s_star(:) according to the order of
!     phi_s_star(:). 
!
         do k = 1,klen
            ztmp(k) = s_star(array_n(k)+1) 
         enddo
!         
         do k = 1,klen
            s_star(k+1) = ztmp(k)  
         enddo
!
      endif 
!
!@Deleting possible poles with phi_s_star=0 ...
!
      len_s_star = 0
!
      do while ( phi_s_star(2) .le. eps .and. len_s_star .lt. klen )
!
         len_s_star = len_s_star + 1 
!
!        Shift the arrays from-right-to-left for one index 
!
         do k = 2,klen 
            phi_s_star(k) = phi_s_star(k+1) 
            s_star(k)     = s_star(k+1) 
         enddo
!
!        until the fisrt element phi_s_star(1) > eps. In this case,
!        we have phi_s_star(k) >=  phi_s_star(1) > eps, for k>1,
!        since phi_s_star(:) is increased.
!
      enddo
!
      len_s_star = klen - len_s_star 
!
      len_s_star_p1 = len_s_star + 1
!
      phi_s_star(1) = 0.0_rk 
      s_star(1)     = cmplx(0.0_rk,0.0_rk,rk)
!
!@Strength of the singularities ...
!
!     Now the efficient length of phi_s_star = len_s_star+2, because 
!        phi_s_star = [ 0, phi_s_star, +Inf] ;
!
      p(1) = max( 0.0_rk, -2*(afa*gma-bta+1) )
!
      do k = 2,len_s_star_p1
         p(k) = gma
      enddo
!
      do k = 1,len_s_star
         q(k) = gma
      enddo
!
      q(len_s_star_p1) = rinf 
!
      phi_s_star(len_s_star_p1+1) = rinf 
!
!     So far, length(phi_s_star) must be >= len_s_star_p1 + 1 = klen+2
!
!@Looking for the admissible regions with respect to round-off errors
!
      do j = 1,klp1
         array_n(j) = -1
      enddo

      local_logeps = log_eps 
 
      rtm1 = ( local_logeps - lnep ) / t
!
      n_admissible_regions = 0 
!
      do k = 1, len_s_star_p1
         if (  phi_s_star(k) .lt. rtm1  .and.     &
               phi_s_star(k) .lt. phi_s_star(k+1) ) then 
            array_n(k) = k
            n_admissible_regions = n_admissible_regions + 1
         endif 
      enddo
!
      if ( n_admissible_regions .gt. 0 ) then 

         j = 0 

         do k = 1, len_s_star_p1
            if ( array_n(k) .gt. 0 ) then
               j = j + 1
               admissible_regions(j) = array_n(k)
            endif 
         enddo
!
         jj1 = admissible_regions( n_admissible_regions ) 
!
      else 
!
         n_admissible_regions = 1
         jj1 = 1 
!
      endif 
!
      do j = 1,klp1
         array_mu(j) = rinf  
         array_h(j)  = rinf  
         array_n(j)  = iinf  
      enddo
!
!@Evaluation of parameters for inversion of LT in each admissible region
!
      not_found_region = .true.

      do while ( not_found_region )

         do j = 1, n_admissible_regions

            j1 = admissible_regions (j) 

            if ( j1 < len_s_star_p1 ) then 
               call primsub_optimalparam_rb (             & 
                    t, phi_s_star(j1), phi_s_star(j1+1),  &
                    p(j1), q(j1), local_logeps, epsmul10, &
                    muj, hj, nj                           ) 

            else
               call primsub_optimalparam_ru (      &
                    t, phi_s_star(j1),             &
                    p(j1), local_logeps, epsmul10, &
                    muj, hj, nj  ) 

            endif 

            array_mu(j1) = muj  
            array_h (j1) = hj  
            array_n (j1) = nj 

         enddo
!
         n = minval(array_n)
!
         if ( n > 200 ) then 
            local_logeps = local_logeps + log(10.0_rk)
         else
            not_found_region = .false.
         endif 

      enddo
!
!@Selection of the admissible region for integration which involves the
! minimum number of nodes 
!
      idmin = 0
      n     = iinf 
!
      do k = 1, jj1
         if ( n .gt. array_n (k) ) then 
            n = array_n (k)
            idmin = k 
         endif 
      enddo
!
      mu = array_mu(idmin) 
      h  = array_h (idmin) 
!
!     Alright, from now on everything is transparent.
!
!@Evaluation of the inverse Laplace transform (herein z:=lambda)
!
      e_integral = cmplx(0.0_rk,0.0_rk,rk) 
!
      do k = -n, n
!
         uk = h*k 
         zk = mu*(iz*uk + 1)**2 
         zd = 2*mu*(iz-uk)
!
         fk = ( zk**(afa*gma-bta) / ((zk**afa - z)**gma) )*zd
!
         e_integral = e_integral + fk*exp(zk*t)
!
      enddo
!
      e_integral = h* e_integral / pimul2 / iz 
!
!@Evaluation of residues
!
      e_residues = cmplx(0.0_rk,0.0_rk,rk) 

      do k = idmin+1, len_s_star_p1
         e_residues = e_residues + &
                    (1/afa)*( s_star(k)**(1-bta) )*exp(t*s_star(k))
      enddo
!
!@Evaluation of the ML function
!
      E = e_integral + e_residues 
!
      if (aimag(z) .eq. 0.0_rk) E = cmplx(real(E,rk),0.0_rk,rk)
!
      return 
      end subroutine 
!=======================================================================
!
!     This routine sorts the array only by its index, i.e. idord. 
!     After calling, the ARRAY IS STILL UNCHANGED. To sort it, swap
!        array(j) <-> array(idord(j)) for all j.
!
!=====
!
      subroutine primsub_qsort_r8_idx ( n, array, idord )
!     
!     use mod_mlf_garrappa
!      
      implicit none 
!
      integer(ik),intent(in)   :: n 
      real(rk),dimension(n)    :: array
      integer(ik),dimension(n) :: idord
!
!     Dependence: 
!
!     external :: primsub_r8_quick_sort
!
!     local variables:
!
      integer(ik) :: i
!
      if ( n .lt. 2 ) return 
!
      do i = 1, n
         idord(i) = i
      enddo
! 
      call primsub_r8_quick_sort ( n, array, idord, 1, n )
!
      return 
      end subroutine
!=====
      recursive subroutine primsub_r8_quick_sort (      &
                           n, array, idord, left, right )
!                           
!     use mod_mlf_garrappa
!
      implicit none 
!
      integer(ik),intent(in)   :: n 
      real(rk),dimension(n)    :: array
      integer(ik),dimension(n) :: idord
      integer(ik),intent(in)   :: left
      integer(ik),intent(in)   :: right
!
!     Dependence: 
!
!     real(4) :: auxifun_uniran
!
!     local variables:
!
      integer(ik) :: i, last, itmp 
      integer(ik) :: ichoose
!
!     GO 
!
      if ( left .ge. right ) return
! 
!     ichoose returns an integer ramdomly in the set 
!     of integers {left,left+1,...,right}
!     (Make sure that we initialized the RNG already, in DRIVERs)
!      
      ichoose = left + int(real(right-left+1)*auxifun_uniran())
!      
      itmp           = idord(left) 
      idord(left)    = idord(ichoose) 
      idord(ichoose) = itmp 
!
      last = left
! 
      do i = left+1, right
         if (  array(idord(i)) .lt. array(idord(left)) ) then
            last = last + 1
            itmp        = idord(last)
            idord(last) = idord(i)
            idord(i)    = itmp 
         endif
      enddo
!
      itmp        = idord(left)
      idord(left) = idord(last) 
      idord(last) = itmp 
!     
      call primsub_r8_quick_sort( n, array, idord, left,  last-1 )
      call primsub_r8_quick_sort( n, array, idord, last+1, right )
!     
      end subroutine 
!======================================================================
!
      subroutine auxisub_initrandom
      implicit none 
      integer,dimension(:),allocatable :: gieo
      integer :: i, n
      logical,save :: rng_not_yet_init = .true.
!
      if (rng_not_yet_init) then 
         call random_seed(size=n)
         allocate( gieo(n) )
         call system_clock(count=gieo(1))           
         do i = 2,n
            gieo(i) = gieo(i-1) + iand(gieo(i-1),31) + 1 
         enddo
         call random_seed(put=gieo)
         deallocate( gieo )
         rng_not_yet_init = .false.
      endif 
      return 
      end subroutine  
!=====
      function auxifun_uniran() result(r) 
      implicit none 
      real(4) :: r
!
      call random_number( r )
!
      return 
      end function 
!======================================================================
      function specfun_gamma (x)  result(f)
!
!     use mod_mlf_garrappa 
!
!----------------------------------------------------------------------
!
! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants.  Let
!
! BETA   - radix for the floating-point representation
! MAXEXP - the smallest positive power of beta that overflows
!
! Then the following machine-dependent constants must be declared 
!   in DATA statements.  IEEE values are provided as a default.
!
! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - the largest machine representable floating-point number;
!          approximately BETA**MAXEXP, XINF = huge(XINF)
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0, EPS = epsilon(EPS)
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable, XMININ = tiny(XMININ)
!
!     Approximate values for some important machines are:
!
!                            beta       maxexp        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624 <--*
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489
!
!
!***                       huge(.)     epsilon(.)   tiny(.)
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! Error returns (This is what I need: ERROR TRAPPING)
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!
!  Intrinsic functions required are:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.
!
!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.
!
!  Latest modification: March 12, 1992
!
!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!  Modified by Tran Quoc Viet, tranquocviet@tdt.edu.vn
!              Sat Apr 28 17:13:38 +07 2018
!
!----------------------------------------------------------------------
!
!     Input & Output:
!
      real(rk) :: f 
      real(rk),intent(in) :: x 
!
!     Local variables:
!
      integer(ik) :: i, n
      logical     :: logipar
      real(rk)    :: fact, res, sss, xden, xnum, y, y1, ysq, z
!
!----------------------------------------------------------------------
!     Mathematical constants  (CLEANED)
!----------------------------------------------------------------------
!s    data one,half,twelve,two,zero/1.0e0,0.5e0,12.0e0,2.0e0,0.0e0/, &
!s         sqrtpi/0.9189385332046727417803297e0/,                    &
!s         pi/3.1415926535897932384626434e0/
!
!     real(rk),parameter ::                                    &
!     picons= 3.141592653589793238462643383279502884197e+0_rk, &!pi
!     sqrtpi= 1.772453850905516027298167483341145182798e+0_rk   !pi^0.5
!
!----------------------------------------------------------------------
!     Machine dependent parameters  (DEFINED ABOVE)
!----------------------------------------------------------------------
!s    data xbig,xminin,eps/35.040e0,1.18e-38,1.19e-7/, & 
!s         xinf/3.4e38/
!
!d    data xbig,xminin,eps/171.624_rk,2.23d-308,2.22d-16/, & 
!d         xinf/1.79d308/
!----------------------------------------------------------------------
!     Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!----------------------------------------------------------------------
!S    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
!S   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
!S   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
!S   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
!S    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
!S   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
!S   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
!S   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!
      real(rk),dimension(8),parameter ::          &
         p = (/ -1.71618513886549492533811e+0_rk, &
                 2.47656508055759199108314e+1_rk, &
                -3.79804256470945635097577e+2_rk, &
                 6.29331155312818442661052e+2_rk, &
                 8.66966202790413211295064e+2_rk, &
                -3.14512729688483675254357e+4_rk, &
                -3.61444134186911729807069e+4_rk, &
                 6.64561438202405440627855e+4_rk /), &  
         q = (/ -3.08402300119738975254353e+1_rk, &
                 3.15350626979604161529144e+2_rk, &
                -1.01515636749021914166146e+3_rk, &
                -3.10777167157231109440444e+3_rk, &
                 2.25381184209801510330112e+4_rk, &
                 4.75584627752788110767815e+3_rk, &
                -1.34659959864969306392456e+5_rk, &
                -1.15132259675553483497211e+5_rk /)
!
!----------------------------------------------------------------------
!     Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
!s    data c/-1.910444077728e-03,8.4171387781295e-04,
!s   1     -5.952379913043012e-04,7.93650793500350248e-04,
!s   2     -2.777777777777681622553e-03,8.333333333333333331554247e-02,
!s   3      5.7083835261e-03/
!
      real(rk),dimension(7),parameter ::            &
         c = (/ -1.910444077728e-03_rk,             &
                 8.4171387781295e-04_rk,            &
                -5.952379913043012e-04_rk,          &
                 7.93650793500350248e-04_rk,        &
                -2.777777777777681622553e-03_rk,    &
                 8.333333333333333331554247e-02_rk, &
                 5.7083835261e-03_rk /)
!  
!----------------------------------------------------------------------
!     GO
!
      logipar = .false.
      fact = 1.0_rk 
      n = 0
      y = x
!      
      if (y .le. 0.0_rk ) then
!
!        Argument is negative
!
         y = -x
         y1 = aint(y,kind=rk)
         res = y - y1
         if (res .ne. 0.0_rk) then
            if (y1 .ne. aint(y1*0.5_rk,kind=rk)*2.0_rk) logipar = .true.
            fact = -picons / sin(picons*res)
            y = y + 1.0_rk
         else
            res = rinf
            goto 900
         endif
      endif
!
!     Argument is positive
!
      if ( y .lt. reps ) then
!
!        Argument .LT. EPS
!
         if ( y .ge. rtin ) then
            res = 1.0_rk / y
         else
            res = rinf
            goto 900
         endif
      else if (y .lt. 12.0_rk ) then
         y1 = y
         if (y .lt. 1.0_rk ) then
!
!           0.0 .LT. argument .LT. 1.0
!
            z = y
            y = y + 1.0_rk 
         else
!
!           1.0 .LT. argument .LT. 12.0, reduce argument if necessary
!
            n = int(y) - 1
            y = y - dble(n)   ! conv
            z = y - 1.0_rk 
         endif
!
!        Evaluate approximation for 1.0 .LT. argument .LT. 2.0
!
         xnum = 0.0_rk 
         xden = 1.0_rk
         do i = 1, 8
            xnum = (xnum + p(i)) * z
            xden = xden * z + q(i)
         enddo
         res = xnum / xden + 1.0_rk
         if (y1 .lt. y) then
!
!           Adjust result for case  0.0 .LT. argument .LT. 1.0
!
            res = res / y1
         else if (y1 .gt. y) then
!
!           Adjust result for case  2.0 .LT. argument .LT. 12.0
!
            do i = 1, n
               res = res * y
               y = y + 1.0_rk 
            enddo
         endif
      else
!
!        Evaluate for argument .GE. 12.0,
!
         if ( y .le. xbig ) then
            ysq = y * y
            sss = c(7)
            do i = 1, 6
               sss = sss / ysq + c(i)
            enddo
            sss = sss/y - y + sqrtpi
            sss = sss + (y-0.5_rk)*log(y)
            res = exp(sss)
         else
            res = rinf
            goto 900
         endif
      endif
!
!     Final adjustments and return
!
      if ( logipar ) res = -res
      if ( fact .ne. 1.0_rk ) res = fact / res
!      
  900 f = res
!  
      return
      end function specfun_gamma
!=======================================================================
!
!     Finding optimal parameters in a right-bounded region
!
!=====
!
      subroutine primsub_optimalparam_rb  (        & 
                 t, phi_s_star_j, phi_s_star_j1,   &
                 pj, qj, log_eps, epsmul10,        &
                 muj, hj, nj                       )
!
!     use mod_mlf_garrappa
!      
      implicit none 
!
!     Input & Output:
!
      real(rk),intent(in)     :: t, phi_s_star_j, phi_s_star_j1, pj, qj
      real(rk),intent(in)     :: log_eps, epsmul10
      real(rk),intent(out)    :: muj, hj
      integer(ik),intent(out) :: nj
!
!     Local variables: 
!        We shall copy value of the input log_eps to a local variable, 
!        say local_logeps. Since we shall correct this value internally.
!        The value of log_eps outside this scope is unchanged.
!
      real(rk) :: local_logeps
      real(rk) :: fac, f_max, f_min, f_bar, fq, fp, den, w 
      logical  :: conservative_error_analysis
      real(rk) :: sq_phi_star_j, threshold, sq_phi_star_j1
      real(rk) :: sq_phibar_star_j, sq_phibar_star_j1
      logical  :: adm_region
!
!     GO 
! 
      local_logeps = log_eps
 
      fac = 1.01_rk
 
      conservative_error_analysis = .true. 

!@Maximum value of fbar as the ration between tolerance and round-off unit

      f_max = exp(local_logeps - lnep)

!@Evaluation of the starting values for sq_phi_star_j and sq_phi_star_j1

      sq_phi_star_j  = sqrt(phi_s_star_j)
      threshold      = 2*sqrt((local_logeps - lnep)/t)
      sq_phi_star_j1 = min(sqrt(phi_s_star_j1), threshold-sq_phi_star_j)

!@Zero or negative values of pj and qj

      if ( pj < epsmul10 ) then
!
         if ( qj < epsmul10 ) then 
!
            sq_phibar_star_j  = sq_phi_star_j 
            sq_phibar_star_j1 = sq_phi_star_j1 
            adm_region = .true.
!
         else 
!
            sq_phibar_star_j = sq_phi_star_j 
!          
            if ( sq_phi_star_j > 0.0_rk ) then 
              f_min = fac*( sq_phi_star_j /                  &
                            (sq_phi_star_j1 - sq_phi_star_j) )**qj 
            else
               f_min = fac
            endif 
!          
            if ( f_min < f_max ) then 
               f_bar = f_min + (f_min/f_max)*(f_max-f_min) 
               fq = f_bar**(-1.0_rk/qj) 
               sq_phibar_star_j1 = ( 2 *sq_phi_star_j1 - &
                                     fq*sq_phi_star_j    ) / (2+fq)
               adm_region = .true.
            else
               adm_region = .false.
            endif 
!
         endif 
!
      else 
!
         if ( qj < epsmul10 ) then 
!
            sq_phibar_star_j1 = sq_phi_star_j1 
            f_min = fac*( sq_phi_star_j1 /               &
                         (sq_phi_star_j1-sq_phi_star_j) )**pj 
!
            if ( f_min < f_max ) then 
               f_bar = f_min + (f_min/f_max)*(f_max-f_min) 
               fp = f_bar**(-1.0_rk/pj) 
               sq_phibar_star_j = ( 2 *sq_phi_star_j + &
                                    fp*sq_phi_star_j1  )/(2-fp)
               adm_region = .true.
            else
               adm_region = .false.
            endif 
!
         else 
!
            f_min = fac * (sq_phi_star_j+sq_phi_star_j1) / ( &
                          (sq_phi_star_j1-sq_phi_star_j)**max(pj,qj) )
!
            if ( f_min < f_max ) then 
!            
               f_min = max(f_min,1.5_rk) 
               f_bar = f_min + (f_min/f_max)*(f_max-f_min) 
!               
               fp = f_bar**(-1.0_rk/pj) 
               fq = f_bar**(-1.0_rk/qj) 
!               
               if ( conservative_error_analysis ) then 
                  w = -phi_s_star_j1  * t / local_logeps 
               else
                  w =-2*phi_s_star_j1*t/(local_logeps-phi_s_star_j1*t)
               endif 
!               
               den = 2 + w - (1+w)*fp + fq 
!               
               sq_phibar_star_j = ( (2+w+fq)*sq_phi_star_j + &
                                    fp      *sq_phi_star_j1  )/den 
               sq_phibar_star_j1 = (-(1+w)*fq    *sq_phi_star_j + &
                                   (2+w-(1+w)*fp)*sq_phi_star_j1  )/den
!                                   
               adm_region = .true.
!               
            else
               adm_region = .false.
            endif
!
         endif 
!
      endif 
!
!
!
      if ( adm_region ) then 
!
         local_logeps = local_logeps  - log(f_bar)
!
         if ( conservative_error_analysis ) then 
            w = -sq_phibar_star_j1**2 * ( t/local_logeps )
         else
            w = -2*sq_phibar_star_j1**2 * t /(     &
               local_logeps - sq_phibar_star_j1**2 * t ) 
         endif 
!
         muj = ( ( (1+w)*sq_phibar_star_j  + &
                         sq_phibar_star_j1   )/(2+w) )**2
!
         hj =  -( pimul2/local_logeps ) * (                &
                   sq_phibar_star_j1 - sq_phibar_star_j ) / ( &
                   (1+w)*sq_phibar_star_j + sq_phibar_star_j1 ) 
!
         nj = ceiling( sqrt(1.0_rk-local_logeps/t/muj)/hj ) 
!
      else
!
         muj = 0.0_rk  
         hj = 0.0_rk
         nj = iinf
!
      endif 
!      
      return 
      end subroutine primsub_optimalparam_rb
!
!=======================================================================
!
!     Finding optimal parameters in a right-unbounded region
!
!=====
!
      subroutine primsub_optimalparam_ru  (  & 
                 t, phi_s_star_j,            &
                 pj, log_eps, epsmul10,      &
                 muj, hj, nj                 )
!                 
!     use mod_mlf_garrappa
!      
      implicit none 
!
!     Input & Output:
!
      real(rk),intent(in)     :: t, phi_s_star_j, pj
      real(rk),intent(in)     :: log_eps, epsmul10
      real(rk),intent(out)    :: muj, hj
      integer(ik),intent(out) :: nj
!
!     Local varaibles:
!  
      logical  :: istop 
      real(rk) :: f_min, f_max, f_tar, fbar, A, Q, w, u 
      real(rk) :: sq_muj, log_eps_phi_t, phi_t, threshold 
      real(rk) :: sq_phi_s_star_j, sq_phibar_star_j, phibar_star_j
!
!     GO 
!
!@Evaluation of the starting values for sq_phi_star_j
!
      sq_phi_s_star_j = sqrt(phi_s_star_j) 

      if ( phi_s_star_j > 0.0_rk ) then 
          phibar_star_j = phi_s_star_j*1.01_rk 
      else
          phibar_star_j = 1.0e-02_rk
      endif 

      sq_phibar_star_j = sqrt(phibar_star_j) 

!@Definition of some constants

      f_min = 1.0_rk 
      f_max = 10.0_rk 
      f_tar = 5.0_rk

!@Iterative process to look for fbar in [f_min,f_max]

      istop = .false. 
      do 

         phi_t = phibar_star_j * t 
         log_eps_phi_t = log_eps / phi_t 

         nj = ceiling((phi_t/picons)*( 1 - 1.5_rk*log_eps_phi_t + &
                                       sqrt(1-2*log_eps_phi_t) )  ) 
         A  = picons * nj / phi_t 

         sq_muj = sq_phibar_star_j * abs(4-A) / abs(7-sqrt(1+12*A))
         fbar   = ((sq_phibar_star_j-sq_phi_s_star_j)/sq_muj)**(-pj) 
         istop  = (pj < epsmul10) .or. (f_min < fbar .and. fbar < f_max)

         if ( istop ) then
            exit 
         else 
            sq_phibar_star_j = f_tar**(-1.0_rk/pj) * sq_muj + &
                               sq_phi_s_star_j
            phibar_star_j = sq_phibar_star_j**2 
         endif 

      enddo 

      muj = sq_muj**2 
      hj  = ( -3*A - 2 + 2*sqrt(1+12*A) )/(4-A)/real(nj,rk) 

!@Adjusting integration parameters to keep round-off errors under control
 
      threshold = (log_eps - lnep) / t 
!
      if ( muj > threshold ) then 
!
          if ( abs(pj) < epsmul10 ) then 
            Q = 0.0_rk 
          else 
            Q = f_tar**(-1.0_rk/pj) * sqrt(muj) 
          endif 
!
          phibar_star_j = ( Q + sqrt(phi_s_star_j) )**2 
!
          if ( phibar_star_j < threshold ) then 
!
              w = sqrt(lnep/(lnep-log_eps)) 
              u = sqrt(-phibar_star_j*t/lnep) 
!
              muj = threshold 
              nj = ceiling( w*log_eps/ pimul2 / (u*w-1.0_rk) )
              hj = sqrt(lnep/(lnep - log_eps))/real(nj,rk)
!
          else
!
              nj = iinf 
              hj = 0 
!              
          endif 
!
      endif 
!
      return 
      end subroutine primsub_optimalparam_ru
!=======================================================================
!
!     Derivative of Mittag-Leffer function with two parameters. 
!
!     Applying THEOREM 4.1 of the Gorenflo's paper. Make sure that afa>0
!
!=====
      subroutine sub_dermlf_multishot (        &
                 afa, bta, n, incz, z, ince, f )
!
!     use mod_mlf_garrappa
!      
      implicit none 
!
!     Input & Output: using explicit-shape array
!
      real(rk),intent(in)     :: afa, bta 
      integer(ik),intent(in)  :: n, incz, ince 
      complex(rk),dimension(incz,n),intent(in)  :: z
      complex(rk),dimension(ince,n),intent(out) :: f
!
!     Dependencies:
!
!     external :: sub_genmlf_multishot 
!     real(rk) :: specfun_gamma
!
!     Local varaibles:
!
      real(rk),parameter :: qconst = 0.99e0_rk
      complex(rk) :: zloc(4)
      real(rk)    :: d, w, absz, rho 
      integer     :: k, k0, k1, nz, j  
!
!     Check if the arguments afa>0. Otherwise, do nothing and report. 
!
!     if ( afa .le. 0.0_rk ) then 
!        write(012,'(a)') 'ERROR: Input to fmld1 wrong, afa<=0,',afa
!        return 
!     endif 
! 
      do j = 1,n
!
         zloc(1) = z(1,j) 
         absz    = abs(zloc(1))
! 
         if ( absz .eq. 0.0_rk ) then 
!
!           For |z|=0: use Eq. (38) only for k=0, the remainder 
!           as k>0 is zero
!
           zloc(2) = 1.0_rk / specfun_gamma( afa + bta )
 
         else if ( absz .le. qconst )  then 
!
!           Applying THEOREM 4.1 exactly:
!
!           For 0<|z|<=q, where we choose q=qconst,
!
!           +  Calculating k1 from (40) for cases of afa and D:
!
            if ( afa .le. 1.0_rk ) then 
!
!              For afa<=1: from (40) 
!
               d = afa*(afa - 4.0_rk*bta + 6.0_rk) + 1.0_rk
!
               if ( d .gt. 0.0_rk ) then 
                  w  = afa + bta - 1.5_rk
                  k1 = max( floor((3.0_rk - afa - bta)/afa) + 1,   &
                            floor((1.0_rk - 2.0_rk*w*afa + sqrt(d) &
                                  )/(2.0_rk*afa*afa)) + 1          )
               else
                  k1 = floor( (3.0_rk-afa-bta)/afa ) + 1
               endif 

            else 
!
!              For afa>1: from (40) 
!
               k1 = floor( (2.0_rk-afa-bta)/(afa-1.0_rk) ) + 1
 
            endif 
!
!           +  Calculating k0 from k1, for computing (39):
!
!Hint:      You may want to edit something right here to estimate 
!           the truncation errors by the relative error estimate, 
!           instead of the absolute error estimate. Check, plz!
!
            rho = present_epsilon
            k0  = max( k1, floor(log(rho*(1.0_rk-absz))/log(absz)) )
!
!           +  Calculating E' from (39): summing up directly for k=0,k0
!
            zloc(2) = 1/ specfun_gamma( afa + bta )
!            
            do k = 1,k0
               zloc(2) = zloc(2) + (k+1)*zloc(1)**k /         &
                                   specfun_gamma(afa+bta+afa*k)
            enddo
! 
         else 
!
!           For |z|>q, where q=0.1, use (43) with Mittag-Leffler 
!           function
!
            call sub_genmlf_multishot (                             &
                 afa, bta-1.0_rk, 1.0_rk, 1, 1, zloc(1), 1, zloc(3) )
! 
            call sub_genmlf_multishot (                             &
                 afa, bta,        1.0_rk, 1, 1, zloc(1), 1, zloc(4) )
! 
            zloc(2) = ( zloc(3) - zloc(4)*(bta-1.0_rk) )/(afa*zloc(1))
!
         endif 
! 
         f(1,j) =  zloc(2)
! 
      enddo
      return 
      end subroutine 
!=======================================================================
      end module mod_mlf_garrappa
!======================================================================

!     __END__
