program phantom

   use gif_module

   implicit none

! Input parameters that define the problem:
   character                   :: arg*80, value*80

   integer                     :: n, i, k
   real(kind=8)                :: stddev, err, lambda
   real(kind=8), allocatable   :: f_phantom(:,:), f_noisy(:,:), f_filtered(:,:)

! For timing:
   integer                     :: t0, t1, tb, te, clock_rate, clock_max

   write(*,*) 
   call system_clock ( t0, clock_rate, clock_max )

! Defaults:
   n = 512
   stddev = 0.
   lambda = 0.

! Read parameters from command line:
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)

      select case(arg)
      case('-npx')
         call get_command_argument(i+1, value)
         read(value,'(i6)') n
      case('-sigma')
         call get_command_argument(i+1, value)
         read(value,'(e9.4)') stddev
      case('-lambda')
         call get_command_argument(i+1, value)
         read(value,'(e9.4)') lambda
      end select
   end do

   allocate( f_phantom(n,n), f_noisy(n,n), f_filtered(n,n) )

! Generate phantom image
   call shepp_logan( f_phantom, n )

! Timing for adding noise and filtering (part of the program that will be parallelised)

! Add noise
   f_noisy = imnoise( f_phantom, stddev )
   err = sqrt( sum( (f_noisy - f_phantom)**2 ) )/n
   write(*,'(a,e8.2)') 'Noisy image, Sigma = ',err

   call system_clock ( tb, clock_rate, clock_max )

! Filter the image
   f_filtered = cg( f_noisy, lambda )
   err = sqrt( sum( (f_filtered - f_phantom)**2 ) )/n
   write(*,'(a,e8.2)') 'Filtered image, Sigma = ',err

   call system_clock ( te, clock_rate, clock_max )

! Plot the results
   call gif_images( f_phantom, f_noisy, f_filtered )

! Report timings
   call system_clock ( t1, clock_rate, clock_max )
   write(*,'(a,e8.2,a)') 'Elapsed time for filtering = ', real(te-tb)/real(clock_rate), 's.'
   write(*,'(a,e8.2,a)') 'Total elapsed time     = ', real(t1-t0)/real(clock_rate), 's.'
   write(*,'(a)') '_________________________________________________________________________' 

contains

function cg( f_noisy, lambda ) result(f_filtered)
!
! CG algorithm 
!
   implicit none

   real(kind=8), intent(in)    :: lambda
   real(kind=8), intent(in)    :: f_noisy(:,:)
   real(kind=8)                :: f_filtered(size(f_noisy,1),size(f_noisy,2))

   real(kind=8)                :: b(size(f_noisy,1),size(f_noisy,2))
   real(kind=8)                :: x(size(f_noisy,1),size(f_noisy,2))
   real(kind=8)                :: r(size(f_noisy,1),size(f_noisy,2))
   real(kind=8)                :: p(size(f_noisy,1),size(f_noisy,2))
   real(kind=8)                :: Ap(size(f_noisy,1),size(f_noisy,2))
   real(kind=8)                :: alpha, beta, gamma, xi
   real(kind=8)                :: tol, normr, normb
   integer                     :: it, maxit

   if ( lambda == 0. ) then
      f_filtered = f_noisy
      return
   end if

   maxit = 100
   tol   = 1.e-8
   b = lambda*f_noisy
   x = f_noisy
   r = b - mv( x, lambda )
   p = r
   Ap = mv( p, lambda )

   gamma = inprod(r,r)
   normr = sqrt( gamma )
   normb = sqrt( inprod(b,b) )

   it = 0
   do while ( normr/normb > tol .and. it < maxit )

      xi = inprod(p,Ap)
      alpha = gamma/xi

      x = x + alpha*p
      r = r - alpha*Ap;

      beta  = 1/gamma
      gamma = inprod(r,r)
      normr = sqrt(gamma)
      beta  = gamma*beta
      p     = r + beta*p
      Ap    = mv( p, lambda )

      it = it+1
   end do

   f_filtered = x
   write(*,'(a,i4,1x,a,e9.3)') 'CG terminated after ',it, &
           'iterations with relative residual norm ', normr/normb

end function cg

function inprod( v, w ) result(vw)
!
! Inner product. Note that we have stored the vectors in a two-dimensional
! array, so we can not use the intrinsic dot_product
! 
   implicit none

   real(kind=8), intent(in)    :: v(:,:), w(:,:)
   real(kind=8)                :: vw
   integer                     :: nx, ny
   integer                     :: i, j
!
   nx = size(v,1)
   ny = size(v,2)

   vw = 0.
   do j = 1, ny
      do i = 1, nx
         vw = vw + v(i,j)*w(i,j)
      end do
   end do
end function inprod

function mv( v, lambda ) result(w)
!
! Matrix vector multiplication: w = lambda v + A v
! in which A represents -Laplacian with Neumann boundary conditions
! 
   implicit none

   real(kind=8), intent(in)    :: lambda
   real(kind=8), intent(in)    :: v(:,:)
   real(kind=8)                :: w(size(v,1),size(v,2))
   real(kind=8)                :: v_halo(0:size(v,1)+1,0:size(v,2)+1)
   integer                     :: nx, ny
   integer                     :: i, j

   nx = size(v,1)
   ny = size(v,2)
! Add extra row/column around the domain to account for boundary conditions
! or data on other processor. Note that these points are outside the domain.
! They can be eliminated by the homogeneous Neumann condition. This condition
! says that the points on a row/column outside the boundary are equal to the 
! neighbouring points on the row/column inside the boundary.
   v_halo = 0.
   v_halo(1:nx,1:ny) = v                   ! Points in the interior
   v_halo(0,1:ny)    = v_halo(1,1:ny)      ! Neumann condition on the left
   v_halo(nx+1,1:ny) = v_halo(nx,1:ny)     ! Neumann condition on the right
   v_halo(1:nx,0)    = v_halo(1:nx,1)      ! Neumann condition on the bottom
   v_halo(1:nx,ny+1) = v_halo(1:nx,ny)     ! Neumann condition on the top

! Now we can just perform the stencil operation:
   w = lambda*v
   do j = 1, ny
      do i = 1, nx
         w(i,j) = w(i,j) &
             + 4.*v_halo(i,j) -v_halo(i-1,j) -v_halo(i+1,j) -v_halo(i,j-1) -v_halo(i,j+1)
      end do
   end do

end function mv

function imnoise( f_phantom, stddev ) result(f_noisy)
!
! Add noise to the image
!

   implicit none

   real(kind=8), intent(in)  :: f_phantom(:,:)
   real(kind=8)              :: f_noisy(size(f_phantom,1),size(f_phantom,2))
   real(kind=8)              :: fn(size(f_phantom,1),size(f_phantom,2)), stddev
   real(kind=8)              :: stddev_rand
   integer                   :: nx, ny

! Add noise to image
   nx = size(f_phantom,1)
   ny = size(f_phantom,2)

   call RANDOM_SEED
   call RANDOM_NUMBER(fn)
   fn = (fn - 0.5)
   stddev_rand = sqrt(sum( fn**2 )/(nx*ny))
   f_noisy = f_phantom + (stddev/stddev_rand)*fn

end function imnoise

subroutine gif_images( model_image, noisy_image, filtered_image )
!
! Create gif-images and write to file
!
   implicit none

   real(kind=8), allocatable     :: model_image(:,:), noisy_image(:,:), &
                                    filtered_image(:,:)
   integer, parameter            :: n_colours = 256
   integer, allocatable          :: gif_image(:,:,:)
   integer                       :: map(1:3,0:n_colours-1)
   character(len=23)             :: gif_name
   integer                       :: i,j,nx,ny

   nx = size(noisy_image,1)
   ny = size(noisy_image,2)
!
! Colour map grey
   do i = 0,n_colours-1
      map(:,i ) = i
   end do

! Model image
   allocate(gif_image(1,nx,ny))
   gif_name = 'phantom.gif'
   do j = 1, ny
      do i = 1, nx
         gif_image(1,i,j) = int(model_image(i,n-j+1)*n_colours)
      end do
   end do
   where ( gif_image > n_colours-1 ) gif_image = n_colours-1
   where ( gif_image < 0 ) gif_image = 0
   call write_animated_gif( gif_name, gif_image, map )

! Noisy image
   gif_name = 'noisy.gif'
   do j = 1, ny
      do i = 1, nx
         gif_image(1,i,j) = int(noisy_image(i,n-j+1)*n_colours)
      end do
   end do
   where ( gif_image > n_colours-1 ) gif_image = n_colours-1
   where ( gif_image < 0 ) gif_image = 0
   call write_animated_gif( gif_name, gif_image, map )

! Filtered image
   gif_name = 'filtered.gif'
   do j = 1, ny
      do i = 1, nx
         gif_image(1,i,j) = int(filtered_image(i,n-j+1)*n_colours)
      end do
   end do
   where ( gif_image > n_colours-1 ) gif_image = n_colours-1
   where ( gif_image < 0 ) gif_image = 0
   call write_animated_gif( gif_name, gif_image, map )

end subroutine gif_images
   
subroutine shepp_logan( f, n )

!
!    Larry Shepp, Ben Logan,
!    The Fourier reconstruction of a head section,
!    IEEE Transactions on Nuclear Science,
!    Volume  NS-21, June 1974, pages 21-43.
!
!  Parameters:
!
!    Input, integer n, the number of points in each direction.
!
!    Output, real f(n,n), the image values.
!
!  Local parameters:
!
!    Local, integer CHOICE:
!    1, use Archibald's (and Shepp and Logan's) level values;
!    2, use Matlab's level values;
!    3, use Matlab's enhanced contrast level values.
!

   implicit none

   integer ( kind = 4 ) i,j,n
   real ( kind = 8 ) c(4)
   real ( kind = 8 ), dimension ( 4 ) :: c1 = (/ 2.0E+00, -0.98E+00, -0.02E+00, +0.01E+00 /) 
   real ( kind = 8 ), dimension ( 4 ) :: c2 = (/ 1.0E+00, -0.98E+00, -0.02E+00, +0.01E+00 /) 
   real ( kind = 8 ), dimension ( 4 ) :: c3 = (/ 1.0E+00, -0.8E+00,  -0.2E+00,  +0.1E+00  /) 
   integer ( kind = 4 ) choice
   real ( kind = 8 ) eta1
   real ( kind = 8 ) eta2
   real ( kind = 8 ) f(n,n)
   real ( kind = 8 ), parameter :: pi = 3.141593E+00
   real ( kind = 8 ) xi1
   real ( kind = 8 ) xi2
   real ( kind = 8 ) x, x_max, x_min, dx
   real ( kind = 8 ) y, y_max, y_min, dy

   x_min = -1.0E+00
   x_max = +1.0E+00
   dx = (x_max-x_min)/n
   y_min = -1.0E+00
   y_max = +1.0E+00
   dy = (y_max-y_min)/n
 
   choice = 3

   if ( choice == 1 ) then
     c = c1
   else if ( choice == 2 ) then
     c = c2
   else
     c = c3
   end if


   y = y_min + dy/2
   do j = 1, n
      x = x_min + dx/2
      do i = 1, n

         if ( ( x / 0.69E+00 )**2 + ( y / 0.92E+00 )**2 <= 1.0E+00 ) then
            f(i,j) = f(i,j) + c(1)
         end if
 
         if ( ( x / 0.6624E+00 )**2 + ( ( y + 0.0184E+00 ) / 0.874E+00 )**2 <= 1.0E+00 ) then
            f(i,j) = f(i,j) + c(2)
         end if

         xi1  =   ( x - 0.22E+00 ) * cos ( 0.4E+00 * pi ) &
                  + y              * sin ( 0.4E+00 * pi )
         eta1 = - ( x - 0.22E+00 ) * sin ( 0.4E+00 * pi ) &
                  + y              * cos ( 0.4E+00 * pi )

         xi2  =   ( x + 0.22E+00 ) * cos ( 0.6E+00 * pi ) &
                  + y              * sin ( 0.6E+00 * pi )
         eta2 = - ( x + 0.22E+00 ) * sin ( 0.6E+00 * pi ) &
                  + y              * cos ( 0.6E+00 * pi )
 
         if ( ( xi1 / 0.31E+00 )**2 + ( eta1 / 0.11E+00 )**2 <= 1.0E+00 .or. &
              ( xi2 / 0.41E+00 )**2 + ( eta2 / 0.16E+00 )**2 <= 1.0E+00 ) then
            f(i,j) = f(i,j) + c(3)
         end if

         if ( (( x              /0.21E+00  )**2 + (( y - 0.35E+00 )/0.25E+00  )**2 <= 1.0E+00 ) .or. &
              (( x              /0.046E+00 )**2 + (( y - 0.10E+00 )/0.046E+00 )**2 <= 1.0E+00 ) .or. &
              (( x              /0.046E+00 )**2 + (( y + 0.1E0+00 )/0.046E+00 )**2 <= 1.0E+00 ) .or. &
              (((x + 0.080E+00) /0.046E+00 )**2 + (( y + 0.605E+00)/0.023E+00 )**2 <= 1.0E+00 ) .or. &
              (( x              /0.023E+00 )**2 + (( y + 0.605E+00)/0.023E+00 )**2 <= 1.0E+00 ) .or. &
              (((x - 0.06E+00 ) /0.023E+00 )**2 + (( y + 0.605E+00)/0.023E+00 )**2 <= 1.0E+00 ) ) then
            f(i,j) = f(i,j) + c(4)
         end if
         x = x + dx
      end do
      y = y + dy
   end do

   return
end subroutine shepp_logan 

end program phantom
