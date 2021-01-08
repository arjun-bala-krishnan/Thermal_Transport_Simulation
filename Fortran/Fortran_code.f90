program rw_test

  use functions_module

  implicit none

  integer, parameter  ::  nt = 500, np_t = 100, n_iter = 1 !ok
  real (kind = 8), parameter  ::  sig_therm = 0.5d0, dx = 0.1 !ok
  real (kind = 8), dimension(2), parameter  ::  box = [10.d0,10.d0], box2 = box/2.d0 !ok
  real (kind = 8), dimension(2,np_t*2*nt)   ::  x !ok
  real (kind = 8), dimension(nt,-nint(box(1)/dx)+1:nint(box(1)/dx)+1)   ::  his 
  real (kind = 8), dimension(2*np_t*nt)   ::  u
  real (kind=8), dimension(2)  ::  x_avg
  integer  ::  i,k,ix,iy,it,i_iter,n_part,k1,k2

	
  n_part = np_t
  x_avg = 0.d0

! np_t is the number of aprticles, injected every timestep PER species

  x(1,1:nt*np_t)           = -box2(1)
  x(1,nt*np_t+1:nt*2*np_t) =  box2(1)
  x(2,:)                   =  -box2(2) + fib_rnd()*box(2)
  
  u(1:nt*np_t)             =  1.d0
  u(nt*np_t+1:nt*2*np_t)   = -1.d0

  i_iter = 0
  his = 0.d0

  k1=0
  k2 = 0
  
!   iteration over independent experiments
  do while(i_iter < n_iter)

     n_part = 2*np_t
     x(1,1:nt*np_t)           = -box2(1)
     x(1,nt*np_t+1:nt*2*np_t) =  box2(1)
     x(2,:)                   = -box2(2) + fib_rnd() *box(2)
     
     ! iteration over timesteps (corresponding to 1 experiment)
     do it = 1, nt

        print*,it,i_iter
        k=0
        ! loop over particles with positive energies
        do i = 1, n_part/2
           x(1,i) = x(1,i) + sig_therm * gasdev()
           x(2,i) = x(2,i) + sig_therm * gasdev()

           ! reflective boundary conditions
           if( x(1,i) < -box(1)/2 )then
              x(1,i) = -box(1) - x(1,i)
           else if( x(1,i) > box(1)/2 )then
              x(1,i) = box(1) - x(1,i)
           endif
    
           ! perodic boundary conditions
           if( x(2,i)  > box(2)/2 )then
	      x(2,i) = x(2,i) - box(2)
           else if( x(2,i) < -box(2)/2 )then
	      x(2,i) = x(2,i) + box(2)
           endif
           
           k1=k1+1
	   if(i> n_part/2 - np_t)then
	      x_avg(1) = (real(k1-1)*x_avg(1)+x(1,i)+box2(1))/real(k1)
           endif
           ix = nint(x(1,i) / dx)
!print*,"first part",ix
	   his(it,ix) = his(it,ix) + u(i)
        enddo
	

        ! loop over particles with negtive energies
        do i = nt*np_t+1, nt*np_t+n_part/2
           x(1,i) = x(1,i) + sig_therm * gasdev()
           x(2,i) = x(2,i) + sig_therm * gasdev()

           ! reflective boun dary conditions
           if( x(1,i) < -box(1)/2 )then
              x(1,i) = -box(1) - x(1,i)
           else if( x(1,i) > box(1)/2 )then
              x(1,i) = box(1) - x(1,i)
           endif
    
           ! perodic boundary conditions
           if( x(2,i)  > box(2)/2 )then
	      x(2,i) = x(2,i) - box(2)
           else if( x(2,i) < -box(2)/2 )then
	      x(2,i) = x(2,i) + box(2)
           endif

           k2=k2+1
	   if(i> nt*np_t+n_part/2 - np_t)then
	      x_avg(2) = (real(k2-1)*x_avg(2)+box2(1)-x(1,i))/real(k2)
           endif

           ix = nint(x(1,i) / dx)
!print*,"2nd part",ix
           his(it,ix) = his(it,ix) + u(i)
        enddo

        n_part = n_part + 2*np_t
     enddo

     i_iter = i_iter+1
  enddo

    his = his / real(n_iter)
  
  print*," average initial displacement: ",x_avg

  open(unit=99,file="rw_energy_histogram.dat")
  
  do it = 1, nt
     do i = -nint(box(1)/dx)+1,nint(box(1)/dx)+1
        write(99,"(3es14.5)") real(it), real(i)*dx, his(it,i)
     enddo
     write(99,*)
  enddo
  close(99)
	
end program rw_test
