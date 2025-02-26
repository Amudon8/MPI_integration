program mpi_integration
implicit none
include 'mpif.h'

double precision, parameter :: pi = 3.14159265359d0
double precision, allocatable :: x(:), local_a(:), global_a(:)
double precision :: local_area, global_area, h
integer :: n, local_n, i
integer :: rank, comm, comm_size, ierr



! Initialize MPI
call MPI_Init(ierr)
comm = MPI_COMM_WORLD
call MPI_Comm_rank(comm, rank, ierr)
call MPI_Comm_size(comm, comm_size, ierr)

n = 4096
local_n = n/comm_size
h = pi/dble(n-1)
allocate(x(local_n))
allocate(local_a(local_n))
allocate(global_a(n))

do i = 1, local_n
   x(i) = (i-1)*h + rank*local_n*h
   local_a(i) = sin(x(i))
enddo
call MPI_GATHER(local_a, local_n, MPI_DOUBLE_PRECISION, global_a, local_n, MPI_DOUBLE_PRECISION, 0, comm, ierr)

call cross_node_integration(comm, ierr, local_n, local_a, h, global_area)
if (rank == 0) then
   print *,'approx. area = ', global_area
endif


local_area = 0.0d0
global_area = 0.0d0
call integration_1d(local_a, local_n, h, local_area)
call MPI_REDUCE(local_area, global_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
if (rank == 0) then
   print *,'Reduced area = ', global_area

   call integration_1d(global_a, n, h, global_area)
   print *,'area = ', global_area
endif

call MPI_Finalize(ierr)




contains
subroutine integration_1d(f, n, dy, area)
    implicit none
    integer :: i
    integer, intent(in) :: n
    double precision, intent(in) :: dy, f(n)
    double precision, intent(out) :: area
    area = 0.0d0
    do i = 1, n-1
        area = area + 0.5 * dy * (f(i) + f(i+1))
    end do
end subroutine integration_1d

subroutine cross_node_integration(comm, ierr, local_n, local_a, h, global_area)
   implicit none
   integer, intent(in) :: local_n, comm, ierr
   double precision, intent(in) :: h, local_a(local_n)
   double precision, intent(out) :: global_area
   integer :: i, rank

   call MPI_Comm_rank(comm, rank, ierr)
   do i = 1, local_n
      if (rank == 0 .and. i == 1) then
         local_area = local_area + local_a(i)
      elseif (rank == comm_size - 1 .and. i == local_n) then
         local_area = local_area + local_a(i)
      else
         local_area = local_area + 2.0d0*local_a(i)
      endif
   enddo
   local_area = 0.5d0 * h * local_area
   call MPI_REDUCE(local_area, global_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
end subroutine    
end program mpi_integration

