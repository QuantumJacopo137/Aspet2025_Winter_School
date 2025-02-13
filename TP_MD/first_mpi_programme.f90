program first_prog_mpi

implicit none

include 'mpif.h'

integer           :: ierr
integer           :: nprocs         ! total number of processes
integer           :: ME             ! CPU index
integer,parameter :: source=0
character(len=1)  :: numero

call MPI_Init(ierr)

call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr) ! Nombre of procs dispo

call MPI_Comm_rank(MPI_COMM_WORLD,ME,ierr)     ! Id de chaque proc

print*, "nprocs", nprocs, "ME", ME

write(numero,'(I1)') ME

if(ME==0) then
print*, 'I am', ME
endif

open(6,file='TEST'//numero)

write(6,'(A,I2)'), "Hello from processe number ", ME

close(6)

!call MPI_Bcast(pi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!call MPI_Reduce(pi,bob,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr);

!call MPI_AllReduce(pi,bob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr);

!call MPI_Finalize(ierr)

end program
