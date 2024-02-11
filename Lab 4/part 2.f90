PROGRAM MPI
INCLUDE 'mpif.h'
PARAMETER (N=180)
INTEGER :: rank, SIZE, ierr, temp_1, temp_2, temp_3,  next, prev, STATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION:: A_temp(N), B_temp(N,2), C(N,N), time
DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), C_temp(:)

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ierr)

IF (rank == 0) THEN
  ALLOCATE(A(N,N),B(N,N))
  DO I=1,N
     DO J=1,N
        A(J,I)=1.0
        B(J,I)=1.0
     ENDDO
  ENDDO
  time = MPI_WTIME()
ENDIF

ALLOCATE(C_temp(SIZE))
temp_1 = 1
DO I = rank+1, N, size
  K = MIN(temp_1+size-1,N)
  CALL MPI_SCATTER(A(:,temp_1:K),N,MPI_DOUBLE_PRECISION,A_temp,N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  temp_2 = 1
  DO I1= rank+1, N, size
     DO I2 = 1,N
CALL MPI_SCATTER(B(temp_2:MIN(temp_2+size-1,N),I2),1,MPI_DOUBLE_PRECISION,B_temp(I2,1),&
1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     ENDDO

     prev = rank - 1
     next = rank + 1
     IF (rank == 0) prev = SIZE-1
     IF (rank == SIZE-1) next = 0
     
	 temp_3 = rank+1
	 DO I2 = 1, SIZE
	   C_temp(temp_3) = DOT_PRODUCT(A_temp(:),B_temp(:,1))
       CALL MPI_SEND(B_temp(:,1),N,MPI_DOUBLE_PRECISION,next,1,MPI_COMM_WORLD,ierr)
       CALL MPI_RECV(B_temp(:,2),N,MPI_DOUBLE_PRECISION,prev,1,MPI_COMM_WORLD,STATUS,ierr) 
	   B_temp(:,1) = B_temp(:,2)
    
	   IF (temp_3 == 1) THEN 
	   	temp_3 = SIZE
	   ELSE 
	   	temp_3 = temp_3-1
	   END IF
	   
     ENDDO
CALL MPI_GATHER(C_temp,SIZE,MPI_DOUBLE_PRECISION,C(temp_2:MIN(temp_2+size-1,N),temp_1:K),&
SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
	 temp_2=temp_2+SIZE
  ENDDO
  temp_1=temp_1+SIZE
ENDDO

IF (rank == 0) THEN
	PRINT *, 'Time:', MPI_WTIME() - time
	!DO I=1,N
	!	PRINT*,C(:,I)
	!ENDDO
ENDIF

CALL MPI_FINALIZE(ierr)
END PROGRAM 
