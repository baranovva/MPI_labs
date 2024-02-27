      PROGRAM OPENMP_GZ
      INCLUDE 'omp_lib.h'
      INTEGER, PARAMETER :: N = 101
      INTEGER I,J,L,MAXIT,T,S
      REAL(8), ALLOCATABLE, DIMENSION(:,:):: W,X,Y
      REAL(8) ERROR, EPS, T1, T2, dX, dY, C
      MAXIT = 100000 !Максимальное колво итераций
      EPS = 1E-6 !Точность
      ERROR = 1E6 !Ошибка
      L = 1
      ALLOCATE (W(N,N),X(N,N),Y(N,N))
      
      OPEN(1, FILE='101.msh')!чтение сетки с файла
      DO J=1,N
          DO I=1,N
              READ(1,*) X(I,J), Y(I,J)
          ENDDO
      ENDDO
      CLOSE(1)
          
      W(2:N-1,2:N-1) = 0.0! задание граничных условий
      DO J=1,N
         W(1,J) = EXP(-1.5*(X(1,J))**2 + 2*(Y(1,J))**2)
         W(N,J) = EXP(-1.5*(X(N,J))**2 + 2*(Y(N,J))**2)
      ENDDO
      DO I=1,N
         W(I,1) = EXP(-1.5*(X(I,1))**2 + 2*(Y(I,1))**2)
         W(I,N) = EXP(-1.5*(X(I,N))**2 + 2*(Y(I,N))**2)
      ENDDO

      T1 = OMP_GET_WTIME()
      DO S=1,100000000
      DO WHILE (L .LE. MAXIT .AND. ERROR .GT. EPS)
        ERROR = 0.0
!$OMP PARALLEL PRIVATE(I,J,dX,dY,C)!Директива PARALLEL DO позволяет указать,
какие циклы компилятор должен распараллелить.
!$OMP DO REDUCTION(MAX:ERROR)  ! цикл DO,
C помощью  REDUCTION находится максимальная ошибка
        DO J = 2, N-1
          DO I = 2, N-1
            dY = Y(I,J+1) - Y(I,J)
            dX = X(I+1,J) - X(I,J)
            C = W(I,J)
            W(I,J) = 0.5 * (dY*(W(I-1,J) + W(I+1,J)) + dX*(W(I,J-1) + W(I,J+1)) - 0.2*dX*dY)/(dX + dY)
            ERROR = MAX(ERROR, ABS(W(I,J) - C))
            
          END DO
        ENDDO
!$OMP END DO!Директива END PARALLEL DO позволяет указать конец цикла DO,
указанного в директиве PARALLEL DO.
!$OMP END PARALLEL
        L = L + 1
      END DO
      END DO
      T2 = OMP_GET_WTIME()
      PRINT *, 'computational time: ', T2-T1, 's'
      PRINT *, L
      
      !OPEN(2,FILE='Res.dat')  !вывод данных в файл
      !DO J=1,N
      !  DO I=1,N
      !      WRITE(2,*) X(I,J), Y(I,J), W(I,J)
      !  END DO
      !END DO
      !CLOSE(2) 
      
      DEALLOCATE (W,X,Y)
      END PROGRAM

