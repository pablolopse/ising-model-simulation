      ! ===========================================
      ! 
      !  SIM. MODELO DE ISING 2D - METROPOLIS
      !
      ! ===========================================
      PROGRAM METROPOLIS_2D
      IMPLICIT NONE

      ! DEFINIR LAS VARIABLES DEL PROGRAMA
      INTEGER SIZE,DIM,N,JMAG,STEPS,i,j,k
      PARAMETER(SIZE=100, DIM=2, N=SIZE*SIZE, JMAG=1,STEPS=20000)
      INTEGER matrix(SIZE, SIZE) ! NO ADAPTABLE SI SE CAMBIA DIM
      ! MATRIZ DE VECINOS PARA PRECALCULAR CADA POSICION VECINA
      REAL*8 exponentials(0:DIM), E, M,energies(STEPS),mags(STEPS),
     &       avg_e,avg_m,avg_e2
      REAL itemp, dtemp, ftemp

      itemp = 0.2
      dtemp = 0.2
      ftemp = 5.0


      ! Inicializar la matriz una vez, luego bajar en temperatura
      call INITIALIZE(matrix, SIZE, E, M, JMAG)

      ! Abrir o crear el fichero para guardar los resultados
      open(10,File="metropolis_2d.txt")

      ! Bucle de temperatura de alto a bajo
      do i=int((ftemp-itemp)/dtemp),0,-1
        print*, "Performing simulation at T ",itemp+i*dtemp
        ! Calcular las exponenciales para no tener que hacer
        ! operaciones en coma flotante para cada temperatura.
        do k = 0, DIM
          exponentials(k) = exp(-2.0*JMAG*2.0*k/(itemp+i*dtemp))
        end do
        
        ! Inicializar arrays de energia y magnetizacion a cero
        energies = 0
        mags = 0
        
        ! Realizar STEPS pasos Monte Carlo
        do j=1, STEPS
          call MONTECARLO_STEP(matrix,SIZE,E,M,N,DIM,
     &                         exponentials)
          ! Guardar el valor de E y M de cada paso
          energies(j) = E
          mags(j) = M
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        do j=6001,STEPS
          avg_e = avg_e + energies(j)
          avg_m = avg_m + abs(mags(j))
          avg_e2 = avg_e2 + (energies(j)**2)
        end do
        write(10,*) itemp+i*dtemp,avg_e/14000/N,avg_m/14000/N,
     &              (avg_e2/14000 - (avg_e/14000)**2)
     &              /(N*dble(itemp+i*dtemp)**2)
      end do
      close(10)

      END PROGRAM


      SUBROUTINE INITIALIZE(matrix, SIZE,E,M,JMAG)
      IMPLICIT NONE
      INTEGER i,j,SIZE,matrix(SIZE,SIZE),JMAG
      REAL*8 sum,E,M

      ! Inicializar la matriz con condiciones T=inf, orientacion aleatoria
      do i=1, SIZE
        do j=1, SIZE
          call random_number(sum)
          matrix(i,j)=2*nint(sum)-1
        end do
      end do

      ! Calcular magnetizacion y energia iniciales
      E = 0
      M = 0
      do i=1, SIZE ! Iterar sobre cada elemento de la matriz
        do j=1, SIZE
          if (i+1.gt.SIZE) then
            sum = matrix(1,j)
          else
            sum = matrix(i+1,j)
          end if
          if (j+1.gt.SIZE) then
            sum = sum + matrix(i,1)
          else
            sum = sum + matrix(i,j+1)
          end if
          
          ! Sumar cada elemento a M y E
          M = M + matrix(i,j)
          E = E - JMAG*matrix(i,j)*sum
        end do
      end do

      END SUBROUTINE


      SUBROUTINE MONTECARLO_STEP(matrix, SIZE, E, M, N, 
     &            DIM,exponentials)
      IMPLICIT NONE
      INTEGER SIZE,DIM,N,i,j,k,matrix(SIZE, SIZE)
      REAL*8 E,M,dE,sum,exponentials(0:DIM),randv
      
      do k=1, N
        ! Seleccionamos una posicion aleatoria en la matriz
        call random_number(randv)
        j = int(randv*SIZE) + 1
        call random_number(randv)
        i = int(randv*SIZE) + 1

        sum = 0.0
        ! Sumar vecinos izquierda y derecha
        if (i+1.gt.SIZE) then
          sum = sum + matrix(1,j)
          sum = sum + matrix(i-1,j)
        else if (i-1.eq.0) then
          sum = sum + matrix(SIZE,j)
          sum = sum + matrix(i+1,j)
        else
          sum = sum + matrix(i+1,j)
          sum = sum + matrix(i-1,j)
        end if
        ! Ahora para arriba y abajo
        if (j+1.gt.SIZE) then
          sum = sum + matrix(i,1)
          sum = sum + matrix(i,j-1)
        else if (j-1.eq.0) then
          sum = sum + matrix(i,SIZE)
          sum = sum + matrix(i,j+1)
        else
          sum = sum + matrix(i,j+1)
          sum = sum + matrix(i,j-1)
        end if

        ! Calcular el cambio en la energia
        dE = 2*matrix(i,j)*sum

        ! Aceptar o rechazar la transicion
        if (dE.le.0) then
          matrix(i,j) = -matrix(i,j)
          E = E + dE
          M = M + matrix(i,j)*2
        else
          call random_number(randv)
          if (randv.lt.exponentials(int(dE/(2)/2))) then
            matrix(i,j) = -matrix(i,j)
            E = E + dE
            M = M + matrix(i,j)*2
          end if
        end if
      end do
      
      END SUBROUTINE