      ! ===========================================
      ! 
      !    SIM. MODELO DE ISING 2D - WOLFF
      !
      ! ===========================================
      PROGRAM WOLFF_2D
      IMPLICIT NONE

      ! DEFINIR LAS VARIABLES DEL PROGRAMA
      INTEGER SIZE,DIM,N,JMAG,STEPS,NSKIP,NAVG,i,j
      PARAMETER(SIZE=100, DIM=2, N=SIZE*SIZE, JMAG=1,STEPS=20000)
      PARAMETER(NSKIP=6000, NAVG=STEPS-NSKIP)
      INTEGER smatrix(SIZE, SIZE)
      REAL*8 E, M, energies(STEPS), mags(STEPS),
     &       avg_e, avg_m, avg_e2, padd
      REAL itemp, dtemp, ftemp
      INTEGER csize, flipped

      itemp = 0.2
      dtemp = 0.2
      ftemp = 5.0

      ! Inicializar la matriz una vez, luego bajar en temperatura
      call INITIALIZE(smatrix, SIZE, E, M, JMAG)

      ! Abrir o crear el fichero para guardar los resultados
      open(10,File="wolff_2d.txt")

      ! Bucle de temperatura de alto a bajo
      do i=int((ftemp-itemp)/dtemp),0,-1
        print*, "Performing simulation at T ",itemp+i*dtemp
        padd = 1.0d0 - exp(-2.0d0/(itemp+i*dtemp))
        
        ! Inicializar arrays de energia y magnetizacion a cero
        energies = 0
        mags = 0
        
        ! Realizar STEPS barridos Monte Carlo
        ! Cada barrido = suficientes volteos de cluster para voltear ~N espines
        do j = 1, STEPS
          flipped = 0
          do while (flipped .lt. N)
            call MONTECARLO_STEP(smatrix,SIZE,E,M,N,padd,csize)
            flipped = flipped + csize
          end do
          energies(j) = E
          mags(j) = M
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        do j = NSKIP+1, STEPS
          avg_e = avg_e + energies(j)
          avg_m = avg_m + abs(mags(j))
          avg_e2 = avg_e2 + (energies(j)**2)
        end do
        write(10,*) itemp+i*dtemp, avg_e/NAVG/N, avg_m/NAVG/N,
     &              (avg_e2/NAVG - (avg_e/NAVG)**2)
     &              /(N*dble(itemp+i*dtemp)**2)
      end do
      close(10)

      END PROGRAM


      SUBROUTINE INITIALIZE(smatrix, SIZE,E,M,JMAG)
      IMPLICIT NONE
      INTEGER i,j,SIZE,smatrix(SIZE,SIZE),JMAG
      REAL*8 sum,E,M

      ! Inicializar la matriz con condiciones T=inf, orientacion aleatoria
      do i=1, SIZE
        do j=1, SIZE
          call random_number(sum)
          smatrix(i,j)=2*nint(sum)-1
        end do
      end do

      ! Calcular magnetizacion y energia iniciales
      E = 0
      M = 0
      do i=1, SIZE ! Iterar sobre cada elemento de la matriz
        do j=1, SIZE
          if (i+1.gt.SIZE) then
            sum = smatrix(1,j)
          else
            sum = smatrix(i+1,j)
          end if
          if (j+1.gt.SIZE) then
            sum = sum + smatrix(i,1)
          else
            sum = sum + smatrix(i,j+1)
          end if
          
          ! Sumar cada elemento a M y E
          M = M + smatrix(i,j)
          E = E - JMAG*smatrix(i,j)*sum
        end do
      end do

      END SUBROUTINE


      SUBROUTINE MONTECARLO_STEP(smatrix, SIZE, E, M, N, padd,
     &           csize)
      IMPLICIT NONE
      INTEGER SIZE,N,smatrix(SIZE, SIZE),csize
      REAL*8 E,M,padd

      ! Variables locales
      REAL*8 randv, dE_total
      INTEGER stack_i(N), stack_j(N)     ! Pila BFS
      INTEGER clust_i(N), clust_j(N)     ! Lista de miembros del cluster
      INTEGER visited(SIZE, SIZE)        ! 1 si el espin esta en el cluster
      INTEGER si, sj, ni, nj            ! Indices actual / vecino
      INTEGER seed_spin                  ! Valor del espin semilla
      INTEGER top, nclust, c            ! Puntero de pila, tamano del cluster

      ! 1. Elegir un espin semilla aleatorio
      call random_number(randv)
      si = int(randv*SIZE) + 1
      call random_number(randv)
      sj = int(randv*SIZE) + 1
      seed_spin = smatrix(si, sj)

      ! 2. Inicializar pila BFS con la semilla
      visited = 0
      visited(si, sj) = 1
      top = 1
      stack_i(1) = si
      stack_j(1) = sj
      nclust = 0

      ! 3. Crecer cluster: sacar un espin, intentar anadir sus 4 vecinos
      do while (top .gt. 0)
        si = stack_i(top)
        sj = stack_j(top)
        top = top - 1
        nclust = nclust + 1
        clust_i(nclust) = si
        clust_j(nclust) = sj

        ! Vecino derecho (CC periodicas via MOD)
        ni = mod(si, SIZE) + 1
        nj = sj
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call random_number(randv)
          if (randv.lt.padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if

        ! Vecino izquierdo
        ni = mod(si - 2 + SIZE, SIZE) + 1
        nj = sj
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call random_number(randv)
          if (randv.lt.padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if

        ! Vecino inferior
        ni = si
        nj = mod(sj, SIZE) + 1
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call random_number(randv)
          if (randv.lt.padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if

        ! Vecino superior
        ni = si
        nj = mod(sj - 2 + SIZE, SIZE) + 1
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call random_number(randv)
          if (randv.lt.padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if
      end do

      ! 4. Calcular cambio de energia en enlaces frontera, luego voltear cluster.
      ! Solo los enlaces que cruzan la frontera del cluster contribuyen a dE.
      ! Los enlaces interiores (ambos espines voltean) se cancelan.
      dE_total = 0.0d0
      do c = 1, nclust
        si = clust_i(c)
        sj = clust_j(c)

        ! Enlace frontera derecho
        ni = mod(si, SIZE) + 1
        if (visited(ni,sj).eq.0) then
          dE_total = dE_total + 2.0d0*smatrix(si,sj)*smatrix(ni,sj)
        end if
        ! Enlace frontera izquierdo
        ni = mod(si - 2 + SIZE, SIZE) + 1
        if (visited(ni,sj).eq.0) then
          dE_total = dE_total + 2.0d0*smatrix(si,sj)*smatrix(ni,sj)
        end if
        ! Enlace frontera inferior
        nj = mod(sj, SIZE) + 1
        if (visited(si,nj).eq.0) then
          dE_total = dE_total + 2.0d0*smatrix(si,sj)*smatrix(si,nj)
        end if
        ! Enlace frontera superior
        nj = mod(sj - 2 + SIZE, SIZE) + 1
        if (visited(si,nj).eq.0) then
          dE_total = dE_total + 2.0d0*smatrix(si,sj)*smatrix(si,nj)
        end if

        ! Voltear el espin
        smatrix(si,sj) = -smatrix(si,sj)
      end do

      ! 5. Actualizar energia total y magnetizacion
      E = E + dE_total
      M = M - 2.0d0 * seed_spin * nclust
      csize = nclust
      
      END SUBROUTINE