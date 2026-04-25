      ! ===========================================
      !
      !  MODELO DE ISING 2D - SIMULADOR UNIFICADO
      !  Metropolis, Glauber y Wolff
      !  Para L = 16, 32, 64, 128
      !
      ! ===========================================
      PROGRAM SIMULATOR
      IMPLICIT NONE
      INTEGER iL, L
      INTEGER L_values(4)
      DATA L_values /16, 32, 64, 128/

      ! Abrir un fichero por cada tipo de simulacion
      open(10, File='metropolis_2d.txt')
      open(11, File='glauber_2d.txt')
      open(12, File='wolff_2d.txt')

      do iL = 1, 4
        L = L_values(iL)
        print*, 'L =', L
        print*, '  Metropolis...'
        call RUN_METROPOLIS(L, 10)
        print*, '  Glauber...'
        call RUN_GLAUBER(L, 11)
        print*, '  Wolff...'
        call RUN_WOLFF(L, 12)
      end do

      close(10)
      close(11)
      close(12)

      END PROGRAM


      ! ===========================================
      !   SUBRUTINA DE INICIALIZACION DE LA RED
      ! ===========================================
      SUBROUTINE INITIALIZE(matrix, SIZE, E, M, JMAG)
      IMPLICIT NONE
      INTEGER i, j, SIZE, matrix(SIZE,SIZE), JMAG
      REAL*8 sum, E, M

      do i = 1, SIZE
        do j = 1, SIZE
          call random_number(sum)
          matrix(i,j) = 2*nint(sum) - 1
        end do
      end do

      E = 0
      M = 0
      do i = 1, SIZE
        do j = 1, SIZE
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
          M = M + matrix(i,j)
          E = E - JMAG*matrix(i,j)*sum
        end do
      end do

      END SUBROUTINE


      ! ===========================================
      !  CONSTRUIR ARRAY DE TEMPERATURAS (ALTO A BAJO)
      !  Paso grueso 0.1 fuera de [2.2, 2.3]
      !  Paso fino   0.001 dentro de [2, 2.5]
      ! ===========================================
      SUBROUTINE BUILD_TEMPS(temps, NTEMPS, L)
      IMPLICIT NONE
      INTEGER NTEMPS, L
      REAL*8 temps(2000), T_min, T_max, Tc, t
      INTEGER i, j
      REAL*8 temp_val

      Tc = 2.269d0
      T_min = Tc - 2.0d0 / dble(L)
      T_max = Tc + 3.0d0 / dble(L)

      NTEMPS = 0

      ! High temperatures: 5.0 down to T_max (step 0.1)
      t = 5.0d0
      do while (t .ge. T_max + 0.05d0)
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = t
        t = t - 0.1d0
      end do

      ! Fine temperatures: T_max down to T_min (step 0.005)
      ! Start from the first multiple of 0.005 <= T_max
      t = T_max
      do while (t .ge. T_min - 0.0025d0)
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = t
        t = t - 0.005d0
      end do

      ! Low temperatures: T_min down to 0.2 (step 0.1)
      ! Start from the first multiple of 0.1 <= T_min
      t = int(T_min * 10.0d0) / 10.0d0
      if (t .gt. T_min - 0.001d0) t = t - 0.1d0
      do while (t .ge. 0.15d0)
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = t
        t = t - 0.1d0
      end do

      END SUBROUTINE


      ! ===========================================
      !          EJECUTAR METROPOLIS
      ! ===========================================
      SUBROUTINE RUN_METROPOLIS(L, iunit)
      IMPLICIT NONE
      INTEGER L, iunit, N, DIM, JMAG
      PARAMETER(DIM=2, JMAG=1)
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8 exponentials(0:DIM), E, M
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4
      REAL*8 temps(2000), T, Tc, T_min, T_max
      INTEGER NTEMPS, STEPS, NSKIP, NAVG, i, j, k

      N = L*L
      ALLOCATE(matrix(L,L))

      Tc = 2.269d0
      T_min = Tc - 2.0d0 / dble(L)
      T_max = Tc + 3.0d0 / dble(L)

      call BUILD_TEMPS(temps, NTEMPS, L)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        print*,T

        ! Mas pasos cerca de Tc para reducir ruido
        if (T.ge.2.0d0 .and. T.le.2.5d0) then
          STEPS = 2000000
          NSKIP = 600000
        else
          STEPS = 20000
          NSKIP = 6000
        end if
        NAVG = 0

        do k = 0, DIM
          exponentials(k) =
     &      exp(-2.0d0*JMAG*2.0d0*k/T)
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        avg_m2 = 0
        avg_m4 = 0
        do j = 1, STEPS
          call METROPOLIS_MC_STEP(matrix, L, E, M,
     &                           N, DIM, exponentials)
          if (j .gt. NSKIP.and.(mod(j,50).eq.0)) then
            avg_e = avg_e + E
            avg_m = avg_m + abs(M)
            avg_e2 = avg_e2 + E**2
            avg_m2 = avg_m2 + (M/N)**2
            avg_m4 = avg_m4 + (M/N)**4
            NAVG = NAVG + 1
          end if
        end do
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)
     &    /(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG
      end do

      DEALLOCATE(matrix)

      END SUBROUTINE


      ! ===========================================
      !     PASO MONTE CARLO - METROPOLIS
      ! ===========================================
      SUBROUTINE METROPOLIS_MC_STEP(matrix, SIZE,
     &            E, M, N, DIM, exponentials)
      IMPLICIT NONE
      INTEGER SIZE, DIM, N, i, j, k, matrix(SIZE,SIZE)
      REAL*8 E, M, dE, sum, exponentials(0:DIM), randv

      do k = 1, N
        call random_number(randv)
        j = int(randv*SIZE) + 1
        call random_number(randv)
        i = int(randv*SIZE) + 1

        sum = 0.0
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

        dE = 2*matrix(i,j)*sum
        if (dE.le.0) then
          matrix(i,j) = -matrix(i,j)
          E = E + dE
          M = M + matrix(i,j)*2
        else
          call random_number(randv)
          if (randv.lt.exponentials(int(dE/2/2)))
     &    then
            matrix(i,j) = -matrix(i,j)
            E = E + dE
            M = M + matrix(i,j)*2
          end if
        end if
      end do

      END SUBROUTINE


      ! ===========================================
      !          EJECUTAR GLAUBER
      ! ===========================================
      SUBROUTINE RUN_GLAUBER(L, iunit)
      IMPLICIT NONE
      INTEGER L, iunit, N, DIM, JMAG
      PARAMETER(DIM=2, JMAG=1)
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8 exponentials(-DIM:DIM), E, M
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4
      REAL*8 temps(2000), T, Tc, T_min, T_max
      INTEGER NTEMPS, STEPS, NSKIP, NAVG, i, j, k

      N = L*L
      ALLOCATE(matrix(L,L))

      Tc = 2.269d0
      T_min = Tc - 2.0d0 / dble(L)
      T_max = Tc + 3.0d0 / dble(L)

      call BUILD_TEMPS(temps, NTEMPS, L)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        print*,T
        ! Mas pasos cerca de Tc para reducir ruido
        if (T.ge.2.0d0 .and. T.le.2.5d0) then
          STEPS = 2000000
          NSKIP = 600000
        else
          STEPS = 20000
          NSKIP = 6000
        end if
        NAVG = 0

        do k = -DIM, DIM
          exponentials(k) =
     &      1.0d0/(exp(2.0d0*JMAG*2.0d0*k/T)+1.0d0)
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        avg_m2 = 0
        avg_m4 = 0
        do j = 1, STEPS
          call GLAUBER_MC_STEP(matrix, L, E, M,
     &                        N, DIM, exponentials)
          if (j .gt. NSKIP.and.mod(j,50).eq.0) then
            avg_e = avg_e + E
            avg_m = avg_m + abs(M)
            avg_e2 = avg_e2 + E**2
            avg_m2 = avg_m2 + (M/N)**2
            avg_m4 = avg_m4 + (M/N)**4
            NAVG = NAVG + 1
          end if
        end do
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)
     &    /(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG
      end do

      DEALLOCATE(matrix)

      END SUBROUTINE


      ! ===========================================
      !      PASO MONTE CARLO - GLAUBER
      ! ===========================================
      SUBROUTINE GLAUBER_MC_STEP(matrix, SIZE,
     &            E, M, N, DIM, exponentials)
      IMPLICIT NONE
      INTEGER SIZE, DIM, N, i, j, k, matrix(SIZE,SIZE)
      REAL*8 E,M,dE,sum,exponentials(-DIM:DIM),randv

      do k = 1, N
        call random_number(randv)
        j = int(randv*SIZE) + 1
        call random_number(randv)
        i = int(randv*SIZE) + 1

        sum = 0.0
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

        dE = 2*matrix(i,j)*sum
        call random_number(randv)
        if (randv.lt.exponentials(int(dE/2/2)))
     &  then
          matrix(i,j) = -matrix(i,j)
          E = E + dE
          M = M + matrix(i,j)*2
        end if
      end do

      END SUBROUTINE


      ! ===========================================
      !           EJECUTAR WOLFF
      ! ===========================================
      SUBROUTINE RUN_WOLFF(L, iunit)
      IMPLICIT NONE
      INTEGER L, iunit, N, DIM, JMAG, STEPS, NSKIP, NAVG
      INTEGER STEPS_MAX
      PARAMETER(DIM=2, JMAG=1, STEPS_MAX=2000000)
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8, ALLOCATABLE :: energies(:), mags(:)
      REAL*8 E, M, padd, avg_e, avg_m, avg_e2
      REAL*8 avg_m2, avg_m4
      REAL*8 temps(2000), T, Tc, T_min, T_max
      INTEGER NTEMPS, i, j, csize, flipped
      INTEGER, ALLOCATABLE :: stack_i(:), stack_j(:)
      INTEGER, ALLOCATABLE :: clust_i(:), clust_j(:)
      INTEGER, ALLOCATABLE :: visited(:,:)

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(energies(STEPS_MAX), mags(STEPS_MAX))
      ALLOCATE(stack_i(N), stack_j(N))
      ALLOCATE(clust_i(N), clust_j(N))
      ALLOCATE(visited(L,L))

      Tc = 2.269d0
      T_min = Tc - 2.0d0 / dble(L)
      T_max = Tc + 3.0d0 / dble(L)

      call BUILD_TEMPS(temps, NTEMPS, L)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        print*,T
        ! Mas pasos cerca de Tc para reducir ruido
        if (T.ge.2.0d0 .and. T.le.2.5d0) then
          STEPS = 2000000
          NSKIP = 600000
        else
          STEPS = 20000
          NSKIP = 6000
        end if

        if (STEPS.gt.STEPS_MAX) then
          print*, 'ERROR: STEPS > STEPS_MAX en RUN_WOLFF'
          stop
        end if
        
        padd = 1.0d0 - exp(-2.0d0/T)

        energies(1:STEPS) = 0
        mags(1:STEPS) = 0
        ! Cada barrido MC = suficientes volteos de cluster
        ! para voltear ~N espines en total
        do j = 1, STEPS
          flipped = 0
          do while (flipped .lt. N)
            call WOLFF_MC_STEP(matrix, L, E, M,
     &                         N, padd, csize,
     &                         stack_i, stack_j,
     &                         clust_i, clust_j,
     &                         visited)
            flipped = flipped + csize
          end do
          energies(j) = E
          mags(j) = M
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        avg_m2 = 0
        avg_m4 = 0
        NAVG = 0
        do j = NSKIP+1, STEPS
          if (mod(j,5).eq.0) then
            avg_e = avg_e + energies(j)
            avg_m = avg_m + abs(mags(j))
            avg_e2 = avg_e2 + (energies(j)**2)
            avg_m2 = avg_m2 + (mags(j)/N)**2
            avg_m4 = avg_m4 + (mags(j)/N)**4
            NAVG = NAVG + 1
          end if
        end do
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)
     &    /(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG
      end do

      DEALLOCATE(matrix, energies, mags)
      DEALLOCATE(stack_i, stack_j, clust_i, clust_j, visited)

      END SUBROUTINE


      ! ===========================================
      !       PASO MONTE CARLO - WOLFF
      ! ===========================================
      SUBROUTINE WOLFF_MC_STEP(smatrix, SIZE,
     &            E, M, N, padd, csize,
     &            stack_i, stack_j,
     &            clust_i, clust_j,
     &            visited)
      IMPLICIT NONE
      INTEGER SIZE, N, smatrix(SIZE, SIZE), csize
      REAL*8 E, M, padd
      REAL*8 randv, dE_total
      INTEGER stack_i(N), stack_j(N)
      INTEGER clust_i(N), clust_j(N)
      INTEGER visited(SIZE, SIZE)
      INTEGER si, sj, ni, nj
      INTEGER seed_spin
      INTEGER top, nclust, c

      call random_number(randv)
      si = int(randv*SIZE) + 1
      call random_number(randv)
      sj = int(randv*SIZE) + 1
      seed_spin = smatrix(si, sj)

      visited = 0
      visited(si, sj) = 1
      top = 1
      stack_i(1) = si
      stack_j(1) = sj
      nclust = 0

      do while (top .gt. 0)
        si = stack_i(top)
        sj = stack_j(top)
        top = top - 1
        nclust = nclust + 1
        clust_i(nclust) = si
        clust_j(nclust) = sj

        ! Vecino derecho
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

      ! Calcular cambio de energia y voltear cluster
      dE_total = 0.0d0
      do c = 1, nclust
        si = clust_i(c)
        sj = clust_j(c)

        ni = mod(si, SIZE) + 1
        if (visited(ni,sj).eq.0) then
          dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(ni,sj)
        end if
        ni = mod(si - 2 + SIZE, SIZE) + 1
        if (visited(ni,sj).eq.0) then
          dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(ni,sj)
        end if
        nj = mod(sj, SIZE) + 1
        if (visited(si,nj).eq.0) then
          dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(si,nj)
        end if
        nj = mod(sj - 2 + SIZE, SIZE) + 1
        if (visited(si,nj).eq.0) then
          dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(si,nj)
        end if

        smatrix(si,sj) = -smatrix(si,sj)
      end do

      E = E + dE_total
      M = M - 2.0d0 * seed_spin * nclust
      csize = nclust

      END SUBROUTINE