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
        ! call RUN_METROPOLIS(L, 10)
        print*, '  Glauber...'
        ! call RUN_GLAUBER(L, 11)
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
      !  Paso grueso 0.2 fuera de [2, 2.5]
      !  Paso fino   0.01 dentro de [2, 2.5]
      ! ===========================================
      SUBROUTINE BUILD_TEMPS(temps, NTEMPS)
      IMPLICIT NONE
      INTEGER NTEMPS, it
      REAL*8 temps(200)

      NTEMPS = 0
      ! Grueso alto: 5.0, 4.8, ..., 2.6
      do it = 0, 12
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = 5.0d0 - it*0.2d0
      end do
      ! Fino: 2.50, 2.49, ..., 2.00
      do it = 0, 50
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = 2.50d0 - it*0.01d0
      end do
      ! Grueso bajo: 1.8, 1.6, ..., 0.2
      do it = 0, 8
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = 1.8d0 - it*0.2d0
      end do

      END SUBROUTINE


      ! ===========================================
      !          EJECUTAR METROPOLIS
      ! ===========================================
      SUBROUTINE RUN_METROPOLIS(L, iunit)
      IMPLICIT NONE
      INTEGER L, iunit, N, DIM, JMAG, STEPS, NSKIP, NAVG
      PARAMETER(DIM=2, JMAG=1, STEPS=20000, NSKIP=6000)
      PARAMETER(NAVG=STEPS-NSKIP)
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8, ALLOCATABLE :: energies(:), mags(:)
      REAL*8 exponentials(0:DIM), E, M
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4
      REAL*8 temps(200), T
      INTEGER NTEMPS, i, j, k

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(energies(STEPS), mags(STEPS))

      call BUILD_TEMPS(temps, NTEMPS)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)

        do k = 0, DIM
          exponentials(k) =
     &      exp(-2.0d0*JMAG*2.0d0*k/T)
        end do

        energies = 0
        mags = 0
        do j = 1, STEPS
          call METROPOLIS_MC_STEP(matrix, L, E, M,
     &                           N, DIM, exponentials)
          energies(j) = E
          mags(j) = M
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        avg_m2 = 0
        avg_m4 = 0
        do j = NSKIP+1, STEPS
          avg_e = avg_e + energies(j)
          avg_m = avg_m + abs(mags(j))
          avg_e2 = avg_e2 + (energies(j)**2)
          avg_m2 = avg_m2 + (mags(j)/N)**2
          avg_m4 = avg_m4 + (mags(j)/N)**4
        end do
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)
     &    /(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG
      end do

      DEALLOCATE(matrix, energies, mags)

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
      INTEGER L, iunit, N, DIM, JMAG, STEPS, NSKIP, NAVG
      PARAMETER(DIM=2, JMAG=1, STEPS=20000, NSKIP=6000)
      PARAMETER(NAVG=STEPS-NSKIP)
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8, ALLOCATABLE :: energies(:), mags(:)
      REAL*8 exponentials(-DIM:DIM), E, M
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4
      REAL*8 temps(200), T
      INTEGER NTEMPS, i, j, k

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(energies(STEPS), mags(STEPS))

      call BUILD_TEMPS(temps, NTEMPS)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)

        do k = -DIM, DIM
          exponentials(k) =
     &      1.0d0/(exp(2.0d0*JMAG*2.0d0*k/T)+1.0d0)
        end do

        energies = 0
        mags = 0
        do j = 1, STEPS
          call GLAUBER_MC_STEP(matrix, L, E, M,
     &                        N, DIM, exponentials)
          energies(j) = E
          mags(j) = M
        end do

        avg_e = 0
        avg_m = 0
        avg_e2 = 0
        avg_m2 = 0
        avg_m4 = 0
        do j = NSKIP+1, STEPS
          avg_e = avg_e + energies(j)
          avg_m = avg_m + abs(mags(j))
          avg_e2 = avg_e2 + (energies(j)**2)
          avg_m2 = avg_m2 + (mags(j)/N)**2
          avg_m4 = avg_m4 + (mags(j)/N)**4
        end do
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)
     &    /(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG
      end do

      DEALLOCATE(matrix, energies, mags)

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
      PARAMETER(DIM=2, JMAG=1, STEPS=20000, NSKIP=6000)
      PARAMETER(NAVG=STEPS-NSKIP)
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8, ALLOCATABLE :: energies(:), mags(:)
      REAL*8 E, M, padd, avg_e, avg_m, avg_e2
      REAL*8 avg_m2, avg_m4
      REAL*8 temps(200), T
      INTEGER NTEMPS, i, j, csize, flipped
      INTEGER, ALLOCATABLE :: stack_i(:), stack_j(:)
      INTEGER, ALLOCATABLE :: clust_i(:), clust_j(:)
      INTEGER, ALLOCATABLE :: visited(:,:)

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(energies(STEPS), mags(STEPS))
      ALLOCATE(stack_i(N), stack_j(N))
      ALLOCATE(clust_i(N), clust_j(N))
      ALLOCATE(visited(L,L))

      call BUILD_TEMPS(temps, NTEMPS)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        print*,T
        padd = 1.0d0 - exp(-2.0d0/T)

        energies = 0
        mags = 0
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
        do j = NSKIP+1, STEPS
          avg_e = avg_e + energies(j)
          avg_m = avg_m + abs(mags(j))
          avg_e2 = avg_e2 + (energies(j)**2)
          avg_m2 = avg_m2 + (mags(j)/N)**2
          avg_m4 = avg_m4 + (mags(j)/N)**4
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