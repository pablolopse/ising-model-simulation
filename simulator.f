      ! ===========================================
      !
      !  MODELO DE ISING 2D - SIMULADOR UNIFICADO
      !  Metropolis, Glauber y Wolff
      !
      !  Uso: ./simulator <algoritmo> <L> [STEPS [NSKIP]]
      !       algoritmo = metropolis | glauber | wolff
      !       L         = tamano de red (ej. 16, 32, 64, 128)
      !       STEPS     = pasos MC por temperatura (opcional)
      !       NSKIP     = pasos de termalización (opcional)
      !
      ! ===========================================
      PROGRAM SIMULATOR
      IMPLICIT NONE

      ! =============================================
      !   CONFIGURACION: editar solo este bloque
      ! =============================================

      REAL*8 TC_EXACT
      PARAMETER (TC_EXACT = 2.269d0)

      REAL*8 T_START, T_MIN, DT_COARSE, WIN_LO, WIN_HI
      PARAMETER (T_START = 5.0d0, T_MIN = 0.15d0)
      PARAMETER (DT_COARSE = 0.1d0)
      PARAMETER (WIN_LO = 1.0d0, WIN_HI = 3.0d0)

      INTEGER STEPS_DEF, NSKIP_DEF
      PARAMETER (STEPS_DEF = 20000, NSKIP_DEF = 6000)
      ! =============================================

      INTEGER*8 prng_s
      COMMON /PRNG/ prng_s

      CHARACTER*20 arg_alg
      CHARACTER*10 arg_L_str, arg_steps, arg_nskip
      CHARACTER*40 fname
      INTEGER L, steps_run, nskip_run, nargs
      INTEGER*8 clk8

      ! Parse arguments first so we can use them in the seed
      call get_command_argument(1, arg_alg)
      call get_command_argument(2, arg_L_str)
      read(arg_L_str, *) L

      ! Optional: override step counts from command line
      nargs = command_argument_count()
      if (nargs .ge. 3) then
        call get_command_argument(3, arg_steps)
        read(arg_steps, *) steps_run
      else
        steps_run = STEPS_DEF
      end if
      if (nargs .ge. 4) then
        call get_command_argument(4, arg_nskip)
        read(arg_nskip, *) nskip_run
      else
        nskip_run = NSKIP_DEF
      end if

      ! Seed: clock provides entropy across runs; L and algorithm
      ! first letter guarantee different seeds for each of the 12 tasks
      ! even if the clock reads the same value.
      call system_clock(clk8)
      prng_s = clk8 + int(L, 8) * 2654435761_8
      prng_s = ieor(prng_s,
     &  int(ichar(arg_alg(1:1)), 8) * 6364136223846793005_8)
      if (prng_s .eq. 0_8) prng_s = 1_8

      write(fname, '("tmp_",A,"_",I0,".dat")') trim(arg_alg), L
      open(10, File=fname)

      if (trim(arg_alg) .eq. 'metropolis') then
        call RUN_METROPOLIS(L, 10,
     &    TC_EXACT, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE,
     &    steps_run, nskip_run)
      else if (trim(arg_alg) .eq. 'glauber') then
        call RUN_GLAUBER(L, 10,
     &    TC_EXACT, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE,
     &    steps_run, nskip_run)
      else if (trim(arg_alg) .eq. 'wolff') then
        call RUN_WOLFF(L, 10,
     &    TC_EXACT, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE,
     &    steps_run, nskip_run)
      else
        write(*,*) 'Algoritmo desconocido: ', trim(arg_alg)
        write(*,*) 'Uso: ./simulator metropolis|glauber|wolff <L>'
        stop 1
      end if

      close(10)

      END PROGRAM


      ! ===========================================
      !   GENERADOR ALEATORIO RAPIDO (xorshift64)
      !
      !   Periodo 2^64-1. 3 operaciones vs las ~700
      !   del Mersenne Twister de random_number().
      ! ===========================================
      SUBROUTINE FAST_RAND(r)
      IMPLICIT NONE
      REAL*8 r
      INTEGER*8 s
      COMMON /PRNG/ s

      s = ieor(s, ishft(s,  13))
      s = ieor(s, ishft(s,  -7))
      s = ieor(s, ishft(s,  17))
      ! Map lower 53 bits to [0, 1) exactly (double precision mantissa)
      r = dble(iand(s, 9007199254740991_8)) * 1.1102230246251565d-16

      END SUBROUTINE


      ! ===========================================
      !   INICIALIZACION DE LA RED
      !
      !   Red ordenada (todos espines +1).
      ! ===========================================
      SUBROUTINE INITIALIZE(matrix, SIZE, E, M, JMAG)
      IMPLICIT NONE
      INTEGER i, j, SIZE, matrix(SIZE,SIZE), JMAG
      REAL*8 E, M

      do i = 1, SIZE
        do j = 1, SIZE
          matrix(i,j) = 1
        end do
      end do

      E = -dble(JMAG) * 2.0d0 * dble(SIZE) * dble(SIZE)
      M = dble(SIZE) * dble(SIZE)

      END SUBROUTINE


      ! ===========================================
      !   REJILLA DE TEMPERATURAS (BAJO A ALTO)
      ! ===========================================
      SUBROUTINE BUILD_TEMPS(temps, NTEMPS, L,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE)
      IMPLICIT NONE
      INTEGER NTEMPS, L
      REAL*8 temps(2000), TC, WIN_LO, WIN_HI
      REAL*8 T_START, T_MIN, DT_COARSE
      REAL*8 T_lo, T_hi, t, DT_FINE

      T_lo    = TC - (WIN_LO * TC) / dble(L)
      T_hi    = TC + (WIN_HI * TC) / dble(L)
      DT_FINE = (DT_COARSE * TC) / dble(L)

      NTEMPS = 0

      t = T_MIN
      do while (t .le. T_lo - 0.5d0*DT_COARSE)
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = t
        t = t + DT_COARSE
      end do

      t = T_lo
      do while (t .le. T_hi + 0.5d0*DT_FINE)
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = t
        t = t + DT_FINE
      end do

      t = int(T_hi / DT_COARSE) * DT_COARSE
      if (t .lt. T_hi + 0.001d0) t = t + DT_COARSE
      do while (t .le. T_START + 0.5d0*DT_COARSE)
        NTEMPS = NTEMPS + 1
        temps(NTEMPS) = t
        t = t + DT_COARSE
      end do

      END SUBROUTINE


      ! ===========================================
      !   ACUMULAR OBSERVABLES EN EL PASO j
      ! ===========================================
      SUBROUTINE ACCUMULATE_OBS(j, NSKIP, NSAMPLE,
     &  matrix, SIZE, N, E, M,
     &  cos_tab, sin_tab,
     &  avg_e, avg_m, avg_e2,
     &  avg_m2, avg_m4, avg_mk2, NAVG)
      IMPLICIT NONE
      INTEGER j, NSKIP, NSAMPLE, SIZE, N, NAVG
      INTEGER matrix(SIZE,SIZE)
      REAL*8 E, M
      REAL*8 cos_tab(SIZE), sin_tab(SIZE)
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4, avg_mk2
      REAL*8 mk_re, mk_im
      INTEGER ix, iy

      if (j .gt. NSKIP .and. mod(j, NSAMPLE) .eq. 0) then
        avg_e  = avg_e  + E
        avg_m  = avg_m  + abs(M)
        avg_e2 = avg_e2 + E**2
        avg_m2 = avg_m2 + (M/N)**2
        avg_m4 = avg_m4 + (M/N)**4
        mk_re = 0.0d0
        mk_im = 0.0d0
        do iy = 1, SIZE
          do ix = 1, SIZE
            mk_re = mk_re + matrix(ix,iy) * cos_tab(ix)
            mk_im = mk_im + matrix(ix,iy) * sin_tab(ix)
          end do
        end do
        avg_mk2 = avg_mk2 +
     &    (mk_re**2 + mk_im**2) / dble(N)**2
        NAVG = NAVG + 1
      end if

      END SUBROUTINE


      ! ===========================================
      !          EJECUTAR METROPOLIS
      ! ===========================================
      SUBROUTINE RUN_METROPOLIS(L, iunit,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE,
     &  STEPS, NSKIP)
      IMPLICIT NONE
      INTEGER L, iunit, N, DIM, JMAG
      PARAMETER (DIM=2, JMAG=1)
      REAL*8 TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE
      INTEGER STEPS, NSKIP

      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8 exponentials(0:DIM)
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4, avg_mk2
      REAL*8, ALLOCATABLE :: cos_tab(:), sin_tab(:)
      REAL*8, ALLOCATABLE :: m_series(:)
      REAL*8 temps(2000), T, E, M, pi, tau_int
      INTEGER NTEMPS, NAVG, i, j, k, ix, nsamp
      INTEGER, PARAMETER :: NSAMPLE = 50

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(m_series((STEPS-NSKIP)/NSAMPLE + 10))
      pi = 4.0d0 * atan(1.0d0)
      ALLOCATE(cos_tab(L), sin_tab(L))
      do ix = 1, L
        cos_tab(ix) = cos(2.0d0*pi*dble(ix)/dble(L))
        sin_tab(ix) = sin(2.0d0*pi*dble(ix)/dble(L))
      end do

      call BUILD_TEMPS(temps, NTEMPS, L,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        write(0, '("[",I4,"/",I4,"] T=",F6.3)') i, NTEMPS, T
        call flush(0)

        do k = 0, DIM
          exponentials(k) =
     &      exp(-2.0d0*JMAG*2.0d0*dble(k)/T)
        end do

        avg_e   = 0.0d0
        avg_m   = 0.0d0
        avg_e2  = 0.0d0
        avg_m2  = 0.0d0
        avg_m4  = 0.0d0
        avg_mk2 = 0.0d0
        NAVG    = 0
        nsamp   = 0
        do j = 1, STEPS
          call METROPOLIS_MC_STEP(matrix, L, E, M,
     &      N, DIM, exponentials)
          call ACCUMULATE_OBS(j, NSKIP, NSAMPLE,
     &      matrix, L, N, E, M, cos_tab, sin_tab,
     &      avg_e, avg_m, avg_e2, avg_m2, avg_m4,
     &      avg_mk2, NAVG)
          if (j .gt. NSKIP .and. mod(j, NSAMPLE) .eq. 0) then
            nsamp = nsamp + 1
            m_series(nsamp) = abs(M) / dble(N)
          end if
        end do
        call COMPUTE_TAU_INT(m_series, NAVG, tau_int)
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)/(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG, avg_mk2/NAVG,
     &    tau_int * NSAMPLE
      end do

      DEALLOCATE(matrix, cos_tab, sin_tab, m_series)

      END SUBROUTINE


      ! ===========================================
      !     PASO MONTE CARLO - METROPOLIS
      ! ===========================================
      SUBROUTINE METROPOLIS_MC_STEP(matrix, SIZE,
     &  E, M, N, DIM, exponentials)
      IMPLICIT NONE
      INTEGER SIZE, DIM, N, i, j, k, matrix(SIZE,SIZE)
      REAL*8 E, M, dE, sum, exponentials(0:DIM), randv
      INTEGER ni_r, ni_l, nj_d, nj_u, kE

      do k = 1, N
        call FAST_RAND(randv)
        i = int(randv*SIZE) + 1
        call FAST_RAND(randv)
        j = int(randv*SIZE) + 1

        ni_r = mod(i,        SIZE) + 1
        ni_l = mod(i-2+SIZE, SIZE) + 1
        nj_d = mod(j,        SIZE) + 1
        nj_u = mod(j-2+SIZE, SIZE) + 1
        sum = dble(matrix(ni_r,j) + matrix(ni_l,j)
     &           + matrix(i,nj_d) + matrix(i,nj_u))

        dE = 2.0d0 * matrix(i,j) * sum
        if (dE .le. 0.0d0) then
          matrix(i,j) = -matrix(i,j)
          E = E + dE
          M = M + 2*matrix(i,j)
        else
          kE = nint(dE / 4.0d0)
          call FAST_RAND(randv)
          if (randv .lt. exponentials(kE)) then
            matrix(i,j) = -matrix(i,j)
            E = E + dE
            M = M + 2*matrix(i,j)
          end if
        end if
      end do

      END SUBROUTINE


      ! ===========================================
      !          EJECUTAR GLAUBER
      ! ===========================================
      SUBROUTINE RUN_GLAUBER(L, iunit,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE,
     &  STEPS, NSKIP)
      IMPLICIT NONE
      INTEGER L, iunit, N, DIM, JMAG
      PARAMETER (DIM=2, JMAG=1)
      REAL*8 TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE
      INTEGER STEPS, NSKIP

      INTEGER, PARAMETER :: NSAMPLE = 50
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8 exponentials(-DIM:DIM)
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4, avg_mk2
      REAL*8, ALLOCATABLE :: cos_tab(:), sin_tab(:)
      REAL*8, ALLOCATABLE :: m_series(:)
      REAL*8 temps(2000), T, E, M, pi, tau_int
      INTEGER NTEMPS, NAVG, i, j, k, ix, nsamp

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(m_series((STEPS-NSKIP)/NSAMPLE + 10))
      pi = 4.0d0 * atan(1.0d0)
      ALLOCATE(cos_tab(L), sin_tab(L))
      do ix = 1, L
        cos_tab(ix) = cos(2.0d0*pi*dble(ix)/dble(L))
        sin_tab(ix) = sin(2.0d0*pi*dble(ix)/dble(L))
      end do

      call BUILD_TEMPS(temps, NTEMPS, L,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        write(0, '("[",I4,"/",I4,"] T=",F6.3)') i, NTEMPS, T
        call flush(0)

        do k = -DIM, DIM
          exponentials(k) =
     &      1.0d0 / (exp(2.0d0*JMAG*2.0d0*dble(k)/T)
     &               + 1.0d0)
        end do

        avg_e   = 0.0d0
        avg_m   = 0.0d0
        avg_e2  = 0.0d0
        avg_m2  = 0.0d0
        avg_m4  = 0.0d0
        avg_mk2 = 0.0d0
        NAVG    = 0
        nsamp   = 0
        do j = 1, STEPS
          call GLAUBER_MC_STEP(matrix, L, E, M,
     &      N, DIM, exponentials)
          call ACCUMULATE_OBS(j, NSKIP, NSAMPLE,
     &      matrix, L, N, E, M, cos_tab, sin_tab,
     &      avg_e, avg_m, avg_e2, avg_m2, avg_m4,
     &      avg_mk2, NAVG)
          if (j .gt. NSKIP .and. mod(j, NSAMPLE) .eq. 0) then
            nsamp = nsamp + 1
            m_series(nsamp) = abs(M) / dble(N)
          end if
        end do
        call COMPUTE_TAU_INT(m_series, NAVG, tau_int)
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)/(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG, avg_mk2/NAVG,
     &    tau_int * NSAMPLE
      end do

      DEALLOCATE(matrix, cos_tab, sin_tab, m_series)

      END SUBROUTINE


      ! ===========================================
      !      PASO MONTE CARLO - GLAUBER
      ! ===========================================
      SUBROUTINE GLAUBER_MC_STEP(matrix, SIZE,
     &  E, M, N, DIM, exponentials)
      IMPLICIT NONE
      INTEGER SIZE, DIM, N, i, j, k, matrix(SIZE,SIZE)
      REAL*8 E, M, dE, sum, exponentials(-DIM:DIM), randv
      INTEGER ni_r, ni_l, nj_d, nj_u, kE

      do k = 1, N
        call FAST_RAND(randv)
        i = int(randv*SIZE) + 1
        call FAST_RAND(randv)
        j = int(randv*SIZE) + 1

        ni_r = mod(i,        SIZE) + 1
        ni_l = mod(i-2+SIZE, SIZE) + 1
        nj_d = mod(j,        SIZE) + 1
        nj_u = mod(j-2+SIZE, SIZE) + 1
        sum = dble(matrix(ni_r,j) + matrix(ni_l,j)
     &           + matrix(i,nj_d) + matrix(i,nj_u))

        dE = 2.0d0 * matrix(i,j) * sum
        kE = nint(dE / 4.0d0)
        call FAST_RAND(randv)
        if (randv .lt. exponentials(kE)) then
          matrix(i,j) = -matrix(i,j)
          E = E + dE
          M = M + 2*matrix(i,j)
        end if
      end do

      END SUBROUTINE


      ! ===========================================
      !           EJECUTAR WOLFF
      ! ===========================================
      SUBROUTINE RUN_WOLFF(L, iunit,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE,
     &  STEPS, NSKIP)
      IMPLICIT NONE
      INTEGER L, iunit, N, JMAG
      PARAMETER (JMAG=1)
      REAL*8 TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE
      INTEGER STEPS, NSKIP

      INTEGER, PARAMETER :: NSAMPLE = 5
      INTEGER, ALLOCATABLE :: matrix(:,:)
      REAL*8 E, M, padd
      REAL*8 avg_e, avg_m, avg_e2, avg_m2, avg_m4, avg_mk2
      REAL*8, ALLOCATABLE :: cos_tab(:), sin_tab(:)
      REAL*8, ALLOCATABLE :: m_series(:)
      REAL*8 temps(2000), T, pi, tau_int
      INTEGER NTEMPS, NAVG, nsamp
      INTEGER i, j, csize, flipped, ix
      INTEGER, ALLOCATABLE :: stack_i(:), stack_j(:)
      INTEGER, ALLOCATABLE :: clust_i(:), clust_j(:)
      INTEGER, ALLOCATABLE :: visited(:,:)

      N = L*L
      ALLOCATE(matrix(L,L))
      ALLOCATE(stack_i(N), stack_j(N))
      ALLOCATE(clust_i(N), clust_j(N))
      ALLOCATE(visited(L,L))
      ALLOCATE(m_series((STEPS-NSKIP)/NSAMPLE + 10))
      visited = 0
      pi = 4.0d0 * atan(1.0d0)
      ALLOCATE(cos_tab(L), sin_tab(L))
      do ix = 1, L
        cos_tab(ix) = cos(2.0d0*pi*dble(ix)/dble(L))
        sin_tab(ix) = sin(2.0d0*pi*dble(ix)/dble(L))
      end do

      call BUILD_TEMPS(temps, NTEMPS, L,
     &  TC, WIN_LO, WIN_HI, T_START, T_MIN, DT_COARSE)
      call INITIALIZE(matrix, L, E, M, JMAG)

      do i = 1, NTEMPS
        T = temps(i)
        write(0, '("[",I4,"/",I4,"] T=",F6.3)') i, NTEMPS, T
        call flush(0)

        padd = 1.0d0 - exp(-2.0d0*JMAG/T)

        avg_e   = 0.0d0
        avg_m   = 0.0d0
        avg_e2  = 0.0d0
        avg_m2  = 0.0d0
        avg_m4  = 0.0d0
        avg_mk2 = 0.0d0
        NAVG    = 0
        nsamp   = 0
        do j = 1, STEPS
          flipped = 0
          do while (flipped .lt. N)
            call WOLFF_MC_STEP(matrix, L, E, M,
     &        N, padd, csize,
     &        stack_i, stack_j,
     &        clust_i, clust_j, visited)
            flipped = flipped + csize
          end do
          call ACCUMULATE_OBS(j, NSKIP, NSAMPLE,
     &      matrix, L, N, E, M, cos_tab, sin_tab,
     &      avg_e, avg_m, avg_e2, avg_m2, avg_m4,
     &      avg_mk2, NAVG)
          if (j .gt. NSKIP .and. mod(j, NSAMPLE) .eq. 0) then
            nsamp = nsamp + 1
            m_series(nsamp) = abs(M) / dble(N)
          end if
        end do
        call COMPUTE_TAU_INT(m_series, NAVG, tau_int)
        write(iunit,*) L, T, avg_e/NAVG/N,
     &    avg_m/NAVG/N,
     &    (avg_e2/NAVG - (avg_e/NAVG)**2)/(N*T**2),
     &    avg_m2/NAVG, avg_m4/NAVG, avg_mk2/NAVG,
     &    tau_int * NSAMPLE
      end do

      DEALLOCATE(matrix, stack_i, stack_j)
      DEALLOCATE(clust_i, clust_j, visited)
      DEALLOCATE(cos_tab, sin_tab, m_series)

      END SUBROUTINE


      ! ===========================================
      !       PASO MONTE CARLO - WOLFF
      ! ===========================================
      SUBROUTINE WOLFF_MC_STEP(smatrix, SIZE,
     &  E, M, N, padd, csize,
     &  stack_i, stack_j,
     &  clust_i, clust_j, visited)
      IMPLICIT NONE
      INTEGER SIZE, N, smatrix(SIZE, SIZE), csize
      REAL*8 E, M, padd, randv, dE_total
      INTEGER stack_i(N), stack_j(N)
      INTEGER clust_i(N), clust_j(N)
      INTEGER visited(SIZE, SIZE)
      INTEGER si, sj, ni, nj, seed_spin, top, nclust, c

      call FAST_RAND(randv)
      si = int(randv*SIZE) + 1
      call FAST_RAND(randv)
      sj = int(randv*SIZE) + 1
      seed_spin = smatrix(si, sj)

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
        ni = mod(si,        SIZE) + 1
        nj = sj
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call FAST_RAND(randv)
          if (randv .lt. padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if

        ! Vecino izquierdo
        ni = mod(si-2+SIZE, SIZE) + 1
        nj = sj
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call FAST_RAND(randv)
          if (randv .lt. padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if

        ! Vecino inferior
        ni = si
        nj = mod(sj,        SIZE) + 1
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call FAST_RAND(randv)
          if (randv .lt. padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if

        ! Vecino superior
        ni = si
        nj = mod(sj-2+SIZE, SIZE) + 1
        if (visited(ni,nj).eq.0 .and.
     &      smatrix(ni,nj).eq.seed_spin) then
          call FAST_RAND(randv)
          if (randv .lt. padd) then
            visited(ni,nj) = 1
            top = top + 1
            stack_i(top) = ni
            stack_j(top) = nj
          end if
        end if
      end do

      dE_total = 0.0d0
      do c = 1, nclust
        si = clust_i(c)
        sj = clust_j(c)

        ni = mod(si,        SIZE) + 1
        if (visited(ni,sj).eq.0)
     &    dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(ni,sj)

        ni = mod(si-2+SIZE, SIZE) + 1
        if (visited(ni,sj).eq.0)
     &    dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(ni,sj)

        nj = mod(sj,        SIZE) + 1
        if (visited(si,nj).eq.0)
     &    dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(si,nj)

        nj = mod(sj-2+SIZE, SIZE) + 1
        if (visited(si,nj).eq.0)
     &    dE_total = dE_total +
     &      2.0d0*smatrix(si,sj)*smatrix(si,nj)

        smatrix(si,sj) = -smatrix(si,sj)
      end do

      E = E + dE_total
      M = M - 2.0d0 * seed_spin * nclust
      csize = nclust

      ! Reset only the sites that were marked, not the whole array
      do c = 1, nclust
        visited(clust_i(c), clust_j(c)) = 0
      end do

      END SUBROUTINE


      ! ===========================================
      !   TIEMPO DE AUTOCORRELACION INTEGRADO
      !
      !   Metodo de ventana de Madras-Sokal:
      !   suma rho(k) hasta que k >= 6 * tau_int.
      !   Devuelve tau en unidades de muestras.
      ! ===========================================
      SUBROUTINE COMPUTE_TAU_INT(series, n, tau_int)
      IMPLICIT NONE
      INTEGER n, k, j
      REAL*8 series(n), tau_int
      REAL*8 xmean, c0, ck, sum_rho

      xmean = 0.0d0
      do k = 1, n
        xmean = xmean + series(k)
      end do
      xmean = xmean / dble(n)

      c0 = 0.0d0
      do k = 1, n
        c0 = c0 + (series(k) - xmean)**2
      end do
      c0 = c0 / dble(n)

      if (c0 .le. 0.0d0) then
        tau_int = 0.5d0
        return
      end if

      sum_rho = 0.0d0
      do k = 1, n/2
        ck = 0.0d0
        do j = 1, n - k
          ck = ck + (series(j) - xmean) * (series(j+k) - xmean)
        end do
        ck = ck / dble(n - k)
        sum_rho = sum_rho + ck / c0
        if (dble(k) .ge. 6.0d0 * (0.5d0 + sum_rho)) goto 10
      end do
 10   tau_int = 0.5d0 + sum_rho

      END SUBROUTINE
