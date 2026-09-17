#!/bin/bash
# Uso: ./run.sh [--smoke]
#   --smoke   Solo L=16,32 y pocos pasos (verificacion rapida)
set -e

ALGOS=(metropolis glauber wolff)

if [ "${1:-}" = "--smoke" ]; then
  L_VALUES=(16 32)
  EXTRA_ARGS="20000 6000"   # STEPS NSKIP
else
  L_VALUES=(16 32 64 128)
  EXTRA_ARGS=""
fi

N_TASKS=$(( ${#ALGOS[@]} * ${#L_VALUES[@]} ))

T_START=$SECONDS

gfortran -O3 -march=native -ffast-math simulator.f -o simulator

# Clean up any leftover files from previous runs
rm -f status_*.tmp done_*.flag tmp_*.dat

echo "Ejecutando $N_TASKS tareas en paralelo..."
echo ""

# Print initial status block (N_TASKS lines the monitor will overwrite)
for alg in "${ALGOS[@]}"; do
  for L in "${L_VALUES[@]}"; do
    printf "%-12s L=%4d: starting...\n" "$alg" "$L"
  done
done

# Background monitor: redraws the status block every 0.5 s
monitor_progress() {
  while true; do
    sleep 0.5
    printf "\033[${N_TASKS}A"   # move cursor up N_TASKS lines
    for alg in "${ALGOS[@]}"; do
      for L in "${L_VALUES[@]}"; do
        if [ -f "done_${alg}_${L}.flag" ]; then
          status="[HECHO]"
        else
          line=$(tail -1 "status_${alg}_${L}.tmp" 2>/dev/null)
          status="${line:-starting...}"
        fi
        printf "%-12s L=%4d: %-45s\n" "$alg" "$L" "$status"
      done
    done
  done
}
monitor_progress &
MONITOR_PID=$!

# Launch all 12 tasks
SIM_PIDS=()
for alg in "${ALGOS[@]}"; do
  for L in "${L_VALUES[@]}"; do
    ( ./simulator "$alg" "$L" $EXTRA_ARGS \
        2>"status_${alg}_${L}.tmp" \
      && touch "done_${alg}_${L}.flag" ) &
    SIM_PIDS+=($!)
  done
done
wait "${SIM_PIDS[@]}"

kill "$MONITOR_PID" 2>/dev/null

# Final status: overwrite with all-done display
printf "\033[${N_TASKS}A"
for alg in "${ALGOS[@]}"; do
  for L in "${L_VALUES[@]}"; do
    printf "%-12s L=%4d: [HECHO]%-38s\n" "$alg" "$L" ""
  done
done

# Merge temp files into final output files in L order
echo ""
echo "Combinando resultados..."
> metropolis_2d.txt
> glauber_2d.txt
> wolff_2d.txt
for L in "${L_VALUES[@]}"; do
  cat "tmp_metropolis_${L}.dat" >> metropolis_2d.txt
  cat "tmp_glauber_${L}.dat"    >> glauber_2d.txt
  cat "tmp_wolff_${L}.dat"      >> wolff_2d.txt
  rm  "tmp_metropolis_${L}.dat" \
      "tmp_glauber_${L}.dat"    \
      "tmp_wolff_${L}.dat"
done

rm -f status_*.tmp done_*.flag

ELAPSED=$(( SECONDS - T_START ))
echo "Hecho. ($(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s)"
