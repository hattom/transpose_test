# test_transpose
Testing array transpose in modern fortran. To be used as a proxy for the pre-FFT local data transpose in the gyrokinetic code [ORB5](https://www.epfl.ch/research/domains/swiss-plasma-center/research/theory/codes/research_theory_codes_orb5/).


## Building
* Intel: `mpiifort -DCACHE_BLOCK_SIZE=32 -DMKL -xHOST -O2 -fopenmp test_loop.F90 -L${MKLROOT}/lib/intel64 -lmkl_rt -Wl,-rpath,${MKLROOT}/lib/intel64`
* Remove `-DMKL` if not linking against MKL

TODO: Makefile

## Testing
```
partition=${...}
qos=${...}
for OMP_NUM_THREADS in 1 2 4 8 16 32; do ((export OMP_NUM_THREADS; srun -l -p ${partition} -q ${qos} -t 00:05:00 -N1 -n $((32/OMP_NUM_THREADS)) -c ${OMP_NUM_THREADS} ./a.out | tee logs.${OMP_NUM_THREADS}.out) &); done
```
TODO: proper pinning
