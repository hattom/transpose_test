#ifndef CACHE_BLOCK_SIZE
#define CACHE_BLOCK_SIZE 32
#endif

program test_transpose
  use iso_fortran_env, only: REAL64
  implicit none

  real(kind=REAL64), dimension(:,:,:), allocatable :: rhs_bspl
  real, dimension(:,:,:), allocatable :: r2c_rhs

  integer :: ns, nchi, nphi, nidbas, nvp_cart, nphip
  integer :: islw, isup, jclw, jcup, kplw, kpup
  integer :: iter, niter
  integer :: i, j, k, ib, jb

  ns = 1024
  nchi = 1024
  nphi = 1
  nidbas = 3
  niter = 100

  nvp_cart = nphi/1
  nphip = nphi / nvp_cart + nidbas + 0

  if(nidbas == 3) then
    islw  = -1; isup  = ns+1
  else
    islw  = -1; isup  = ns+1
  endif
  jclw = 0; jcup = nchi-1
  kplw = 0; kpup = nphip-1

  allocate(rhs_bspl(islw:isup, jclw:jcup, kplw:kpup))
  allocate(r2c_rhs(jclw:jcup, islw:isup, kplw:kpup))

  rhs_bspl(:,:,:) = 0.0
  r2c_rhs(:,:,:) = 0.0

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    do k=kplw, kpup
      r2c_rhs(jclw:jcup, islw:isup, k) = transpose(rhs_bspl(islw:isup, jclw:jcup, k))
    enddo
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'transpose time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    !$omp parallel do
    do k=kplw, kpup
      r2c_rhs(jclw:jcup, islw:isup, k) = transpose(rhs_bspl(islw:isup, jclw:jcup, k))
    enddo
    !$omp end parallel do
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'transpose omp time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    DO k = kplw, kpup
      DO jb = jclw, jcup, CACHE_BLOCK_SIZE
        DO ib = islw, isup, CACHE_BLOCK_SIZE
          do j=jb, min(jcup, jb+CACHE_BLOCK_SIZE)
            do i=ib, min(isup, ib+CACHE_BLOCK_SIZE)
              r2c_rhs(j, i, k) = rhs_bspl(i, j, k)
          end do
        end do
      END DO
    END DO
  END DO
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'blocked time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    !$omp parallel do
    DO k = kplw, kpup
      DO jb = jclw, jcup, CACHE_BLOCK_SIZE
        DO ib = islw, isup, CACHE_BLOCK_SIZE
          do j=jb, min(jcup, jb+CACHE_BLOCK_SIZE)
            do i=ib, min(isup, ib+CACHE_BLOCK_SIZE)
              r2c_rhs(j, i, k) = rhs_bspl(i, j, k)
            end do
          end do
        END DO
      END DO
    END DO
    !$omp end parallel do
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'blocked omp time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    !$omp parallel do collapse(2)
    DO k = kplw, kpup
      DO jb = jclw, jcup, CACHE_BLOCK_SIZE
        DO ib = islw, isup, CACHE_BLOCK_SIZE
          do j=jb, min(jcup, jb+CACHE_BLOCK_SIZE)
            do i=ib, min(isup, ib+CACHE_BLOCK_SIZE)
              r2c_rhs(j, i, k) = rhs_bspl(i, j, k)
            end do
          end do
        END DO
      END DO
    END DO
    !$omp end parallel do
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'blocked omp (c2) time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    !$omp parallel do collapse(3)
    DO k = kplw, kpup
      DO jb = jclw, jcup, CACHE_BLOCK_SIZE
        DO ib = islw, isup, CACHE_BLOCK_SIZE
          do j=jb, min(jcup, jb+CACHE_BLOCK_SIZE)
            do i=ib, min(isup, ib+CACHE_BLOCK_SIZE)
              r2c_rhs(j, i, k) = rhs_bspl(i, j, k)
            end do
          end do
        END DO
      END DO
    END DO
    !$omp end parallel do
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'blocked omp (c3) time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    DO k = kplw, kpup
      DO j = jclw, jcup
        DO i = islw, isup
          r2c_rhs(j, i, k) = rhs_bspl(i, j, k)
      END DO
    END DO
  END DO
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'loop time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    !$omp parallel do
    DO k = kplw, kpup
      DO j = jclw, jcup
        DO i = islw, isup
          r2c_rhs(j, i, k) = rhs_bspl(i, j, k)
        END DO
      END DO
    END DO
    !$omp end parallel do
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'loop omp time', t1-t0
  end block

#ifdef MKL
  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    do k=kplw, kpup
      call mkl_domatcopy('C', 'T', (jcup-jclw+1), (isup-islw+1), 1., &
            r2c_rhs(jclw:jcup, islw:isup, k), (jcup-jclw+1), &
            rhs_bspl(islw:isup, jclw:jcup, k), (isup-islw+1) )
    enddo
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'mkl time', t1-t0
  end block

  block
  !$ use omp_lib
  real(kind=REAL64) :: t1, t0
  !$ t0 = omp_get_wtime()
  do iter = 1,niter
    if(mod(iter, 100) == -1) print *, iter, '/', niter
    !$omp parallel do
    do k=kplw, kpup
      call mkl_domatcopy('C', 'T', (jcup-jclw+1), (isup-islw+1), 1., &
            r2c_rhs(jclw:jcup, islw:isup, k), (jcup-jclw+1), &
            rhs_bspl(islw:isup, jclw:jcup, k), (isup-islw+1) )
    enddo
    !$omp end parallel do
  end do
  !$ t1 = omp_get_wtime()
  !$ print *, 'mkl omp time', t1-t0
  end block
#endif
end program test_transpose
