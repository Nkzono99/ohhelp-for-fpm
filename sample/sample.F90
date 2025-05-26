!  File: sample.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
module sample
  use ohhelp3

  implicit none
  integer,parameter     :: MAXFRAC=20
  integer,parameter     :: FEB=1,FCD=2
  integer,parameter     :: EX=1,EY=2,EZ=3,BX=4,BY=5,BZ=6
  integer,parameter     :: JX=1,JY=2,JZ=3

  integer               :: sdid(2)
  integer,allocatable   :: nphgram(:,:,:)
  integer,allocatable   :: totalp(:,:)
  type(oh_particle),allocatable&
                        :: pbuf(:)
  integer               :: pbase(3)
  type(oh_mycomm)       :: mycomm
  integer               :: nbor(3,3,3)
  integer,allocatable   :: sdoms(:,:,:)
  integer               :: bcond(2,OH_DIMENSION)
  integer,allocatable   :: bounds(:,:,:)
  integer               :: ftypes(7,3)
  integer               :: cfields(3)
  integer               :: ctypes(3,2,1,2)
  integer               :: fsizes(2,OH_DIMENSION,2)
  real*8,allocatable    :: eb(:,:,:,:,:)
  real*8,allocatable    :: cd(:,:,:,:,:)

  interface
    subroutine initialize_eb(eb, sdom)
      implicit none
      real*8            :: eb(:,:,:,:)
      integer           :: sdom(:,:)
    end subroutine
    subroutine initialize_particles(pbuf, nspec, nphgram)
      use oh_type
      implicit none
      type(oh_particle) :: pbuf(:)
      integer           :: nspec
      integer           :: nphgram(:,:)
    end subroutine
    subroutine lorentz(eb, x, y, z, s, acc)
      implicit none
      real*8            :: eb(:,:,:,:)
      real*8            :: x, y, z
      integer           :: s
      real*8            :: acc(OH_DIMENSION)
    end subroutine
    subroutine scatter(p, s, c)
      use oh_type
      implicit none
      type(oh_particle) :: p
      integer           :: s
      real*8            :: c(3,2,2,2)
    end subroutine
    subroutine rotate_b(eb, x, y, z, rot)
      implicit none
      real*8            :: eb(:,:,:,:)
      integer           :: x, y, z
      real*8            :: rot(OH_DIMENSION)
    end subroutine
    subroutine rotate_e(eb, x, y, z, rot)
      implicit none
      real*8            :: eb(:,:,:,:)
      integer           :: x, y, z
      real*8            :: rot(OH_DIMENSION)
    end subroutine
  end interface

  contains
  subroutine pic(nspec, pcoord, scoord, npmax, nstep)
    implicit none
    integer             :: nspec
    integer             :: pcoord(OH_DIMENSION)
    integer             :: scoord(2,OH_DIMENSION)
    integer*8           :: npmax
    integer             :: nstep

    integer             :: n, t, maxlocalp, currmode

    allocate(totalp(nspec,2))
    n = pcoord(1) * pcoord(2) * pcoord(3)
    allocate(nphgram(n, nspec, 2))
    allocate(sdoms(2, OH_DIMENSION, n))
    allocate(bounds(2, OH_DIMENSION, n))

    maxlocalp = oh_max_local_particles(npmax, MAXFRAC, 0)
    allocate(pbuf(maxlocalp))

    nbor(1,1,1) = -1
    sdoms(1,1,1) = 0;  sdoms(2,1,1) = -1
    bcond(:,:) = reshape((/1,1, 1,1, 1,1/), (/2,OH_DIMENSION/))
    ftypes(:,FEB) = (/6, 0,0, -1,1,  0,0/)                      ! for eb()
    ftypes(:,FCD) = (/3, 0,0,  0,0, -1,2/)                      ! for cd()
    ftypes(1,FCD+1) = -1                                        ! terminator
    cfields(:) = (/FEB,FCD,0/)
    ctypes(:,:,1,FEB) = reshape((/ 0,0,2,  -1,-1,1/), (/3,2/))  ! for eb()
    ctypes(:,:,1,FCD) = reshape((/-1,2,3,  -1,-4,3/), (/3,2/))  ! for cd()

    call oh_init(sdid(:), nspec, MAXFRAC, nphgram(:,:,:), totalp(:,:), &
                 pbuf(:), pbase(:), maxlocalp, mycomm, nbor(:,:,:), &
                 pcoord(:), sdoms(:,:,:), scoord(:,:), 1, bcond(:,:), &
                 bounds(:,:,:), ftypes(:,:), cfields(:), ctypes(:,:,:,:), &
                 fsizes(:,:,:), 0, 0, 0)

    allocate(eb(6, fsizes(1,1,FEB):fsizes(2,1,FEB), &
                   fsizes(1,2,FEB):fsizes(2,2,FEB), &
                   fsizes(1,3,FEB):fsizes(2,3,FEB), 2))
    allocate(cd(3, fsizes(1,1,FCD):fsizes(2,1,FCD), &
                   fsizes(1,2,FCD):fsizes(2,2,FCD), &
                   fsizes(1,3,FCD):fsizes(2,3,FCD), 2))

    call initialize_eb(eb(:,:,:,:,1), sdoms(:,:,sdid(1)))
    call initialize_particles(pbuf(:), nspec, nphgram(:,:,1))

    currmode = oh_transbound(0, 0)
    if (currmode.lt.0) then
      call oh_bcast_field(eb(1,0,0,0,1), eb(1,0,0,0,2), FEB)
      currmode = 1
    end if
    call oh_exchange_borders(eb(1,0,0,0,1), eb(1,0,0,0,2), FEB, currmode)

    do t=1, nstep
      call particle_push(pbuf(pbase(1):), nspec, totalp(:,1), &
                         eb(:,:,:,:,1), sdoms(:,:,sdid(1)), sdid(1), 0, &
                         nphgram(:,:,1))
      if (sdid(2).ge.0) &
        call particle_push(pbuf(pbase(2):), nspec, totalp(:,2), &
                           eb(:,:,:,:,2), sdoms(:,:,sdid(2)), sdid(2), 1, &
                           nphgram(:,:,2))
      currmode = oh_transbound(currmode, 0)
      if (currmode.lt.0) then
        call oh_bcast_field(eb(1,0,0,0,1), eb(1,0,0,0,2), FEB)
        currmode = 1
      end if

      call current_scatter(pbuf(pbase(1):), nspec, totalp(:,1), &
                           cd(:,:,:,:,1), sdoms(:,:,sdid(1)), &
                           ctypes(:,:,1,FCD))
      if (sdid(2).ge.0) &
        call current_scatter(pbuf(pbase(2):), nspec, totalp(:,2), &
                             cd(:,:,:,:,2), sdoms(:,:,sdid(2)), &
                             ctypes(:,:,1,FCD))
      if (currmode.ne.0) &
        call oh_allreduce_field(cd(1,0,0,0,1), cd(1,0,0,0,2), FCD)
      call oh_exchange_borders(cd(1,0,0,0,1), cd(1,0,0,0,2), FCD, currmode)
      call add_boundary_current(cd(:,:,:,:,1), sdoms(:,:,sdid(1)), &
                                ctypes(:,:,1,FCD))
      if (sdid(2).ge.0) &
        call add_boundary_current(cd(:,:,:,:,2), sdoms(:,:,sdid(2)), &
                                  ctypes(:,:,1,FCD))

      call field_solve_e(eb(:,:,:,:,1), cd(:,:,:,:,1), sdoms(:,:,sdid(1)))
      call field_solve_b(eb(:,:,:,:,1), sdoms(:,:,sdid(1)))
      if (sdid(2).ge.0) then
        call field_solve_e(eb(:,:,:,:,2), cd(:,:,:,:,2), sdoms(:,:,sdid(2)))
        call field_solve_b(eb(:,:,:,:,2), sdoms(:,:,sdid(2)))
      end if
      call oh_exchange_borders(eb(1,0,0,0,1), eb(1,0,0,0,2), FEB, currmode)
    end do
  end subroutine

  subroutine particle_push(pbuf, nspec, totalp, eb, sdom, n, ps, nphgram)
    implicit none
    type(oh_particle) :: pbuf(:)
    integer           :: nspec
    integer           :: totalp(:)
    real*8            :: eb(:,:,:,:)
    integer           :: sdom(:,:)
    integer           :: n
    integer           :: ps
    integer           :: nphgram(:,:)

    integer           :: xl, yl, zl, xu, yu, zu
    integer           :: s, p, q, m
    real*8            :: acc(OH_DIMENSION)

    xl=sdom(1,1);  yl=sdom(1,2);  zl=sdom(1,3)
    xu=sdom(2,1);  yu=sdom(2,2);  zu=sdom(2,3)

    p = 0
    do s=1, nspec
      nphgram(n+1,s) = totalp(s)
      do q=1, totalp(s)
        p = p + 1
        call lorentz(eb, pbuf(p)%x-xl, pbuf(p)%y-yl, pbuf(p)%z-zl, s, acc)
        pbuf(p)%vx = pbuf(p)%vx + acc(1)
        pbuf(p)%vy = pbuf(p)%vy + acc(2)
        pbuf(p)%vz = pbuf(p)%vz + acc(3)
        pbuf(p)%x = pbuf(p)%x + pbuf(p)%vx
        pbuf(p)%y = pbuf(p)%y + pbuf(p)%vy
        pbuf(p)%z = pbuf(p)%z + pbuf(p)%vz
        if (pbuf(p)%x.lt.xl .or. pbuf(p)%x.ge.xu .or. &
            pbuf(p)%y.lt.yl .or. pbuf(p)%y.ge.yu .or. &
            pbuf(p)%z.lt.zl .or. pbuf(p)%z.ge.zu) then
          m = oh_map_particle_to_neighbor(pbuf(p)%x, pbuf(p)%y, pbuf(p)%z, ps)
          nphgram(n+1,s) = nphgram(n+1,s) - 1
          nphgram(m+1,s) = nphgram(m+1,s) + 1
          pbuf(p)%nid = m
        end if
      end do
    end do
  end subroutine

  subroutine current_scatter(pbuf, nspec, totalp, cd, sdom, ctype)
    implicit none
    type(oh_particle) :: pbuf(:)
    integer           :: nspec
    integer           :: totalp(:)
    real*8            :: cd(:,:,:,:)
    integer           :: sdom(:,:)
    integer           :: ctype(3,2)

    integer           :: xl, yl, zl, xu, yu, zu
    integer           :: s, p, q
    integer           :: i, j, k
    real*8            :: x, y, z
    real*8            :: c(3,2,2,2)

    xl = sdom(1,1);  yl = sdom(1,2);  zl = sdom(1,3)
    xu = sdom(2,1)-xl;  yu = sdom(2,2)-yl;  zu = sdom(2,3)-zl
    do k=ctype(1,1), zu+ctype(1,2)+ctype(1,3)-1
      do j=ctype(1,1), yu+ctype(1,2)+ctype(1,3)-1
        do i=ctype(1,1), xu+ctype(1,2)+ctype(1,3)-1
          cd(JX, i, j, k) = 0.0d0
          cd(JY, i, j, k) = 0.0d0
          cd(JZ, i, j, k) = 0.0d0
    end do;  end do;  end do

    p = 0
    do s=1, nspec
      do q=1, totalp(s)
        p = p + 1
        call scatter(pbuf(p), s, c)
        x = pbuf(p)%x - xl;  y = pbuf(p)%y - yl;  z = pbuf(p)%z - zl
        do k=0,1;  do j=0,1;  do i=0,1
          cd(JX, x+i, y+j, z+k) = cd(JX, x+i, y+j, z+k) + c(JX, i, j, k)
          cd(JY, x+i, y+j, z+k) = cd(JY, x+i, y+j, z+k) + c(JY, i, j, k)
          cd(JZ, x+i, y+j, z+k) = cd(JZ, x+i, y+j, z+k) + c(JZ, i, j, k)
        end do;  end do;  end do;
      end do
    end do
  end subroutine

  subroutine add_boundary_current(cd, sdom, ctype)
    implicit none
    real*8            :: cd(:,:,:,:)
    integer           :: sdom(2,OH_DIMENSION)
    integer           :: ctype(3,2)

    integer           :: xu, yu, zu
    integer           :: sl, dl, nl, su, du, nu

    xu = sdom(2,1) - sdom(1,1)
    yu = sdom(2,2) - sdom(1,2)
    zu = sdom(2,3) - sdom(1,3)
    sl = ctype(2,2);  nl = ctype(3,2);  dl = sl + nl
    su = ctype(2,1);  nu = ctype(3,1);  du = su - nu

    call add_boundary_curr(sl, sl, xu+(su+nu-sl), &
                           sl, sl, yu+(su+nu-sl), &
                           sl, dl, nl, cd)
    call add_boundary_curr(sl, sl, xu+(su+nu-sl), &
                           sl, sl, yu+(su+nu-sl), &
                           zu+su, zu+du, nu, cd)
    call add_boundary_curr(sl, sl, xu+(su+nu-sl), &
                           sl, dl, nl, &
                           dl, dl, zu+(du-dl), cd)
    call add_boundary_curr(sl, sl, xu+(su+nu-sl), &
                           yu+su, yu+du, nu, &
                           dl, dl, zu+(du-dl), cd)
    call add_boundary_curr(sl, dl, nl, &
                           dl, dl, yu+(du-dl), &
                           dl, dl, zu+(du-dl), cd)
    call add_boundary_curr(xu+su, xu+du, nu, &
                           dl, dl, yu+(du-dl), &
                           dl, dl, zu+(du-dl), cd)
  end subroutine

  subroutine add_boundary_curr(xs, xd, nx, ys, yd, ny, zs, zd, nz, cd)
    implicit none
    integer           :: xs, xd, nx, ys, yd, ny, zs, zd, nz
    integer           :: i, j, k
    real*8            :: cd(:,:,:,:)

    do k=0, nz-1;  do j=0, ny-1;  do i=0, nx-1
      cd(JX, xd+i, yd+j, zd+k) = &
           cd(JX, xd+i, yd+j, zd+k) + cd(JX, xs+i, ys+j, zs+k)
      cd(JY, xd+i, yd+j, zd+k) = &
           cd(JY, xd+i, yd+j, zd+k) + cd(JY, xs+i, ys+j, zs+k)
      cd(JZ, xd+i, yd+j, zd+k) = &
           cd(JZ, xd+i, yd+j, zd+k) + cd(JZ, xs+i, ys+j, zs+k)
    end do;  end do;  end do
  end subroutine

  subroutine field_solve_e(eb, cd, sdom)
    implicit none
    real*8            :: eb(:,:,:,:)
    real*8            :: cd(:,:,:,:)
    integer           :: sdom(2,OH_DIMENSION)

    integer           :: xu, yu, zu, x, y, z
    real*8            :: rot(OH_DIMENSION)

    xu = sdom(2,1) - sdom(1,1)
    yu = sdom(2,2) - sdom(1,2)
    zu = sdom(2,3) - sdom(1,3)
    do z=0, zu;  do y=0, yu;  do x=0, xu
      call rotate_b(eb(:,:,:,:), x, y, z, rot)
      eb(EX, x, y, z) = eb(EX, x, y, z) + &
                        (1/EPSILON)*((1/MU)*rot(1) + cd(JX, x, y, z))
      eb(EY, x, y, z) = eb(EY, x, y, z) + &
                        (1/EPSILON)*((1/MU)*rot(2) + cd(JY, x, y, z))
      eb(EZ, x, y, z) = eb(EZ, x, y, z) + &
                        (1/EPSILON)*((1/MU)*rot(3) + cd(JZ, x, y, z))
    end do;  end do;  end do
  end subroutine

  subroutine field_solve_b(eb, sdom)
    implicit none
    real*8            :: eb(:,:,:,:)
    integer           :: sdom(2,OH_DIMENSION)

    integer           :: xu, yu, zu, x, y, z
    real*8            :: rot(OH_DIMENSION)

    xu = sdom(2,1) - sdom(1,1)
    yu = sdom(2,2) - sdom(1,2)
    zu = sdom(2,3) - sdom(1,3)
    do z=0, zu-1;  do y=0, yu-1;  do x=0, xu-1
      call rotate_e(eb(:,:,:,:), x, y, z, rot)
      eb(BX, x, y, z) = eb(BX, x, y, z) + rot(1)
      eb(BY, x, y, z) = eb(BY, x, y, z) + rot(2)
      eb(BZ, x, y, z) = eb(BZ, x, y, z) + rot(3)
    end do;  end do;  end do
  end subroutine
end module
