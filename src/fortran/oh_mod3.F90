!  File: oh_mod3.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module ohhelp3
  use ohhelp2
  implicit none
  interface
    subroutine oh3_init(sdid, nspec, maxfrac, nphgram, totalp, &
                        pbuf, pbase, maxlocalp, mycomm, nbor, pcoord, &
                        sdoms, scoord, nbound, bcond, bounds, ftypes, &
                        cfields, ctypes, fsizes, &
                        stats, repiter, verbose)
      use oh_type
      implicit none
      integer,intent(out)   :: sdid(2)
      integer,intent(in)    :: nspec
      integer,intent(in)    :: maxfrac
      integer,intent(inout) :: nphgram(:,:,:)
      integer,intent(out)   :: totalp(:,:)
      type(oh_particle),intent(inout) :: pbuf(:)
      integer,intent(out)   :: pbase(3)
      integer,intent(in)    :: maxlocalp
      type(oh_mycomm),intent(out) :: mycomm
#if OH_DIMENSION==1
      integer,intent(inout) :: nbor(3)
#else
#if OH_DIMENSION==2
      integer,intent(inout) :: nbor(3,3)
#else
      integer,intent(inout) :: nbor(3,3,3)
#endif
#endif
      integer,intent(in)    :: pcoord(OH_DIMENSION)
      integer,intent(inout) :: sdoms(:,:,:)
      integer,intent(in)    :: scoord(2,OH_DIMENSION)
      integer,intent(in)    :: nbound
      integer,intent(in)    :: bcond(2,OH_DIMENSION)
      integer,intent(inout) :: bounds(:,:,:)
      integer,intent(in)    :: ftypes(:,:)
      integer,intent(in)    :: cfields(:)
      integer,intent(in)    :: ctypes(:,:,:,:)
      integer,intent(out)   :: fsizes(:,:,:)
      integer,intent(in)    :: stats
      integer,intent(in)    :: repiter
      integer,intent(in)    :: verbose
    end subroutine
    subroutine oh13_init(sdid, nspec, maxfrac, nphgram, totalp, &
                         rcounts, scounts, mycomm, nbor, pcoord, &
                         sdoms, scoord, nbound, bcond, bounds, ftypes, &
                         cfields, ctypes, fsizes, &
                         stats, repiter, verbose)
      use oh_type
      implicit none
      integer,intent(out)   :: sdid(2)
      integer,intent(in)    :: nspec
      integer,intent(in)    :: maxfrac
      integer,intent(inout) :: nphgram(:,:,:)
      integer,intent(out)   :: totalp(:,:)
      integer,intent(out)   :: rcounts(:,:,:)
      integer,intent(out)   :: scounts(:,:,:)
      type(oh_mycomm),intent(out) :: mycomm
#if OH_DIMENSION==1
      integer,intent(inout) :: nbor(3)
#else
#if OH_DIMENSION==2
      integer,intent(inout) :: nbor(3,3)
#else
      integer,intent(inout) :: nbor(3,3,3)
#endif
#endif
      integer,intent(in)    :: pcoord(OH_DIMENSION)
      integer,intent(inout) :: sdoms(:,:,:)
      integer,intent(in)    :: scoord(2,OH_DIMENSION)
      integer,intent(in)    :: nbound
      integer,intent(in)    :: bcond(2,OH_DIMENSION)
      integer,intent(inout) :: bounds(:,:,:)
      integer,intent(in)    :: ftypes(:,:)
      integer,intent(in)    :: cfields(:)
      integer,intent(in)    :: ctypes(:,:,:,:)
      integer,intent(out)   :: fsizes(:,:,:)
      integer,intent(in)    :: stats
      integer,intent(in)    :: repiter
      integer,intent(in)    :: verbose
    end subroutine
    subroutine oh3_grid_size(size)
      implicit none
      real*8,intent(in)     :: size(OH_DIMENSION)
    end subroutine
    integer function oh3_transbound(currmode, stats)
      implicit none
      integer,intent(in) :: currmode
      integer,intent(in) :: stats
    end function
#if OH_DIMENSION==1
    integer function oh3_map_particle_to_neighbor(x, ps)
      implicit none
      real*8,intent(inout) :: x
      integer,intent(in)   :: ps
    end function
#else
#if OH_DIMENSION==2
    integer function oh3_map_particle_to_neighbor(x, y, ps)
      implicit none
      real*8,intent(inout) :: x
      real*8,intent(inout) :: y
      integer,intent(in)   :: ps
    end function
#else
    integer function oh3_map_particle_to_neighbor(x, y, z, ps)
      implicit none
      real*8,intent(inout) :: x
      real*8,intent(inout) :: y
      real*8,intent(inout) :: z
      integer,intent(in)   :: ps
    end function
#endif
#endif
#if OH_DIMENSION==1
    integer function oh3_map_particle_to_subdomain(x)
      implicit none
      real*8,intent(in)  :: x
    end function
#else
#if OH_DIMENSION==2
    integer function oh3_map_particle_to_subdomain(x, y)
      implicit none
      real*8,intent(in)  :: x
      real*8,intent(in)  :: y
    end function
#else
    integer function oh3_map_particle_to_subdomain(x, y, z)
      implicit none
      real*8,intent(in)  :: x
      real*8,intent(in)  :: y
      real*8,intent(in)  :: z
    end function
#endif
#endif
    subroutine oh3_bcast_field(pfld, sfld, ftype)
      implicit none
      real*8,intent(in)  :: pfld
      real*8,intent(out) :: sfld
      integer,intent(in) :: ftype
    end subroutine
    subroutine oh3_allreduce_field(pfld, sfld, ftype)
      implicit none
      real*8,intent(inout) :: pfld
      real*8,intent(inout) :: sfld
      integer,intent(in)   :: ftype
    end subroutine
    subroutine oh3_reduce_field(pfld, sfld, ftype)
      implicit none
      real*8,intent(inout) :: pfld
      real*8,intent(in)    :: sfld
      integer,intent(in)   :: ftype
    end subroutine
    subroutine oh3_exchange_borders(pfld, sfld, ctype, bcast)
      implicit none
      real*8,intent(inout) :: pfld
      real*8,intent(out)   :: sfld
      integer,intent(in)   :: ctype
      integer,intent(in)   :: bcast
    end subroutine
  end interface
end module
