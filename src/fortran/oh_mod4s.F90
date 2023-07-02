!  File: oh_mod4s.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module ohhelp4s
  use ohhelp3
  implicit none
  interface
    subroutine oh4s_init(sdid, nspec, maxfrac, npmax, minmargin, maxdensity, &
                         totalp, pbase, maxlocalp, cbufsize, mycomm, nbor, &
                         pcoord, sdoms, scoord, nbound, bcond, bounds, &
                         ftypes, cfields, ctypes, fsizes, zbound, &
                         stats, repiter, verbose)
      use oh_type
      implicit none
      integer,intent(out)   :: sdid(2)
      integer,intent(in)    :: nspec
      integer,intent(in)    :: maxfrac
      integer*8,intent(in)  :: npmax
      integer,intent(in)    :: minmargin
      integer,intent(in)    :: maxdensity
      integer,intent(out)   :: totalp(:,:)
      integer,intent(out)   :: pbase(3)
      integer,intent(out)   :: maxlocalp
      integer,intent(out)   :: cbfsize
      type(oh_mycomm),intent(out) :: mycomm
      integer,intent(inout) :: nbor(3,3,3)
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
      integer,intent(out)   :: zbound(2,2)
      integer,intent(in)    :: stats
      integer,intent(in)    :: repiter
      integer,intent(in)    :: verbose
    end subroutine
    integer subroutine oh4s_particle_buffer(maxlocalp, pbuf)
      use oh_type
      implicit none
      integer,intent(in)   :: maxlocalp
      type(oh_particle),intent(inout) :: pbuf(:)
    end subroutine
    subroutine oh4s_per_grid_histogram(pghgram, pgindex)
      implicit none
      integer,intent(inout) :: pghgram
      integer,intent(inout) :: pgindex
    end subroutine
    integer function oh4s_transbound(currmode, stats)
      implicit none
      integer,intent(in) :: currmode
      integer,intent(in) :: stats
    end function
    subroutine oh4s_exchange_border_data(buf, sbuf, rbuf, type)
      implicit none
      real*8,intent(inout) :: buf
      real*8,intent(out)   :: sbuf
      real*8,intent(out)   :: rbuf
      integer,intent(in)   :: type
    end subroutine
    integer function oh4s_map_particle_to_neighbor(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
    integer function oh4s_map_particle_to_subdomain(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
    integer function oh4s_inject_particle(part, ps)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
    end function
    subroutine oh4s_remove_mapped_particle(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end subroutine
    integer function oh4s_remap_particle_to_neighbor(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
    integer function oh4s_remap_particle_to_subdomain(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
  end interface
end module
