!  File: oh_mod4p.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module ohhelp4p
  use ohhelp3
  implicit none
  interface
    subroutine oh4p_init(sdid, nspec, maxfrac, totalp, pbuf, pbase, &
                         maxlocalp, mycomm, nbor, pcoord, sdoms, scoord, &
                         nbound, bcond, bounds, ftypes, cfields, ctypes, &
                         fsizes, &
                         stats, repiter, verbose)
      use oh_type
      implicit none
      integer,intent(out)   :: sdid(2)
      integer,intent(in)    :: nspec
      integer,intent(in)    :: maxfrac
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
    integer function oh4p_max_local_particles(npmax, maxfrac, minmargin, &
                                              hsthresh)
      implicit none
      integer*8,intent(in) :: npmax
      integer,intent(in)   :: maxfrac
      integer,intent(in)   :: minmargin
      integer,intent(in)   :: hsthresh
    end function
    subroutine oh4p_per_grid_histogram(pghgram)
      implicit none
      integer,intent(inout) :: pghgram
    end subroutine
    integer function oh4p_transbound(currmode, stats)
      implicit none
      integer,intent(in) :: currmode
      integer,intent(in) :: stats
    end function
    integer function oh4p_map_particle_to_neighbor(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
    integer function oh4p_map_particle_to_subdomain(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
    integer function oh4p_inject_particle(part, ps)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
    end function
    subroutine oh4p_remove_mapped_particle(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end subroutine
    integer function oh4p_remap_particle_to_neighbor(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
    integer function oh4p_remap_particle_to_subdomain(part, ps, s)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
      integer,intent(in)   :: ps
      integer,intent(in)   :: s
    end function
  end interface
end module
