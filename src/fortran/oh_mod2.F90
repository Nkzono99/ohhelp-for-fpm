!  File: oh_mod2.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module ohhelp2
  use ohhelp1
  implicit none
  interface
    subroutine oh2_init(sdid, nspec, maxfrac, nphgram, totalp, &
                        pbuf, pbase, maxlocalp, mycomm, nbor, pcoord, &
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
      integer,intent(in)    :: stats
      integer,intent(in)    :: repiter
      integer,intent(in)    :: verbose
    end subroutine
    integer function oh2_max_local_particles(npmax, maxfrac, minmargin)
      implicit none
      integer*8,intent(in) :: npmax
      integer,intent(in)   :: maxfrac
      integer,intent(in)   :: minmargin
    end function
    integer function oh2_transbound(currmode, stats)
      implicit none
      integer,intent(in) :: currmode
      integer,intent(in) :: stats
    end function
    subroutine oh2_inject_particle(part)
      use oh_type
      implicit none
      type(oh_particle),intent(in) :: part
    end subroutine
    subroutine oh2_remap_injected_particle(part)
      use oh_type
      implicit none
      type(oh_particle),intent(in) :: part
    end subroutine
    subroutine oh2_remove_injected_particle(part)
      use oh_type
      implicit none
      type(oh_particle),intent(inout) :: part
    end subroutine
    subroutine oh1_set_total_particles
    end subroutine
  end interface
end module
