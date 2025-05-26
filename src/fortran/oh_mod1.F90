!  File: oh_mod1.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module ohhelp1
  use oh_type
  implicit none
  interface
    subroutine oh1_fam_comm(fortran_comm)
      implicit none
      integer,intent(in)    :: fortran_comm
    end subroutine
    subroutine oh1_init(sdid, nspec, maxfrac, nphgram, totalp, rcounts, &
                        scounts, mycomm, nbor, pcoord, stats, repiter, verbose)
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
      integer,intent(in)    :: stats
      integer,intent(in)    :: repiter
      integer,intent(in)    :: verbose
    end subroutine
    subroutine oh1_neighbors(nbor)
      implicit none
#if OH_DIMENSION==1
      integer,intent(inout) :: nbor(3,3)
#else
#if OH_DIMENSION==2
      integer,intent(inout) :: nbor(3,3,3)
#else
      integer,intent(inout) :: nbor(3,3,3,3)
#endif
#endif
    end subroutine
    subroutine oh1_families(famindex, members)
      implicit none
      integer,intent(inout) :: famindex(:)
      integer,intent(inout) :: members(:)
    end subroutine
    integer function oh1_transbound(currmode, stats)
      implicit none
      integer,intent(in) :: currmode
      integer,intent(in) :: stats
    end function
    integer function oh1_accom_mode()
      implicit none
    end function
    subroutine oh1_broadcast(pbuf, sbuf, pcount, scount, ptype, stype)
      implicit none
      real*8,intent(in)  :: pbuf
      real*8,intent(out) :: sbuf
      integer,intent(in) :: pcount
      integer,intent(in) :: scount
      integer,intent(in) :: ptype
      integer,intent(in) :: stype
    end subroutine
    subroutine oh1_all_reduce(pbuf, sbuf, pcount, scount, ptype, stype, &
                              pop, sop)
      implicit none
      real*8,intent(inout) :: pbuf
      real*8,intent(inout) :: sbuf
      integer,intent(in)   :: pcount
      integer,intent(in)   :: scount
      integer,intent(in)   :: ptype
      integer,intent(in)   :: stype
      integer,intent(in)   :: pop
      integer,intent(in)   :: sop
    end subroutine
    subroutine oh1_reduce(pbuf, sbuf, pcount, scount, ptype, stype, pop, sop)
      implicit none
      real*8,intent(inout) :: pbuf
      real*8,intent(in)    :: sbuf
      integer,intent(in)   :: pcount
      integer,intent(in)   :: scount
      integer,intent(in)   :: ptype
      integer,intent(in)   :: stype
      integer,intent(in)   :: pop
      integer,intent(in)   :: sop
    end subroutine
    subroutine oh1_init_stats(key, ps)
      implicit none
      integer,intent(in) :: key
      integer,intent(in) :: ps
    end subroutine
    subroutine oh1_stats_time(key, ps)
      implicit none
      integer,intent(in) :: key
      integer,intent(in) :: ps
    end subroutine
    subroutine oh1_show_stats(step, currmode)
      implicit none
      integer,intent(in) :: step
      integer,intent(in) :: currmode
    end subroutine
    subroutine oh1_print_stats(nstep)
      implicit none
      integer,intent(in) :: nstep
    end subroutine
    subroutine oh1_verbose(message)
      implicit none
      character(*),intent(in) :: message
    end subroutine
  end interface
end module
