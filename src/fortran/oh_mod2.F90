!  File: oh_mod2.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module ohhelp2
    use oh_type
    implicit none

    private
    public oh2_init
    public oh2_max_local_particles
    public oh2_transbound
    public oh2_inject_particle
    public oh2_remap_injected_particle
    public oh2_remove_injected_particle
    public oh1_set_total_particles

    interface
        subroutine oh2_init_(sdid, nspec, maxfrac, nphgram, totalp, &
                             pbuf, pbase, maxlocalp, mycomm, nbor, pcoord, &
                             stats, repiter, verbose) bind(c, name='oh2_init_')
            use, intrinsic :: iso_c_binding
            import oh_particle
            import oh_mycomm
            implicit none
            integer(c_int), intent(out)   :: sdid
            integer(c_int), intent(in)    :: nspec
            integer(c_int), intent(in)    :: maxfrac
            integer(c_int), intent(inout) :: nphgram
            integer(c_int), intent(out)   :: totalp
            type(oh_particle), intent(inout) :: pbuf
            integer(c_int), intent(out)   :: pbase
            integer(c_int), intent(in)    :: maxlocalp
            type(oh_mycomm), intent(out) :: mycomm
            integer(c_int), intent(inout) :: nbor
            integer(c_int), intent(in)    :: pcoord
            integer(c_int), intent(in)    :: stats
            integer(c_int), intent(in)    :: repiter
            integer(c_int), intent(in)    :: verbose
        end subroutine

        integer function oh2_max_local_particles(npmax, maxfrac, minmargin)
            implicit none
            integer*8, intent(in) :: npmax
            integer, intent(in)   :: maxfrac
            integer, intent(in)   :: minmargin
        end function

        integer function oh2_transbound(currmode, stats)
            implicit none
            integer, intent(in) :: currmode
            integer, intent(in) :: stats
        end function

        subroutine oh2_inject_particle(part)
            import oh_particle
            implicit none
            type(oh_particle), intent(in) :: part
        end subroutine

        subroutine oh2_remap_injected_particle(part)
            import oh_particle
            implicit none
            type(oh_particle), intent(in) :: part
        end subroutine

        subroutine oh2_remove_injected_particle(part)
            import oh_particle
            implicit none
            type(oh_particle), intent(inout) :: part
        end subroutine

        subroutine oh1_set_total_particles
        end subroutine
    end interface

contains

    subroutine oh2_init(sdid, nspec, maxfrac, nphgram, totalp, &
                        pbuf, pbase, maxlocalp, mycomm, nbor, pcoord, &
                        stats, repiter, verbose)
        integer, intent(out)   :: sdid(2)
        integer, intent(in)    :: nspec
        integer, intent(in)    :: maxfrac
        integer, intent(inout) :: nphgram(:, :, :)
        integer, intent(out)   :: totalp(:, :)
        type(oh_particle), intent(inout) :: pbuf(:)
        integer, intent(out)   :: pbase(3)
        integer, intent(in)    :: maxlocalp
        type(oh_mycomm), intent(out) :: mycomm

#if OH_DIMENSION==1
        integer, intent(inout) :: nbor(3)
#else
#if OH_DIMENSION==2
        integer, intent(inout) :: nbor(3, 3)
#else
        integer, intent(inout) :: nbor(3, 3, 3)
#endif
#endif

        integer, intent(in)    :: pcoord(OH_DIMENSION)
        integer, intent(in)    :: stats
        integer, intent(in)    :: repiter
        integer, intent(in)    :: verbose

        call oh2_init_(sdid(1), nspec, maxfrac, nphgram(1, 1, 1), totalp(1, 1), &
                       pbuf(1), pbase(1), maxlocalp, mycomm, &
#if OH_DIMENSION==1
                       nbor(1), &
#else
#if OH_DIMENSION==2
                       nbor(1, 1), &
#else
                       nbor(1, 1, 1), &
#endif
#endif
                       pcoord(1), &
                       stats, repiter, verbose)
    end subroutine

end module
