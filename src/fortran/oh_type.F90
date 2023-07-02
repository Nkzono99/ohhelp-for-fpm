!  File: oh_type.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module oh_type
    type oh_mycomm
        sequence
        integer :: prime, sec, rank, root, black
    end type oh_mycomm

    type oh_particle
        sequence
        real*8 :: x, y, z, vx, vy, vz
        integer   :: pid
        integer :: preside

#ifdef OH_BIG_SPACE
        integer*8 :: nid
#else
        integer :: nid
#endif

        integer :: spec
    end type

    !> preside
    integer, parameter :: OH_PCL_ALIVE = 0
    integer, parameter :: OH_PCL_INJECTED = 1
    integer, parameter :: OH_PCL_TO_BE_ACCUMULATED = -1
    integer, parameter :: OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL = -2
end module
