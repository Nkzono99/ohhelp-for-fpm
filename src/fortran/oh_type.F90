!  File: oh_type.F90
!  Version 1.1.1 (2015/10/23)
!  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
!                           (ACCMS, Kyoto University)
!  This program can be freely used, redistributed and modified for non-
!  commercial purpose providing that the copyright notice above remains
!  unchanged.
#include "oh_config.h"
module oh_type
    use, intrinsic :: iso_c_binding
    implicit none

    private
    public oh_mycomm
    public oh_particle
    public OH_PCL_ALIVE
    public OH_PCL_INJECTED
    public OH_PCL_TO_BE_ACCUMULATED
    public OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL

    type, bind(c) :: oh_mycomm
        integer(c_int) :: prime, sec, rank, root, black
    end type

    type, bind(c) :: oh_particle
        real(c_double) :: x, y, z, vx, vy, vz
        integer(c_int)   :: pid

#ifdef OH_BIG_SPACE
        integer(c_long) :: nid
#else
        integer(c_int) :: nid
#endif

        integer(c_int) :: spec
        integer(c_int) :: preside
    end type

    !> preside
    integer, parameter :: OH_PCL_ALIVE = 0
    integer, parameter :: OH_PCL_INJECTED = 1
    integer, parameter :: OH_PCL_TO_BE_ACCUMULATED = -1
    integer, parameter :: OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL = -2

end module
