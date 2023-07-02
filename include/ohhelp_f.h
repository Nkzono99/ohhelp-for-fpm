/* File: ohhelp_f.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#include "oh_config.h"

#define oh_neighbors(A1) \
        oh1_neighbors(A1)
#define oh_families(A1,A2) \
        oh1_families(A1,A2)
#define oh_accom_mode() \
        oh1_accom_mode()
#define oh_broadcast(A1,A2,A3,A4,A5,A6) \
        oh1_broadcast(A1,A2,A3,A4,A5,A6)
#define oh_all_reduce(A1,A2,A3,A4,A5,A6,A7,A8) \
        oh1_all_reduce(A1,A2,A3,A4,A5,A6,A7,A8)
#define oh_reduce(A1,A2,A3,A4,A5,A6,A7,A8) \
        oh1_reduce(A1,A2,A3,A4,A5,A6,A7,A8)
#define oh_init_stats(A1,A2)    oh1_init_stats(A1,A2)
#define oh_stats_time(A1,A2)    oh1_stats_time(A1,A2)
#define oh_show_stats(A1,A2)    oh1_show_stats(A1,A2)
#define oh_print_stats(A1)      oh1_print_stats(A1)
#define oh_verbose(A1)          oh1_verbose(A1)
#if OH_LIB_LEVEL==1
#define oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13) \
        oh1_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13)
#define oh_transbound(A1,A2) oh1_transbound(A1,A2)
#else
#define oh_set_total_particles oh2_set_total_particles
#if OH_LIB_LEVEL!=4
#define oh_max_local_particles(A1,A2,A3) oh2_max_local_particles(A1,A2,A3)
#define oh_inject_particle(A1)           oh2_inject_particle(A1)
#define oh_remap_injected_particle(A1)   oh2_remap_injected_particle(A1)
#define oh_remove_injected_particle(A1)  oh2_remove_injected_particle(A1)
#endif
#if OH_LIB_LEVEL==2
#define oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14) \
        oh2_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)
#define oh_transbound(A1,A2) oh2_transbound(A1,A2)
#else
#define oh_grid_size(A1)                 oh3_grid_size(A1)
#define oh_bcast_field(A1,A2,A3)         oh3_bcast_field(A1,A2,A3)
#define oh_reduce_field(A1,A2,A3)        oh3_reduce_field(A1,A2,A3)
#define oh_allreduce_field(A1,A2,A3)     oh3_allreduce_field(A1,A2,A3)
#define oh_exchange_borders(A1,A2,A3,A4) oh3_exchange_borders(A1,A2,A3,A4)
#if OH_LIB_LEVEL==3
#if OH_DIMENSION==1
#define oh_map_particle_to_neighbor(A1,A2) \
        oh3_map_particle_to_neighbor(A1,A2)
#define oh_map_particle_to_subdomain(A1) \
        oh3_map_particle_to_subdomain(A1)
#elif OH_DIMENSION==2
#define oh_map_particle_to_neighbor(A1,A2,A3) \
        oh3_map_particle_to_neighbor(A1,A2,A3)
#define oh_map_particle_to_subdomain(A1,A2) \
        oh3_map_particle_to_subdomain(A1,A2)
#else
#define oh_map_particle_to_neighbor(A1,A2,A3,A4) \
        oh3_map_particle_to_neighbor(A1,A2,A3,A4)
#define oh_map_particle_to_subdomain(A1,A2,A3) \
        oh3_map_particle_to_subdomain(A1,A2,A3)
#endif
#define \
oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23) \
oh3_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23)
#define oh_transbound(A1,A2)            oh3_transbound(A1,A2)
#else
#ifdef OH_LIB_LEVEL_4P
#define \
oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22) \
oh4p_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22)
#define oh_max_local_particles(A1,A2,A3,A4) \
        oh4p_max_local_particles(A1,A2,A3,A4)
#define oh_per_grid_histogram(A1) oh4p_per_grid_histogram(A1)
#define oh_transbound(A1,A2)      oh4p_transbound(A1,A2)
#define oh_map_particle_to_neighbor(A1,A2,A3) \
        oh4p_map_particle_to_neighbor(A1,A2,A3)
#define oh_map_particle_to_subdomain(A1,A2,A3) \
        oh4p_map_particle_to_subdomain(A1,A2,A3)
#define oh_inject_particle(A1,A2) oh4p_inject_particle(A1,A2)
#define oh_remove_mapped_particle(A1,A2,A3) \
        oh4p_remove_mapped_particle(A1,A2,A3)
#define oh_remap_particle_to_neighbor(A1,A2,A3) \
        oh4p_remap_particle_to_neighbor(A1,A2,A3)
#define oh_remap_particle_to_subdomain(A1,A2,A3) \
        oh4p_remap_particle_to_subdomain(A1,A2,A3)
#else
#define \
oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26) \
oh4s_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26)
#define oh_particle_buffer(A1,A2) oh4s_per_grid_histogram(A1,A2)
#define oh_per_grid_histogram(A1,A2) oh4s_per_grid_histogram(A1,A2)
#define oh_transbound(A1,A2)      oh4s_transbound(A1,A2)
#define oh_exchange_border_data(A1,A2,A3,A4) \
        oh4s_exchange_border_data(A1,A2,A3,A4)
#define oh_map_particle_to_neighbor(A1,A2,A3) \
        oh4s_map_particle_to_neighbor(A1,A2,A3)
#define oh_map_particle_to_subdomain(A1,A2,A3) \
        oh4s_map_particle_to_subdomain(A1,A2,A3)
#define oh_inject_particle(A1,A2) oh4s_inject_particle(A1,A2)
#define oh_remove_mapped_particle(A1,A2,A3) \
        oh4s_remove_mapped_particle(A1,A2,A3)
#define oh_remap_particle_to_neighbor(A1,A2,A3) \
        oh4s_remap_particle_to_neighbor(A1,A2,A3)
#define oh_remap_particle_to_subdomain(A1,A2,A3) \
        oh4s_remap_particle_to_subdomain(A1,A2,A3)
#endif
#endif
#endif
#endif
