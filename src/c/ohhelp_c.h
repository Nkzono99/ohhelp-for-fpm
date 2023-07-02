/* File: ohhelp_c.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#include <mpi.h>
struct S_mycommc {
  MPI_Comm prime, sec;
  int rank, root, black;
};
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

void oh1_neighbors(int **nbor);
void oh1_families(int **famindex, int **members);
int  oh1_accom_mode();
void oh1_broadcast(void *pbuf, void *sbuf, int pcount, int scount,
                   MPI_Datatype ptype, MPI_Datatype stype);
void oh1_all_reduce(void *pbuf, void *sbuf, int pcount, int scount,
                    MPI_Datatype ptype, MPI_Datatype stype,
                    MPI_Op pop, MPI_Op sop);
void oh1_reduce(void *pbuf, void *sbuf, int pcount, int scount,
                MPI_Datatype ptype, MPI_Datatype stype,
                MPI_Op pop, MPI_Op sop);
void oh1_init_stats(int key, int ps);
void oh1_stats_time(int key, int ps);
void oh1_show_stats(int step, int currmode);
void oh1_print_stats(int nstep);
void oh1_verbose(char *message);

#if OH_LIB_LEVEL==1
#define oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13) \
        oh1_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13)
#define oh_transbound(A1,A2) oh1_transbound(A1,A2)
void oh1_init(int **sdid, int nspec, int maxfrac, int **nphgram,
              int **totalp, int **rcounts, int **scounts, void *mycomm,
              int **nbor, int *pcoord, int stats, int repiter, int verbose);
int  oh1_transbound(int currmode, int stats);

#else
#define oh_set_total_particles() oh2_set_total_particles()
#include "oh_part.h"
void oh2_set_total_particles();

#if OH_LIB_LEVEL!=4
#define oh_max_local_particles(A1,A2,A3) oh2_max_local_particles(A1,A2,A3)
#define oh_inject_particle(A1)           oh2_inject_particle(A1)
#define oh_remap_injected_particle(A1)   oh2_remap_injected_particle(A1)
#define oh_remove_injected_particle(A1)  oh2_remove_injected_particle(A1)
int  oh2_max_local_particles(long long int npmax, int maxfrac, int minmargin);
void oh2_remap_injected_particle(struct S_particle *part);
void oh2_remove_injected_particle(struct S_particle *part);
#endif
void oh2_inject_particle(struct S_particle *part);

#if OH_LIB_LEVEL==2
#define oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14) \
        oh2_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)
#define oh_transbound(A1,A2) oh2_transbound(A1,A2)
void oh2_init(int **sdid, int nspec, int maxfrac, int **nphgram,
              int **totalp, struct S_particle **pbuf, int **pbase,
              int maxlocalp, void *mycomm, int **nbor,
              int *pcoord, int stats, int repiter, int verbose);
int  oh2_transbound(int currmode, int stats);

#else
#define oh_grid_size(A1)                 oh3_grid_size(A1)
#define oh_bcast_field(A1,A2,A3)         oh3_bcast_field(A1,A2,A3)
#define oh_reduce_field(A1,A2,A3)        oh3_reduce_field(A1,A2,A3)
#define oh_allreduce_field(A1,A2,A3)     oh3_allreduce_field(A1,A2,A3)
#define oh_exchange_borders(A1,A2,A3,A4) oh3_exchange_borders(A1,A2,A3,A4)
void oh3_grid_size(double size[OH_DIMENSION]);
void oh3_bcast_field(void *pfld, void *sfld, int ftype);
void oh3_reduce_field(void *pfld, void *sfld, int ftype);
void oh3_allreduce_field(void *pfld, void *sfld, int ftype);
void oh3_exchange_borders(void *pfld, void *sfld, int ctype, int bcast);

#if OH_LIB_LEVEL==3
#if OH_DIMENSION==1
#define oh_map_particle_to_neighbor(A1,A2) \
        oh3_map_particle_to_neighbor(A1,A2)
#define oh_map_particle_to_subdomain(A1) \
        oh3_map_particle_to_subdomain(A1)
int  oh3_map_particle_to_neighbor(double *x, int ps);
int  oh3_map_particle_to_subdomain(double x);
#elif OH_DIMENSION==2
#define oh_map_particle_to_neighbor(A1,A2,A3) \
        oh3_map_particle_to_neighbor(A1,A2,A3)
#define oh_map_particle_to_subdomain(A1,A2) \
        oh3_map_particle_to_subdomain(A1,A2)
int  oh3_map_particle_to_neighbor(double *x, double *y, int ps);
int  oh3_map_particle_to_subdomain(double x, double y);
#else
#define oh_map_particle_to_neighbor(A1,A2,A3,A4) \
        oh3_map_particle_to_neighbor(A1,A2,A3,A4)
#define oh_map_particle_to_subdomain(A1,A2,A3) \
        oh3_map_particle_to_subdomain(A1,A2,A3)
int  oh3_map_particle_to_neighbor(double *x, double *y, double *z, int ps);
int  oh3_map_particle_to_subdomain(double x, double y, double z);
#endif

#define \
oh_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23) \
oh3_init(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23)
#define oh_transbound(A1,A2)            oh3_transbound(A1,A2)
void oh3_init(int **sdid, int nspec, int maxfrac, int **nphgram, int **totalp,
              struct S_particle **pbuf, int **pbase, int maxlocalp,
              void *mycomm, int **nbor, int *pcoord,
              int **sdoms, int *scoord, int nbound, int *bcond, int **bounds,
              int *ftypes, int *cfields, int *ctypes, int **fsizes,
              int stats, int repiter, int verbose);
void oh13_init(int **sdid, int nspec, int maxfrac, int **nphgram,
               int **totalp, int **rcounts, int **scounts,
               void *mycomm, int **nbor, int *pcoord,
               int **sdoms, int *scoord, int nbound, int *bcond, int **bounds,
               int *ftypes, int *cfields, int *ctypes, int **fsizes,
               int stats, int repiter, int verbose);
int  oh3_transbound(int currmode, int stats);

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
void oh4p_init(int **sdid, const int nspec, const int maxfrac, int **totalp,
               struct S_particle **pbuf, int **pbase, const int maxlocalp,
               void *mycomm, int **nbor, int *pcoord, int **sdoms, int *scoord,
               const int nbound, int *bcond, int **bounds, int *ftypes,
               int *cfields, int *ctypes, int **fsizes,
               const int stats, const int repiter, const int verbose);
int  oh4p_max_local_particles(const long long int npmax, const int maxfrac,
                              const int minmargin, const int hsthresh);
void oh4p_per_grid_histogram(int **pghgram);
int  oh4p_transbound(int currmode, int stats);
int  oh4p_map_particle_to_neighbor(struct S_particle *part, const int ps,
                                   const int s);
int  oh4p_map_particle_to_subdomain(struct S_particle *part, const int ps,
                                    const int s);
int  oh4p_inject_particle(const struct S_particle *part, const int ps);
void oh4p_remove_mapped_particle(struct S_particle *part, const int ps,
                                 const int s);
int  oh4p_remap_particle_to_neighbor(struct S_particle *part, const int ps,
                                     const int s);
int  oh4p_remap_particle_to_subdomain(struct S_particle *part, const int ps,
                                      const int s);
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
void oh4s_init(int **sdid, const int nspec, const int maxfrac,
               const long long int npmax, const int minmargin,
               const int maxdensity, int **totalp, int **pbase,
               int *maxlocalp, int *cbufsize, void *mycomm, int **nbor,
               int *pcoord, int **sdoms, int *scoord, const int nbound,
               int *bcond, int **bounds, int *ftypes, int *cfields,
               int *ctypes, int **fsizes, int **zbound,
               const int stats, const int repiter, const int verbose);
void oh4s_particle_buffer(const int maxlocalp, struct S_particle **pbuf);
void oh4s_per_grid_histogram(int **pghgram, int **pgindex);
int  oh4s_transbound(int currmode, int stats);
void oh4s_exchange_border_data(void *buf, void *sbuf, void *rbuf,
                               MPI_Datatype type);
int  oh4s_map_particle_to_neighbor(struct S_particle *part, const int ps,
                                   const int s);
int  oh4s_map_particle_to_subdomain(struct S_particle *part, const int ps,
                                    const int s);
int  oh4s_inject_particle(const struct S_particle *part, const int ps);
void oh4s_remove_mapped_particle(struct S_particle *part, const int ps,
                                 const int s);
int  oh4s_remap_particle_to_neighbor(struct S_particle *part, const int ps,
                                     const int s);
int  oh4s_remap_particle_to_subdomain(struct S_particle *part, const int ps,
                                      const int s);
#endif
#endif
#endif
#endif
