/* File: ohhelp3.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
EXTERN int excludeLevel2;

#define OH_LOWER 0
#define OH_UPPER 1

EXTERN int (*SubDomains)[OH_DIMENSION][2];      /* [N][D][l,u] */
EXTERN double (*SubDomainsFloat)[OH_DIMENSION][2];
struct S_grid {
  int n, coord[2], size;
  double fcoord[2], fsize, gsize, rgsize;
  struct {
    int size, n, thresh;
    double rfsize, rfsizeplus, fthresh;
  } light;
};
EXTERN struct S_grid Grid[3];

struct S_subdomdesc {
  struct {
    int c[2], h, n;
    double fc[2];
  } coord[OH_DIMENSION];
  int id;
};
EXTERN struct S_subdomdesc *SubDomainDesc;

static struct S_message {
  char xyz[4];
  char loup[2][6];
} Message = {
  "xyz",
  {"lower", "upper"}
};

EXTERN int nOfBoundaries;
EXTERN int (*Boundaries)[OH_DIMENSION][2];      /* [N][D][l,u] */
EXTERN int Adjacent[OH_DIMENSION][2];           /* [D][l,u] */

EXTERN int nOfFields;

#define OH_FTYPE_ES 0
#define OH_FTYPE_LO 1
#define OH_FTYPE_UP 2
#define OH_FTYPE_BL 3
#define OH_FTYPE_BU 4
#define OH_FTYPE_RL 5
#define OH_FTYPE_RU 6
#define OH_FTYPE_N  7
EXTERN int (*FieldTypes)[OH_FTYPE_N];           /* [F][es,lo,up,bl,bu,rl,ru] */

struct S_brdesc {
  int base, size[2];
};
struct S_flddesc {
  int esize, ext[2], size[OH_DIMENSION];
  struct S_brdesc bc, red;
};
EXTERN struct S_flddesc *FieldDesc;             /* [F] */

EXTERN int nOfExc;
EXTERN int *BoundaryCommFields;                 /* [C] */

#define OH_CTYPE_FROM 0
#define OH_CTYPE_TO   1
#define OH_CTYPE_SIZE 2
#define OH_CTYPE_N    3
EXTERN int (*BoundaryCommTypes)[2][OH_CTYPE_N]; /* [C][B][d,u][from,to,size] */

struct S_bcomm {
  int buf, count, deriv;
  MPI_Datatype type;
};
struct S_borderexc {
  struct S_bcomm send, recv;
};
EXTERN struct S_borderexc (*BorderExc)[2][OH_DIMENSION][2];
                                                        /* [C][ps][D][l,u] */

/* Prototypes for the functions called from simulator code */
void oh3_grid_size(double size[OH_DIMENSION]);
void oh3_bcast_field(void *pfld, void *sfld, int ftype);
void oh3_reduce_field(void *pfld, void *sfld, int ftype);
void oh3_allreduce_field(void *pfld, void *sfld, int ftype);
void oh3_exchange_borders(void *pfld, void *sfld, int ctype, int bcast);

void oh3_grid_size_(double size[OH_DIMENSION]);
void oh3_bcast_field_(void *pfld, void *sfld, int *ftype);
void oh3_reduce_field_(void *pfld, void *sfld, int *ftype);
void oh3_allreduce_field_(void *pfld, void *sfld, int *ftype);
void oh3_exchange_borders_(void *pfld, void *sfld, int *ctype, int *bcast);

#if OH_DIMENSION==1
int  oh3_map_particle_to_neighbor(double *x, int ps);
int  oh3_map_particle_to_subdomain(double x);
int  oh3_map_region_to_adjacent_node_(double *x, int *ps);
int  oh3_map_particle_to_neighbor_(double *x, int *ps);
int  oh3_map_region_to_node_(double *x);
int  oh3_map_particle_to_subdomain_(double *x);
#elif OH_DIMENSION==2
int  oh3_map_particle_to_neighbor(double *x, double *y, int ps);
int  oh3_map_particle_to_subdomain(double x, double y);
int  oh3_map_region_to_adjacent_node_(double *x, double *y, int *ps);
int  oh3_map_particle_to_neighbor_(double *x, double *y, int *ps);
int  oh3_map_region_to_node_(double *x, double *y);
int  oh3_map_particle_to_subdomain_(double *x, double *y);
#else
int  oh3_map_particle_to_neighbor(double *x, double *y, double *z, int ps);
int  oh3_map_particle_to_subdomain(double x, double y, double z);
int  oh3_map_region_to_adjacent_node_(double *x, double *y, double *z,
                                      int *ps);
int  oh3_map_particle_to_neighbor_(double *x, double *y, double *z,
                                   int *ps);
int  oh3_map_region_to_node_(double *x, double *y, double *z);
int  oh3_map_particle_to_subdomain_(double *x, double *y, double *z);
#endif

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

void oh3_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
               int *totalp, struct S_particle *pbuf, int *pbase,
               int *maxlocalp, struct S_mycommf *mycomm, int *nbor,
               int *pcoord, int *sdoms, int *scoord, int *nbound, int *bcond,
               int *bounds, int *ftypes, int *cfields, int *ctypes,
               int *fsizes, int *stats, int *repiter, int *verbose);
void oh13_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
                int *totalp, int *rcounts, int *scounts,
                struct S_mycommf *mycomm, int *nbor, int *pcoord,
                int *sdoms, int *scoord, int *nbound, int *bcond, int *bounds,
                int *ftypes, int *cfields, int *ctypes, int *fsizes,
                int *stats, int *repiter, int *verbose);
int  oh3_transbound_(int *currmode, int *stats);

/* Prototype for the function called from higher-level library code */
void init3(int **sdid, int nspec, int maxfrac, int **nphgram, int **totalp,
           int **rcounts, int **scounts, struct S_particle **pbuf, int **pbase,
           int maxlocalp, struct S_mycommc *mycommc, struct S_mycommf *mycommf,
           int **nbor, int *pcoord, int **sdoms, int *scoord, int nbound,
           int *bcond, int **bounds, int *ftypes, int *cfields, int cfid,
           int *ctypes, int **fsizes, int stats, int repiter, int verbose,
           int skip2);
void set_field_descriptors(int (*ft)[OH_FTYPE_N], int sd[OH_DIMENSION][2],
                           int ps);
void clear_border_exchange();
int  map_irregular_subdomain(double x, double y, double z);
