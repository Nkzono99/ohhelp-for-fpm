/* File: ohhelp1.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdarg.h>
#include <mpi.h>

#include "oh_config.h"
#include "oh_stats.h"

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifdef  OH_LIB_LEVEL_4PS
#define OH_POS_AWARE
#endif

/* constants for D-dimensional simulation */
#define OH_DIM_X        0
#define OH_DIM_Y        1
#define OH_DIM_Z        2
#if OH_DIMENSION==1
#define OH_NEIGHBORS    3
#elif OH_DIMENSION==2
#define OH_NEIGHBORS    (3*3)
#else
#define OH_NEIGHBORS    (3*3*3)
#endif

MPI_Comm fam_comm;

#define MCW fam_comm      /* shorthand of communicator, MPI_COMM_WORLD or CTCA_subcomm */

typedef long long int dint;     /* shorthand of 64-bit integer */

#ifndef EXTERN
#define EXTERN extern
#endif

/* Basic process configuration variables */
EXTERN int nOfNodes;
EXTERN int myRank;
EXTERN int RegionId[2], *SubdomainId;
#define MODE_NORM_PRI (0)
#define MODE_NORM_SEC (1)
#define MODE_REB_SEC  (-1)
#define MODE_ANY_PRI  (2)
#define MODE_ANY_SEC  (3)
#define Mode_PS(M)       (M&1)
#define Mode_Acc(M)      (M&2)
#define Mode_Set_Pri(M)  (M&2)
#define Mode_Set_Sec(M)  (M|1)
#define Mode_Set_Norm(M) (M&1)
#define Mode_Set_Any(M)  (M|2)
#define Mode_Is_Norm(M)  (M<2)
#define Mode_Is_Any(M)   (M>=2)
EXTERN int currMode, accMode;

/* Number of particles and related variables */
EXTERN int  nOfSpecies;
EXTERN int  maxFraction;
EXTERN int  *NOfPLocal;                 /* [2][nOfSpecies][nOfNodes] */
EXTERN int  *NOfPrimaries;              /* [2][nOfSpecies][nOfNodes] */
EXTERN dint *TotalPGlobal;              /* [nOfNodes+1] */
EXTERN dint nOfParticles;
EXTERN int  nOfLocalPMax;
EXTERN dint *NOfPToStay;                /* [nOfNodes] */
EXTERN int  *TotalP;                    /* [2][nOfSpecies] */
EXTERN int  *TotalPNext;                /* [2][nOfSpecies] */
EXTERN int  primaryParts, totalParts;
EXTERN int  *NOfRecv, *RecvCounts;      /* [2][nOfSpecies][nOfNodes] */
EXTERN int  *NOfSend, *SendCounts;      /* [2][nOfSpecies][nOfNodes] */
EXTERN int  *InjectedParticles;         /* [2][2][nOfSpecies] */
EXTERN int  *TempArray;                 /* [nOfNodes] */
EXTERN MPI_Datatype T_Histogram;

/* Computation node descriptors */
struct S_node {
  struct {int prime, sec;} stay;
  struct {dint prime, sec;} get;
  struct {int prime, sec, black, rank;} comm;
  struct S_node *parent, *sibling, *child;
  int id, parentid;
};
EXTERN struct S_node *Nodes, *NodesNext, **NodeQueue;

/* Heap structure for load rebalancing */
struct S_heap {
  int n, *node, *index;
};
EXTERN struct S_heap LessHeap, GreaterHeap;

/* Structured variables for particle transfer */
struct S_commlist {
  int sid, rid, region, count, tag;     /* tag = spec + nOfSpecies*sec */
};
EXTERN struct S_commlist *CommList, *SecRList;
EXTERN int RLIndex[OH_NEIGHBORS+1];
EXTERN int SLHeadTail[2], SecSLHeadTail[2], SecRLSize;
EXTERN MPI_Datatype T_Commlist;
struct S_commsched_context {
  int neighbor, sender, spec, comidx, dones, donen;
};

/* Structured variables for MPI communicator */
EXTERN MPI_Group GroupWorld;
EXTERN struct {
  int n;
  MPI_Comm *body;       /* [nOfNodes] */
} Comms;
struct S_mycommc {
  MPI_Comm prime, sec;
  int rank, root, black;
};
EXTERN struct S_mycommc *MyComm, *MyCommC;
EXTERN struct S_mycommf {
  int prime, sec;
  int rank, root, black;
} *MyCommF;

/* Neighboring information */
EXTERN int Neighbors[3][OH_NEIGHBORS], SrcNeighbors[OH_NEIGHBORS];
  /* <BSW,BS,BSE,BW,B,BE,BNW,BN,BNE,    : 00..04..08
       SW, S, SE, W,O  E, NW, N, NE,    : 09..13..17
      TSW,TS,TSE,TW,T,TE,TNW,TN,TNE>    : 18..22..26 */
EXTERN int *DstNeighbors;

/* Structures and variables for statistics and verbose messaging */
#define STATS_PART_MOVE_PRI_MIN 0
#define STATS_PART_MOVE_PRI_MAX 1
#define STATS_PART_MOVE_PRI_AVE 2
#define STATS_PART_GET_PRI_MIN  3
#define STATS_PART_GET_PRI_MAX  4
#define STATS_PART_PUT_PRI_MIN  5
#define STATS_PART_PUT_PRI_MAX  6
#define STATS_PART_PG_PRI_AVE   7
#define STATS_PART_MOVE_SEC_MIN 8
#define STATS_PART_MOVE_SEC_MAX 9
#define STATS_PART_MOVE_SEC_AVE 10
#define STATS_PART_GET_SEC_MIN  11
#define STATS_PART_GET_SEC_MAX  12
#define STATS_PART_PUT_SEC_MIN  13
#define STATS_PART_PUT_SEC_MAX  14
#define STATS_PART_PG_SEC_AVE   15
#define STATS_PART_PRIMARY      16
#define STATS_PART_SECONDARY    17
#define STATS_PARTS             (STATS_PART_SECONDARY+1)

#ifdef OH_DEFINE_STATS
static char *StatsPartStrings[STATS_PARTS] = {
  "p2p transfer[pri,min]",
  "p2p transfer[pri,max]",
  "p2p transfer[pri,ave]",
  "get[pri,min]",
  "get[pri,max]",
  "put[pri,min]",
  "put[pri,max]",
  "put&get[pri,ave]",
  "p2p transfer[sec,min]",
  "p2p transfer[sec,max]",
  "p2p transfer[sec,ave]",
  "get[sec,min]",
  "get[sec,max]",
  "put[sec,min]",
  "put[sec,max]",
  "put&get[sec,ave]",
  "transition to pri",
  "transition to sec",
};
#endif

struct S_statscurr {
  struct {
    double value, val[2*STATS_TIMINGS+2];
    int key, ev[2*STATS_TIMINGS+2];
  } time;
  dint part[STATS_PARTS];
};
struct S_statstime {
  double min, max, total;
  int ev;
};
struct S_statspart {
  dint min, max, total;
};
struct S_statstotal {
  struct S_statstime time[2*STATS_TIMINGS];
  struct S_statspart part[STATS_PARTS];
};
EXTERN struct S_stats {
  struct S_statscurr curr;
  struct S_statstotal subtotal, total;
} Stats;

EXTERN MPI_Datatype T_StatsTime;
EXTERN MPI_Op Op_StatsTime, Op_StatsPart;
EXTERN int statsMode, reportIteration, verboseMode;

/* Prototypes for the functions called from simulator code */
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

void oh1_fam_comm(MPI_Comm *fortran_comm);

void oh1_init(int **sdid, int nspec, int maxfrac, int **nphgram,
              int **totalp, int **rcounts, int **scounts, void *mycomm,
              int **nbor, int *pcoord, int stats, int repiter, int verbose);
int  oh1_transbound(int currmode, int stats);

void oh1_neighbors_(int *nbor);
void oh1_families_(int *famindex, int *members);
int  oh1_accom_mode_();
void oh1_broadcast_(void *pbuf, void *sbuf, int *pcount, int *scount,
                    int *ptype, int *stype);
void oh1_all_reduce_(void *pbuf, void *sbuf, int *pcount, int *scount,
                     int *ptype, int *stype, int *pop, int *sop);
void oh1_reduce_(void *pbuf, void *sbuf, int *pcount, int *scount,
                 int *ptype, int *stype, int *pop, int *sop);
void oh1_init_stats_(int *key, int *ps);
void oh1_stats_time_(int *key, int *ps);
void oh1_show_stats_(int *step, int *currmode);
void oh1_print_stats_(int *nstep);
void oh1_verbose_(char *message);
void oh1_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
               int *totalp, int *rcounts, int *scounts,
               struct S_mycommf *mycomm, int *nbor, int *pcoord, int *stats,
               int *repiter, int *verbose);
int  oh1_transbound_(int *currmode, int *stats);

/* Prototypes for the functions called from higher-level library code */
void  init1(int **sdid, int nspec, int maxfrac, int **nphgram,
            int **totalp, int **rcounts, int **scounts,
            struct S_mycommc *mycommc, struct S_mycommf *mycommf, int **nbor,
            int *pcoord, int stats, int repiter, int verbose);
void* mem_alloc(int esize, int count, char* varname);
void  mem_alloc_error(char* varname, size_t size);
void  errstop(char* format, ...);
void  local_errstop(char* format, ...);
void  set_total_particles();
int   transbound1(int currmode, int stats, int level);
int   try_primary1(int currmode, int level, int stats);
int   try_stable1(int currmode, int level, int stats);
void  rebalance1(int currmode, int level, int stats);
void  build_new_comm(int currmode, int level, int nbridx, int stats);
void  vprint(char* format, ...);
void  dprint(char* format, ...);

/* Macro for verbose messaging. */
#define Verbose(L,VP) {\
  if (verboseMode>=L) {\
    MPI_Barrier(MCW);\
    if (myRank==0 || verboseMode>=3) VP;\
  }\
}
