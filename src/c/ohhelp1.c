/* File: ohhelp1.c
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#define EXTERN
#define OH_DEFINE_STATS
#include "ohhelp1.h"

/* Prototypes for private functions. */
static void  count_stay();
static dint  assign_particles(dint npr, dint npt, struct S_node *ch, int incgp,
                              int *nget);
static int   compare_int(const void* x, const void* y);
static void  schedule_particle_exchange(int reb);
static int   count_real_stay(int *np);
static void  sched_comm(int toget, int rid, int tag, int reb,
                        struct S_commsched_context *context);
static void  make_comm_count(int currmode, int level, int reb, int oldparent,
                             int stats);
static void  make_recv_count(struct S_commlist* rlist, int rlsize);
static void  make_send_count(struct S_commlist* slist, int slsize);
static void  count_next_particles(struct S_commlist* rlist, int rlsize);
static void  push_heap(int id, struct S_heap* heap, int greater);
static int   pop_heap(struct S_heap* heap, int greater);
static void  remove_heap(struct S_heap* heap, int greater, int rem);
static void  clear_stats(struct S_statstotal *stotal);
static void  stats_primary_comm(int currmode);
static void  stats_secondary_comm(int currmode, int reb);
static void  stats_comm(int* nrecv, int* nsend, dint* scp, int ns);
static void  update_stats(struct S_statstotal *stotal, int step, int currmode);
static void  stats_reduce_part(void* inarg, void* ioarg, int* len,
                               MPI_Datatype* type);
static void  print_stats(struct S_statstotal *stotal, int cstep, int n);
static void  stats_reduce_time(void* inarg, void* ioarg, int* len,
                               MPI_Datatype* type);

void
oh1_fam_comm_(MPI_Comm *fortran_comm) {
  fam_comm = *fortran_comm;
}

void
oh1_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
          int *totalp, int *rcounts, int *scounts, struct S_mycommf *mycomm,
          int *nbor, int *pcoord, int *stats, int *repiter, int *verbose) {
  init1(&sdid, *nspec, *maxfrac, &nphgram, &totalp, &rcounts, &scounts,
        NULL, mycomm, &nbor, pcoord, *stats, *repiter, *verbose);
}
void
oh1_init(int **sdid, int nspec, int maxfrac, int **nphgram,
         int **totalp, int **rcounts, int **scounts, void *mycomm,
         int **nbor, int *pcoord, int stats, int repiter, int verbose) {
  init1(sdid, nspec, maxfrac, nphgram, totalp, rcounts, scounts,
        (struct S_mycommc*)mycomm, NULL, nbor, pcoord, stats, repiter,
        verbose);
}
static int (*NeighborsShadow)[OH_NEIGHBORS] = NULL;
static int *NeighborsTemp = NULL;
void
init1(int **sdid, int nspec, int maxfrac, int **nphgram,
      int **totalp, int **rcounts, int **scounts, struct S_mycommc *mycommc,
      struct S_mycommf *mycommf, int **nbor, int *pcoord,
      int stats, int repiter, int verbose) {

  int nn, ns, me, i, s, clsize;
  int *nb = *nbor;
  int bl[2]={1,1};
  MPI_Datatype tmptype[2]={MPI_DATATYPE_NULL, MPI_UB};
  MPI_Aint disp[2]={0, sizeof(int)};

  MPI_Comm_size(MCW, &nn);  nOfNodes = nn;
  MPI_Comm_rank(MCW, &me);  myRank = me;
  currMode = MODE_NORM_PRI;  accMode = 0;

  verboseMode = verbose;
  Verbose(1, vprint("oh_init"));

  if (!*sdid)  *sdid = (int*)mem_alloc(sizeof(int), 2, "SubdomainID");
  SubdomainId = *sdid;
  (*sdid)[0] = RegionId[0] = me;  (*sdid)[1] = RegionId[1] = -1;
  ns = nOfSpecies = nspec;
  maxFraction = maxfrac;

  if (!*nphgram)
    *nphgram = (int*)mem_alloc(sizeof(int), 2*ns*nn, "NOfPLocal");
  NOfPLocal = *nphgram;
  if (!*totalp)  *totalp = (int*)mem_alloc(sizeof(int), 2*ns, "TotalP");
  TotalPNext = *totalp;
  TotalP = NULL;
  for(i=0; i<2*ns*nn; i++)  NOfPLocal[i] = 0;

  if (rcounts) {
    if (!*rcounts)
      *rcounts = (int*)mem_alloc(sizeof(int), 2*ns*nn, "RecvCounts");
    RecvCounts = *rcounts;
  }
  if (scounts) {
    if (!*scounts)
      *scounts = (int*)mem_alloc(sizeof(int), 2*ns*nn, "SendCounts");
    SendCounts = *scounts;
  }
  NOfRecv = (int*)mem_alloc(sizeof(int), 2*ns*nn, "NOfRecv");
  NOfSend = (int*)mem_alloc(sizeof(int), 2*ns*nn, "NOfSend");

  NOfPrimaries = (int*) mem_alloc(sizeof(int),  2*ns*nn, "NOfPrimaries");
  TotalPGlobal = (dint*)mem_alloc(sizeof(dint), nn+1, "TotalPGlobal");
  NOfPToStay   = (dint*)mem_alloc(sizeof(dint), nn, "NOfPToStay");
  InjectedParticles = (int*)mem_alloc(sizeof(int), 4*ns, "InjectedParticles");
  for (s=0; s<ns*2; s++)  InjectedParticles[s] = 0;
#ifndef OH_POS_AWARE
  TempArray    = (int*) mem_alloc(sizeof(int),  nn, "TempArray");
#endif

  MPI_Type_vector(2*ns, 1, nn, MPI_INT, tmptype);
  MPI_Type_struct(2, bl, disp, tmptype, &T_Histogram);
  MPI_Type_commit(&T_Histogram);

  Nodes = (struct S_node*)mem_alloc(sizeof(struct S_node), nn, "Nodes");
  NodesNext = (struct S_node*)mem_alloc(sizeof(struct S_node), nn,
                                        "NodesNext");
  for (i=0; i<nn; i++) Nodes[i].id = i;
  NodeQueue =(struct S_node**)mem_alloc(sizeof(struct S_node*), nn,
                                        "NodeQueue");

  LessHeap.node     = (int*)mem_alloc(sizeof(int), nn*2, "LessHeap") - 1;
  LessHeap.index    = LessHeap.node + nn + 1;
  GreaterHeap.node  = (int*)mem_alloc(sizeof(int), nn*2, "GreaterHeap") - 1;
  GreaterHeap.index = GreaterHeap.node + nn + 1;

  clsize = 2*OH_NEIGHBORS*(nn*ns+1)+nn*(ns+3);
  if (clsize<(14+4*ns)*nn)  clsize = (14+4*ns)*nn;
  CommList = (struct S_commlist*)mem_alloc(sizeof(struct S_commlist), clsize,
                                           "CommList");
  MPI_Type_contiguous(sizeof(struct S_commlist), MPI_BYTE, &T_Commlist);
  MPI_Type_commit(&T_Commlist);

  Comms.body = (MPI_Comm*)mem_alloc(sizeof(MPI_Comm), nn, "Comms");
  Comms.n = 0;
  MPI_Comm_group(MCW, &GroupWorld);
  MyComm = (struct S_mycommc*)mem_alloc(sizeof(struct S_mycommc), 1, "MyComm");
  MyComm->prime = MyComm->sec  = MPI_COMM_NULL;
  MyComm->rank  = MyComm->root = MyComm->black = 0;
  if ((MyCommC=mycommc))  *mycommc = *MyComm;
  if ((MyCommF=mycommf)) {
    MyCommF->prime = MyCommF->sec = MPI_Comm_c2f(MPI_COMM_NULL);
    MyCommF->rank = MyCommF->root = MyCommF->black = 0;
  }
  if (!nb) {
    if (NeighborsShadow) {
      nb = *nbor = (int*)mem_alloc(sizeof(int), OH_NEIGHBORS, "Neighbors");
    } else {
      nb = *nbor = (int*)mem_alloc(sizeof(int), 3*OH_NEIGHBORS, "Neighbors");
    }
    nb[0] = -1;
  }
  if (nb[0]==-1) {
    int p=pcoord[0];
    int q=(OH_DIMENSION>1)?pcoord[1]:1, r=(OH_DIMENSION>2)?pcoord[2]:1;
    int j, k, l;
    int yplus=(OH_DIMENSION>1)?2:0, zplus=(OH_DIMENSION>2)?2:0;
    int xoff, yoff, zoff;
    if (nn!=p*q*r || p<0 || q<0 || r<0)
      errstop("<# of x-nodes>(%d) * <# of y-nodes>(%d) * <# of z-nodes>(%d) "
              "should be equal to <# of nodes>(%d)", p, q, r, nn);
    i = me % p;  j = (me/p) % q;  k = me / (p*q);
    for (l=0,zoff=-1; zoff<zplus; zoff++) {
      for (yoff=-1; yoff<yplus; yoff++) {
        for (xoff=-1; xoff<2; xoff++,l++) {
          nb[l] = (i+xoff+p)%p + (((j+yoff+q)%q) + ((k+zoff+r)%r)*q)*p;
        }
      }
    }
  } else {
    for (i=0; i<OH_NEIGHBORS; i++) {
      int n=nb[i], m=nb[(OH_NEIGHBORS-1)-i], k;
      MPI_Status st;
      if (m>=0) {
        if (n>=0)
          MPI_Sendrecv(&i, 1, MPI_INT, n, 0, &k, 1, MPI_INT, m, 0, MCW, &st);
        else
          MPI_Recv(&k, 1, MPI_INT, m, 0, MCW, &st);
        if (k!=i)
          local_errstop("rank-%d's %d-th neighbor rank-%d says "
                        "rank-%d is not %d-th neighbor but %d-th",
                        me, (OH_NEIGHBORS-1)-i, m, me, i, k);
      } else if (n>=0) {
        MPI_Send(&i, 1, MPI_INT, n, 0, MCW);
      }
    }
  }
  NeighborsTemp = nb;
  if (NeighborsShadow && nb!=(int*)NeighborsShadow)
    for (i=0; i<OH_NEIGHBORS; i++)  NeighborsShadow[0][i] = nb[i];

  DstNeighbors = Neighbors[0];
  for (i=0; i<nn; i++) TempArray[i] = 0;
  for (i=0; i<OH_NEIGHBORS; i++) {
    int dst=nb[i],  src=nb[(OH_NEIGHBORS-1)-i];
    if (dst<0)  DstNeighbors[i] = -(nn+1);
    else {
      DstNeighbors[i] = (TempArray[dst]&1) ? -(dst+1) : dst;
      TempArray[dst] |= 1;
    }
    if (src<0)
      SrcNeighbors[i] = -(nn+1);
    else {
      SrcNeighbors[i] = (TempArray[src]&2) ? -(src+1) : src;
      TempArray[src] |= 2;
    }
  }
  statsMode = stats;
  reportIteration = repiter;
}
void*
mem_alloc(int esize, int count, char* varname) {

  size_t size = (size_t)esize*(size_t)count;
  void* ptr = malloc(size);
  if (!ptr) mem_alloc_error(varname, size);
  return(ptr);
}
void
mem_alloc_error(char* varname, size_t size) {
  errstop("out of virtual memory for %s(%lld)", varname, size);
}
void
errstop(char* format, ...) {
  va_list v;
  va_start(v, format);

  if (myRank==0) {
    vfprintf(stderr, format, v);
    fprintf(stderr, "\n");
  }
  va_end(v);
  MPI_Finalize();  exit(1);
}
void
local_errstop(char* format, ...) {
  va_list v;
  va_start(v, format);

  vfprintf(stderr, format, v);
  fprintf(stderr, "\n");
  va_end(v);
  MPI_Abort(MCW, 1);
}
void
oh1_neighbors_(int *nbor) {
  oh1_neighbors(&nbor);
}
void
oh1_neighbors(int **nbor) {
  int *nb = *nbor;
  int i;

  if (!nb)
    nb = *nbor = (int*)mem_alloc(sizeof(int), 3*OH_NEIGHBORS, "Neighbors");
  if (NeighborsTemp && nb!=NeighborsTemp)
    for (i=0; i<OH_NEIGHBORS; i++)  nb[i] = NeighborsTemp[i];
  NeighborsShadow = (int(*)[OH_NEIGHBORS])nb;
}
static int *FamIndex = NULL;
static int *FamMembers = NULL;
void
oh1_families_(int *famindex, int *members) {
  oh1_families(&famindex, &members);
}
void
oh1_families(int **famindex, int **members) {
  int *fidx = *famindex,  *fmem = *members;
  int nn, i;

  MPI_Comm_size(MCW, &nn);
  if (!fidx)
    fidx = *famindex = (int*)mem_alloc(sizeof(int), nn+1, "FamIndex");
  if (!fmem)
    fmem = *members = (int*)mem_alloc(sizeof(int), nn*2, "FamMembers");
  for (i=0; i<nn; i++)  fidx[i] = fmem[i] = i;
  fidx[nn] = nn;
  FamIndex = fidx;  FamMembers = fmem;
}
void
set_total_particles() {
  int ns=nOfSpecies, nn=nOfNodes, nnns=nn*ns;
  int cm=(Mode_PS(currMode))&&(RegionId[1]>=0);
  int s, i, j, tpp, tps;

  if (!TotalP)  TotalP = (int*)mem_alloc(sizeof(int), 2*ns, "TotalP");
  primaryParts = 0;  totalParts = 0;
  for (s=0,j=0; s<ns; s++) {
    for (i=0,tpp=0,tps=0; i<nn; i++,j++) {
      tpp += NOfPLocal[j];  tps += NOfPLocal[nnns+j];
    }
    if (!cm)  tps = 0;
    TotalP[s] = TotalPNext[s] = tpp;  TotalP[ns+s] = TotalPNext[ns+s] = tps;
    primaryParts += tpp;  totalParts += tps;
  }
  totalParts += primaryParts;
}
int
oh1_transbound_(int *currmode, int *stats) {
  return(transbound1(*currmode, *stats, 1));
}
int
oh1_transbound(int currmode, int stats) {
  return(transbound1(currmode, stats, 1));
}
int
transbound1(int currmode, int stats, int level) {
  int ret=MODE_NORM_SEC, nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns, nnns2=2*nnns;
  int i, j, k, s, p, tp, tpn, *nbor;
  dint nofp;

  Verbose(1,vprint("oh_transbound"));
  if ((stats=statsMode&&stats)) oh1_stats_time(STATS_TRANSBOUND, 0);
  if (!TotalP)  set_total_particles();
  currmode = Mode_PS(currmode);
  if (currmode!=Mode_PS(currMode))
    local_errstop("currmode given to oh_transbound() does not match with "
                  "that the library maintains");

  for (i=0; i<nn; i++) TotalPGlobal[i] = 0;
  for (p=0,j=0,tp=0,tpn=0; p<=currmode; p++) {
    for (s=0; s<ns; s++) {
      for (i=0; i<nn; i++,j++) {
        int np=NOfPLocal[j];
        TotalPGlobal[i] += np;  tp += np;
      }
    }
    for (i=0,nbor=Neighbors[p]; i<OH_NEIGHBORS; i++) {
      int n=nbor[i];
      if (n>=0)
        for (s=0,k=(p==0)?n:n+nnns; s<ns; s++,k+=nn)  tpn+= NOfPLocal[k];
    }
  }
  TotalPGlobal[nn] = (tp==tpn && Mode_Is_Norm(currMode)) ? 0 : 1;

  MPI_Alltoall(NOfPLocal, 1, T_Histogram, NOfPrimaries, 1, T_Histogram, MCW);
#ifndef INTEL_MPI_BUG_FIXED
  for (p=0,k=myRank; p<2; p++)  for (s=0; s<ns; s++,k+=nn)
    NOfPrimaries[k] = NOfPLocal[k];
#endif
  MPI_Allreduce(MPI_IN_PLACE, TotalPGlobal, nn+1, MPI_LONG_LONG_INT, MPI_SUM,
                MCW);
  if (TotalPGlobal[nn])  currmode = Mode_Set_Any(currmode);

  for (i=0,nofp=0; i<nn; i++)  nofp += TotalPGlobal[i];
  nOfParticles = nofp;
  nOfLocalPMax = nofp*(maxFraction+100)/100/nn;
  accMode = Mode_Is_Any(currmode) ? 1 : 0;

  if (level>1) return(currmode);
  if (try_primary1(currmode, 1, stats))  ret = MODE_NORM_PRI;
  else if (!Mode_PS(currmode) || !try_stable1(currmode, 1, stats)) {
    rebalance1(currmode, 1, stats);  ret = MODE_REB_SEC;
  }
  for (i=0; i<nnns2; i++) {
    NOfPLocal[i] = 0;  RecvCounts[i] = NOfRecv[i];  SendCounts[i] = NOfSend[i];
  }
  for (s=0; s<ns*2; s++) TotalP[s] = TotalPNext[s];
  return((currMode=ret));
}
int
try_primary1(int currmode, int level, int stats) {
  int nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns, me=myRank, nlpmax=nOfLocalPMax;
  int i, j, s;

  Verbose(2,vprint("try_primary(%s,%s)",
                   Mode_PS(currmode)?"secondary":"primary",
                   Mode_Acc(currmode)?"anywhere":"normal"));
  for (i=0; i<nn; i++) {
    if (TotalPGlobal[i]>nlpmax) return(FALSE);
  }
  if (stats) stats_primary_comm(currmode);
  Verbose(2,vprint("try_primary=TRUE"));

  SubdomainId[1] = RegionId[1] = -1;
  if (Mode_PS(currmode) && FamIndex) {
    int *fidx = FamIndex,  *fmem = FamMembers;
    for (i=0; i<nn; i++)  fidx[i] = fmem[i] = i;
    fidx[nn] = nn;
  }
  if (level>1) return(TRUE);
  for (i=0,s=0; s<ns; s++) {
    int t = 0;
    for (j=0; j<nn; j++,i++) {
      NOfSend[i] = NOfPLocal[i] + NOfPLocal[i+nnns];
      t += (NOfRecv[i] = NOfPrimaries[i] + NOfPrimaries[i+nnns]);
      NOfSend[i+nnns] = NOfRecv[i+nnns] = 0;
      if (j==me) {
        NOfSend[i] -= NOfPLocal[i];  NOfRecv[i] -= NOfPrimaries[i];
      }
    }
    TotalPNext[s] = t;  TotalPNext[ns+s] = 0;
  }
  return(TRUE);
}
#define Special_Pexc_Sched(LEVEL) (LEVEL<0)
int
try_stable1(int currmode, int level, int stats) {
  int nn=nOfNodes;
  int nlpmax=nOfLocalPMax;
  struct S_node *node, *ch;
  int i;

  if (stats) oh1_stats_time(STATS_TRY_STABLE, 0);
  Verbose(2,vprint("try_stable"));
  count_stay();

  for (i=0; ; i++) {                    /* bottom up traversal of node tree */
    int nid, stayprime, staysec;
    dint putprimemax, floating, putprime, room, getsec=0;
    struct S_node *parent;
    node = NodeQueue[i];
    nid = node->id;
    putprimemax = node->get.prime;      /* 0 for leaf, or max number of p's
                                           children can accommodate for non
                                           leaf */
    stayprime = node->stay.prime;
    staysec = node->stay.sec;
    parent = node->parent;
    floating = TotalPGlobal[nid] - NOfPToStay[nid];
    putprime = putprimemax - floating;
    if (putprime>stayprime)  putprime = stayprime;
    node->get.prime = -putprime;
    room = nlpmax - (stayprime + staysec - putprime);
    if (room<0) {                       /* have to put some secondaries to get
                                           primaries */
      node->get.sec = getsec = room;    /* getsec is negative to mean to put */
      if (room+staysec<0) return(FALSE);
                                        /* nlpmax<stay.prime+getprime */
    }
    if (!parent) break;
    if (getsec) {                       /* getsec is negative to mean to put
                                           and thus number of parant's to-stay
                                           is decremented to make its getprime
                                           larger in the result */
      NOfPToStay[node->parentid] += getsec;
    } else {                            /* getsec is 0 to mean the node has
                                           some room to get secondaries */
      parent->get.prime += room;
    }
  }
  Verbose(2,vprint("try_stable=TRUE"));

  for (i=nn-1; i>=0; i--) {             /* top down traversal of node tree */
    int nid, k, npfrac, incgp;
    dint nproot, floating, nptotal, npave;
    node = NodeQueue[i];
    if (!(ch=node->child))  continue;   /* a leaf may reside below some non-
                                           leaves in NodeQueue when its number
                                           of primaries is equal to the
                                           average */
    nid = node->id;
    floating = TotalPGlobal[nid] - NOfPToStay[nid];
                                        /* # of transboundaries + overflows */
    nproot = node->stay.prime + node->stay.sec + node->get.sec;
    if (nproot>nlpmax) {                /* secondary assignment made primary
                                           overflow */
      dint getprime = nlpmax - nproot;  /* getprime<0 to mean to put */
      node->get.prime = getprime;
      floating -= getprime;
      nproot = nlpmax;
    } else {
      node->get.prime = 0;
    }
    if (floating==0)  continue;
    incgp = 0;
    if ((nptotal=assign_particles(nproot, floating, ch, 0, &k))<0) {
                                        /* try to avoid moving primaries of
                                           children to their children */
      incgp = 1;
      if ((nptotal=assign_particles(nproot, floating, ch, 1, &k))<0)
                                        /* allow moving primaries of
                                           children to their children */
        errstop("SECONDARY PARTICLE ASSIGNMENT STABILITY CHECK ERROR");
    }
    npave = nptotal / k;
    npfrac = nptotal - npave*k;         /* should be faster than nptotal%k */
    if (npfrac==0) {
      npave--;  npfrac = k;
    }                                   /* npave = ceil(average)-1 */
    if (nproot<=npave) {
      if (npfrac-- > 0)  nproot--;      /* get to have ceil(average) */
      node->get.prime = npave - nproot;
    }
    for (ch=node->child; ch; ch=ch->sibling) {
      int npch = ch->stay.prime + ch->stay.sec + ch->get.sec;
      int gp = ch->get.prime;
      if (gp>0 || incgp)  npch += gp;
      if (npch<=npave) {
        if (npfrac-- > 0)  npch--;
        ch->get.sec = npave - npch;
      }
    }
  }
  if (Special_Pexc_Sched(level)) return(TRUE);
  schedule_particle_exchange(currmode==MODE_NORM_SEC ? 0 : -1);
  make_comm_count(currmode, level, 0, Nodes[myRank].parentid, stats);
  return(TRUE);
}
static void
count_stay() {
  int nn=nOfNodes, ns=nOfSpecies, me=myRank;
  int *np, *stay=TempArray;
  struct S_node *node;
  int i, s, sec;

  stay[me] = 0;
  for (s=0,np=NOfPLocal; s<ns; s++,np+=nn)  stay[me] += np[me];
                                                /* NOfPLocal[0][s][me] */
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, stay, 1, MPI_INT, MCW);
  for (i=0,node=Nodes; i<nn; i++,node++) {
      node->stay.prime = NOfPToStay[i] = stay[i];
      node->get.prime = 0;
  }
  sec = RegionId[1];
  stay[me] = 0;
  if (sec>=0)                                   /* np=&NOfPLocal[1][0][0] */
    for (s=0; s<ns; s++,np+=nn)  stay[me] += np[sec];
                                                /* NOfPLocal[1][s][me] */
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, stay, 1, MPI_INT, MCW);
  for (i=0,node=Nodes; i<nn; i++,node++) {
    sec = node->parentid;
    if (sec>=0)  NOfPToStay[sec] += (node->stay.sec = stay[i]);
    else         node->stay.sec = 0;
    node->get.sec = 0;
  }
}
static dint
assign_particles(dint npr, dint npt, struct S_node *ch, int incgp, int *nget) {
  int *np=TempArray;            /* used just for temporary sorting buffer */
  int n, i;
  dint nlpmax = nOfLocalPMax;

  np[0] = npr;
  for (n=1; ch; ch=ch->sibling, n++) {
    int gp=ch->get.prime;
    np[n] = ch->stay.prime + ch->stay.sec + ch->get.sec;
    if (gp>0 || incgp)  np[n] += gp;
  }
  qsort(np, n, sizeof(int), compare_int);       /* sort ascendingly */
  for (i=0; i<n; i++) {
    dint npc=np[i];
    if (npt<=npc*i)  break;
    npt += npc;
  }
  *nget = i;
  return(npt>nlpmax*i ? -1 : npt);
}
static int
compare_int(const void* x, const void* y) {
  int xx=*((int*)x), yy=*((int*)y);

  if (xx<yy)  return(-1);
  if (xx>yy)  return(1);
  return(0);
}
static void
schedule_particle_exchange(int reb) {
  int me=myRank, nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns;
  struct S_node *mynode, *ch;
  int i, slidx;
  struct S_commsched_context context;

  RLIndex[0] = 0;
  context.neighbor = 0;
  context.sender = (reb<0 ||reb==3) ? 0 : DstNeighbors[0];
  context.comidx = 0;
  context.spec = 0;  context.dones = 0;  context.donen = 0;
  for (i=0; i<nn; i++) TempArray[i] = 0;
  TempArray[context.sender] = 1;
  if (reb>0) {
    mynode = NodesNext + me;
    for (ch=mynode->child; ch; ch=ch->sibling)
      ch->get.sec -= (ch->stay.sec=count_real_stay(NOfPrimaries+nnns+ch->id));
                                                /* NOfPrimaries[1][0][cid] */
    sched_comm(mynode->get.prime, me, 0, reb, &context);
    for (ch=mynode->child; ch; ch=ch->sibling)
      sched_comm(ch->get.sec, ch->id, ns, reb, &context);
    if (mynode->parent)
      mynode->get.sec -=
        (mynode->stay.sec=count_real_stay(NOfPLocal+nnns+mynode->parentid));
                                                /* NOfPLocal[1][0][pid] */
  } else {
    mynode = Nodes + me;
    sched_comm(mynode->get.prime, me, 0, reb, &context);
    for (ch=mynode->child; ch; ch=ch->sibling)
      sched_comm(ch->get.sec, ch->id, ns, reb, &context);
    if (reb<0)  reb = 3;
  }
  SLHeadTail[0] = slidx = context.comidx;
  if (reb==3) return;

  for (i=context.neighbor+1; i<=OH_NEIGHBORS; i++)  RLIndex[i] = slidx;
  for (i=0; i<OH_NEIGHBORS; i++) {
    int dst=DstNeighbors[i];
    int src=SrcNeighbors[i];
    int rc;
    MPI_Status st;
    if (dst==me) continue;
    if (src>=0) {
      if (dst>=0)
        MPI_Sendrecv(CommList+RLIndex[i], RLIndex[i+1]-RLIndex[i], T_Commlist,
                     dst, 0,
                     CommList+slidx, nn+nnns, T_Commlist, src, 0, MCW, &st);
      else
        MPI_Recv(CommList+slidx, nn+nnns, T_Commlist, src, 0, MCW, &st);
      MPI_Get_count(&st, T_Commlist, &rc);
      slidx += rc;
    } else if (dst>=0)
      MPI_Send(CommList+RLIndex[i], RLIndex[i+1]-RLIndex[i], T_Commlist,
               dst, 0, MCW);
  }
  SLHeadTail[1] = slidx;
  if (reb==2) return;

  SecSLHeadTail[0] = SecSLHeadTail[1] = 0;
  oh1_broadcast(SLHeadTail, SecSLHeadTail, 2, 2, MPI_INT, MPI_INT);
  oh1_broadcast(CommList, CommList+slidx, slidx, SecSLHeadTail[1],
                T_Commlist, T_Commlist);
}
static int
count_real_stay(int *np) {
  const int ns=nOfSpecies, nn=nOfNodes;
  int stay, s;

  for (s=0,stay=0; s<ns; s++,np+=nn)  stay += *np;
  return(stay);
}
static void
sched_comm(int toget, int rid, int tag, int reb,
           struct S_commsched_context *context) {
  int neighbor = context->neighbor, sid = context->sender;
  int comidx = context->comidx;
  int s=context->spec, havedones=context->dones, havedonen=context->donen;
  int me=myRank;
  int nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns;
  int i=nn*s+sid;                       /* [0][s][sid] */
  struct S_node *nodesnext = reb>0 ? NodesNext : Nodes;

  while (toget>0) {
    struct S_node *snoden=nodesnext+sid;
    int npp = NOfPrimaries[i], nps = NOfPrimaries[i+nnns];
                                                        /* [0/1][s][sid] */
    int toput, hdninc = 0;
    int next=1;
    if (sid<0)  local_errstop("PARTICLE TRANSFER SCHEDULING ERROR");
    if (sid==me) {
      int nput = snoden->get.prime + havedonen;
      nput = nput<0 ? -nput : 0;
      if (nput<npp) npp = nput;
      hdninc = npp;
    }
    else if (snoden->parentid==me) {
      int nput = snoden->get.sec + havedonen;
      nput = nput<0 ? -nput : 0;
      if (nput<nps) nps = nput;
      hdninc = nps;
    }
    toput = npp + nps - havedones;
    if (toput>0) {
      struct S_commlist *cptr = CommList+(comidx++);
      if (toput>toget) {
        havedones += toget;
        toput = toget;
        next = 0;
      }
      cptr->rid = rid;  cptr->sid = sid;  cptr->region = me;
      cptr->count = toput;  cptr->tag = tag + s;
      toget -= toput;
    }
    if (next) {
      havedones = 0;  havedonen += hdninc;
      s++;  i += nn;
      if (s==ns) {
        havedonen = 0;
        s = 0;
        if (reb>=0 && reb!=3) {
          struct S_node *nodes = (neighbor<OH_NEIGHBORS ? Nodes : NodesNext);
          struct S_node *snode = nodes + sid;
          while (sid>=0) {
            if (neighbor==OH_NEIGHBORS) {
              snode = snode->sibling;  sid = snode ? snode - nodes : -1;
            }
            else if (sid==DstNeighbors[neighbor] && reb<2 && snode->child) {
              snode = snode->child;  sid = snode - Nodes;
            }
            else if (sid!=DstNeighbors[neighbor] && reb<2 && snode->sibling) {
              snode = snode->sibling;  sid = snode - Nodes;
            }
            else {
              RLIndex[++neighbor] = comidx;
              while(neighbor<OH_NEIGHBORS && (sid=DstNeighbors[neighbor])<0)
                RLIndex[++neighbor] = comidx;
              if (neighbor==OH_NEIGHBORS) {
                nodes = NodesNext;
                snode = nodes[me].child;
                sid = (snode && reb) ? snode - nodes : -1;
              } else {
                snode = Nodes + sid;
              }
            }
            if (sid>=0 && TempArray[sid]==0) {
              TempArray[sid] = 1;  break;
            }
          }
        } else {
          sid++;
        }
        i = sid;
      }
    }
  }
  context->neighbor = neighbor;  context->sender = sid;
  context->comidx = comidx;
  context->spec = s;  context->dones = havedones;  context->donen = havedonen;
}
static void
make_comm_count(int currmode, int level, int reb, int oldparent, int stats) {
  int nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns, nnns2=nnns*2, me=myRank;
  struct S_node *mynode=Nodes+me;
  int newparent=mynode->parentid;
  int ps, s, i;

  if (reb || currmode==MODE_ANY_SEC) {
    SecRLSize = 0;
    oh1_broadcast(SLHeadTail, &SecRLSize, 1, 1, MPI_INT, MPI_INT);
    if (currmode==MODE_NORM_PRI)
      SecRList = CommList + SLHeadTail[1];
    else if (currmode==MODE_NORM_SEC)
      SecRList = CommList + SLHeadTail[1] + SecSLHeadTail[1];
    else
      SecRList = CommList + SLHeadTail[0];
    oh1_broadcast(CommList, SecRList, SLHeadTail[0], SecRLSize,
                  T_Commlist, T_Commlist);
  } else {
    SecRList = CommList + SLHeadTail[1];
    SecRLSize = SecSLHeadTail[0];
  }
  for (s=0; s<ns*2; s++) TotalPNext[s] = 0;             /* TotalPNext[p][s] */
  if (level==1 || Mode_Is_Any(currmode) || stats) {
    for (i=0; i<nnns2; i++)  NOfRecv[i] = NOfSend[i] = 0;
    make_recv_count(CommList, SLHeadTail[0]);
    if (newparent>=0)
      make_recv_count(SecRList, SecRLSize);
    if (oldparent!=newparent && oldparent>=0)
      make_recv_count(CommList+SLHeadTail[1], SecSLHeadTail[0]);
    if (Mode_Is_Any(currmode)) {
      MPI_Alltoall(NOfRecv, 1, T_Histogram, NOfSend, 1, T_Histogram, MCW);
#ifndef INTEL_MPI_BUG_FIXED
      for (ps=0,i=me; ps<2; ps++)  for (s=0; s<ns; s++,i+=nn)
        NOfSend[i] = NOfRecv[i];
#endif
    } else {
      make_send_count(CommList+SLHeadTail[0], SLHeadTail[1]-SLHeadTail[0]);
      if (oldparent>=0)
        make_send_count(CommList+SLHeadTail[1]+SecSLHeadTail[0],
                        SecSLHeadTail[1]-SecSLHeadTail[0]);
    }
  } else {
    count_next_particles(CommList, SLHeadTail[0]);
    if (newparent>=0)
      count_next_particles(SecRList, SecRLSize);
  }
  if (stats) stats_secondary_comm(currmode, reb);
  if (level==1) {
    for (ps=0,i=0; ps<(newparent<0?1:2); ps++) {
      int putme = ps==0 ? -mynode->get.prime : -mynode->get.sec;
      int *mynps = ps==0 ? NOfPLocal+me : NOfPLocal+nnns+newparent;
      if (putme<0) putme = 0;
      for (s=0; s<ns; s++,i++,mynps+=nn) {
        int stay=*mynps;
        int tpni=TotalPNext[i];
        if (putme<stay) {
          TotalPNext[i] = tpni + stay - putme;  putme = 0;
        }
        else putme -= stay;
      }
    }
  }
}
static void
make_recv_count(struct S_commlist* rlist, int rlsize) {
  int me=myRank, nn=nOfNodes;
  int i;

  for (i=0; i<rlsize; i++) {
    int rid=rlist[i].rid, sid=rlist[i].sid;
    int tag=rlist[i].tag, count=rlist[i].count;
    if (rid==me) {
      NOfRecv[tag*nn+sid] = count;
      TotalPNext[tag] += count;
    }
    if (sid==me)
      NOfSend[tag*nn+rid] = count;
  }
}
static void
make_send_count(struct S_commlist* slist, int slsize) {
  int me=myRank, nn=nOfNodes;
  int i;

  for (i=0; i<slsize; i++) {
    if (slist[i].sid==me)
      NOfSend[slist[i].tag*nn+slist[i].rid] = slist[i].count;
  }
}
static void
count_next_particles(struct S_commlist* rlist, int rlsize) {
  int me=myRank, i;

  for (i=0; i<rlsize; i++) {
    if (rlist[i].rid==me)
      TotalPNext[rlist[i].tag] += rlist[i].count;
  }
}
void
oh1_broadcast_(void* pbuf, void* sbuf, int *pcount, int *scount,
               int *ptype, int *stype) {
  oh1_broadcast(pbuf, sbuf, *pcount, *scount,
                MPI_Type_f2c(*ptype), MPI_Type_f2c(*stype));
}
void
oh1_broadcast(void* pbuf, void* sbuf, int pcount, int scount,
              MPI_Datatype ptype, MPI_Datatype stype) {

  if (MyComm->black) {
    if (MyComm->prime!=MPI_COMM_NULL && pcount)
      MPI_Bcast(pbuf, pcount, ptype, MyComm->rank, MyComm->prime);
    if (MyComm->sec!=MPI_COMM_NULL && scount)
      MPI_Bcast(sbuf, scount, stype, MyComm->root, MyComm->sec);
  } else {
    if (MyComm->sec!=MPI_COMM_NULL && scount)
      MPI_Bcast(sbuf, scount, stype, MyComm->root, MyComm->sec);
    if (MyComm->prime!=MPI_COMM_NULL && pcount)
      MPI_Bcast(pbuf, pcount, ptype, MyComm->rank, MyComm->prime);
  }
}
void
rebalance1(int currmode, int level, int stats) {
  int nn=nOfNodes;
  dint nofp=nOfParticles;
  dint npavefloor=nofp/nn;
  dint npfracin=nofp-npavefloor*nn, npfracout=npfracin;
  dint npavein=npavefloor+(npfracin==0 ? 0 : 1), npaveout=npavein;
  int ns=nOfSpecies;
  int i, j, k, s, bot, pm=Mode_PS(currmode)-1, me=myRank;
  struct S_node *node, *mynode=NodesNext+me, *root;

  if (stats) oh1_stats_time(STATS_REBALANCE, 0);
  Verbose(2,vprint("rebalance"));

  LessHeap.n = GreaterHeap.n = 0;
  for (i=0; i<nn; i++) GreaterHeap.index[i] = 0;

  for (i=0,bot=0,node=NodesNext; i<nn; i++,node++) {
    dint npg=TotalPGlobal[i];
    if (npg<npavein) {
      if (--npfracin==0) npavein--;
      push_heap(i, &LessHeap, 0);
      NodeQueue[bot++] = node;
    } else {
      push_heap(i, &GreaterHeap, 1);
    }
    *node = Nodes[i];
    node->child = NULL;
    if (pm) node->parentid = -1;
  }
  while (LessHeap.n) {
    struct S_node *parent;
    dint npg;
    int get, pid, h;
    j = pop_heap(&LessHeap, 0);
    node = NodesNext + j;
    get = npaveout - TotalPGlobal[j];
    if (--npfracout==0) npaveout--;
    if ((k=node->parentid)>=0 && (h=GreaterHeap.index[k]))
      remove_heap(&GreaterHeap, 1, h);
    else
      k = pop_heap(&GreaterHeap, 1);
    node->get.sec = get;
    parent = NodesNext + k;
    node->parentid = k;  node->parent = parent;
    node->sibling = parent->child;
    parent->child = node;
    npg = (TotalPGlobal[k] -= get);
    if (npg<npavein) {
      if (--npfracin==0) npavein--;
      push_heap(k, &LessHeap, 0);  NodeQueue[bot++] = parent;
    } else {
      push_heap(k, &GreaterHeap, 1);
    }
  }
  root = NodesNext + GreaterHeap.node[1];
  root->parentid = -1;  root->parent = root->sibling = NULL;
  root->get.sec = 0;
  k = root->id;
  for (i=2; i<=GreaterHeap.n; i++) {
    j = GreaterHeap.node[i];
    node = NodesNext + j;
    node->get.sec = 0;
    node->parentid = k;  node->parent = root;
    node->sibling = root->child;  root->child = node;
    NodeQueue[bot++] = node;
  }
  NodeQueue[bot] = root;

  mynode->get.prime = TotalPGlobal[me] -
                      (mynode->stay.prime=count_real_stay(NOfPLocal+me));
  if (Special_Pexc_Sched(level)) return;
  schedule_particle_exchange(currmode==MODE_NORM_SEC ?
                             1 : (currmode==MODE_NORM_PRI ? 2 : 3));
  build_new_comm(currmode, level, 1, stats);
}
void
build_new_comm(int currmode, int level, int nbridx, int stats) {
  int bot=nOfNodes-1, me=myRank;
  struct S_node *node, *ch;
  struct S_node *mynode=NodesNext+me, *root=NodeQueue[bot];
  int oldparent=Mode_PS(currmode) ? Nodes[me].parentid : -1;
  int i, j;
  MPI_Group grpw=GroupWorld, grp;

  node = Nodes;  Nodes = NodesNext;  NodesNext = node;

  if (stats) oh1_stats_time(STATS_REB_COMM, 0);
  for (i=0; i<Comms.n; i++) {
    if (Comms.body[i] != MPI_COMM_NULL)
      MPI_Comm_free(Comms.body+i);
  }
  root->comm.black = 0;  root->comm.sec = -1;
  for (i=0; bot>=0; bot--) {
    int black, rid;
    node = NodeQueue[bot];
    if (!(ch=node->child))  continue;   /* a leaf may reside below some non-
                                           leaves in NodeQueue when its number
                                           of primaries is equal to the
                                           average */
    node->comm.prime = i;
    rid = TempArray[0] = node->id;
    black = 1 - node->comm.black;
    for (j=1; ch; ch=ch->sibling, j++) {
      TempArray[j] = ch->id;
      ch->comm.prime = -1;  ch->comm.sec = i;  ch->comm.black = black;
    }
    MPI_Group_incl(grpw, j, TempArray, &grp);
    MPI_Group_translate_ranks(grpw, 1, &rid, grp, &(node->comm.rank));
    MPI_Comm_create(MCW, grp, Comms.body+i);
    MPI_Group_free(&grp);
    i++;
  }
  Comms.n = i;

  if (FamIndex) {
    int *fidx = FamIndex,  *fmem = FamMembers;
    int nn = nOfNodes, j;
    for (i=0,j=0; i<nn; i++) {
      fidx[i] = j;
      fmem[j++] = i;
      for (ch=Nodes[i].child; ch; ch=ch->sibling, j++)  fmem[j] = ch->id;
    }
    fidx[nn] = j;  fmem[j] = root->id;
  }
  MyComm->prime =
    mynode->comm.prime<0 ? MPI_COMM_NULL : Comms.body[mynode->comm.prime];
  MyComm->sec =
    mynode->comm.sec<0 ? MPI_COMM_NULL : Comms.body[mynode->comm.sec];
  MyComm->rank = mynode->comm.prime<0 ? -1 : mynode->comm.rank;
  MyComm->black = mynode->comm.black;
  if ((node=mynode->parent))  MyComm->root = node->comm.rank;
  else MyComm->root = -1;
  if (MyCommC) *MyCommC = *MyComm;
  if (MyCommF) {
    MyCommF->prime = MPI_Comm_c2f(MyComm->prime);
    MyCommF->sec   = MPI_Comm_c2f(MyComm->sec);
    MyCommF->rank  = MyComm->rank;
    MyCommF->root  = MyComm->root;
    MyCommF->black = MyComm->black;
  }
  oh1_broadcast(Neighbors[0], Neighbors[nbridx], OH_NEIGHBORS, OH_NEIGHBORS,
                MPI_INT, MPI_INT);
  if (NeighborsShadow) {
    int (*nb)[OH_NEIGHBORS] = NeighborsShadow;
    for (i=0; i<OH_NEIGHBORS; i++)  nb[2][i] = nb[1][i];
    oh1_broadcast(nb[0], nb[1], OH_NEIGHBORS, OH_NEIGHBORS, MPI_INT, MPI_INT);
  }
  SubdomainId[1] = RegionId[1] = mynode->parentid;

  if (!Special_Pexc_Sched(level))
    make_comm_count(currmode, level, 1,
                    (Mode_Is_Norm(currmode) ? oldparent : -1), stats);
}
static void
push_heap(int r, struct S_heap* heap, int greater) {
  int n=heap->n, *hnode=heap->node, *index=heap->index;
  dint np=TotalPGlobal[r];
  int m, q, g;

  heap->n = ++n;
  for (; n>1; n=m) {
    m = n>>1;  q = hnode[m];
    g = (np>TotalPGlobal[q]) ? 1 : 0;
    if (g!=greater) break;
    hnode[n] = q;  index[q] = n;
  }
  hnode[n] = r;  index[r] = n;
}
static int
pop_heap(struct S_heap* heap, int greater) {
  int pop=heap->node[1];

  remove_heap(heap, greater, 1);
  return(pop);
}
static void
remove_heap(struct S_heap* heap, int greater, int rem) {
  int n=heap->n, *hnode=heap->node, *index=heap->index;
  int id=hnode[n];
  dint np=TotalPGlobal[id];
  int i;

  heap->n = --n;  index[hnode[rem]] = 0;
  if (rem>n) return;
  for (i=rem; ; ) {
    int left=(i<<1), right=left+1;
    if (right<=n) {
      int lid=hnode[left], rid=hnode[right];
      dint lnp=TotalPGlobal[lid], rnp=TotalPGlobal[rid];
      int cgl=(np>lnp)?1:0, cgr=(np>rnp)?1:0, lgr=(lnp>rnp)?1:0;
      if (cgl==greater) {
        if (cgr==greater) {
          hnode[i] = id;  index[id] = i;  return;
        } else {
          hnode[i] = rid;  index[rid] = i;  i = right;
        }
      } else if (lgr==greater) {
        hnode[i] = lid;  index[lid] = i;  i = left;
      } else {
        hnode[i] = rid;  index[rid] = i;  i = right;
      }
    } else {
      if (left<=n) {
        int lid=hnode[left];
        int cgl=(np>TotalPGlobal[lid])?1:0;
        if (cgl==greater) {
          hnode[i] = id;  index[id] = i;
        } else {                /* we know left node has no children. */
          hnode[i] = lid;  index[lid] = i;
          hnode[left] = id;  index[id] = left;
        }
      } else {
        hnode[i] = id;  index[id] = i;
      }
      return;
    }
  }
}
int
oh1_accom_mode_() {
  return(accMode);
}
int
oh1_accom_mode() {
  return(accMode);
}
void
oh1_all_reduce_(void *pbuf, void *sbuf, int *pcount, int *scount,
                int *ptype, int *stype, int *pop, int *sop) {
  oh1_all_reduce(pbuf, sbuf, *pcount, *scount,
                 MPI_Type_f2c(*ptype), MPI_Type_f2c(*stype),
                 MPI_Op_f2c(*pop), MPI_Op_f2c(*sop));
}
void
oh1_all_reduce(void *pbuf, void *sbuf, int pcount, int scount,
               MPI_Datatype ptype, MPI_Datatype stype,
               MPI_Op pop, MPI_Op sop) {

  if (MyComm->black) {
    if (MyComm->prime!=MPI_COMM_NULL)
      MPI_Allreduce(MPI_IN_PLACE, pbuf, pcount, ptype, pop, MyComm->prime);
    if (MyComm->sec!=MPI_COMM_NULL)
      MPI_Allreduce(MPI_IN_PLACE, sbuf, scount, stype, sop, MyComm->sec);
  } else {
    if (MyComm->sec!=MPI_COMM_NULL)
      MPI_Allreduce(MPI_IN_PLACE, sbuf, scount, stype, sop, MyComm->sec);
    if (MyComm->prime!=MPI_COMM_NULL)
      MPI_Allreduce(MPI_IN_PLACE, pbuf, pcount, ptype, pop, MyComm->prime);
  }
}
void
oh1_reduce_(void *pbuf, void *sbuf, int *pcount, int *scount,
            int *ptype, int *stype, int *pop, int *sop) {
  oh1_reduce(pbuf, sbuf, *pcount, *scount,
             MPI_Type_f2c(*ptype), MPI_Type_f2c(*stype),
             MPI_Op_f2c(*pop), MPI_Op_f2c(*sop));
}
void
oh1_reduce(void *pbuf, void *sbuf, int pcount, int scount,
           MPI_Datatype ptype, MPI_Datatype stype, MPI_Op pop, MPI_Op sop) {

  if (MyComm->black) {
    if (MyComm->prime!=MPI_COMM_NULL)
      MPI_Reduce(MPI_IN_PLACE, pbuf, pcount, ptype, pop, MyComm->rank,
                 MyComm->prime);
    if (MyComm->sec!=MPI_COMM_NULL)
      MPI_Reduce(sbuf, NULL, scount, stype, sop, MyComm->root, MyComm->sec);
  } else {
    if (MyComm->sec!=MPI_COMM_NULL)
      MPI_Reduce(sbuf, NULL, scount, stype, sop, MyComm->root, MyComm->sec);
    if (MyComm->prime!=MPI_COMM_NULL)
      MPI_Reduce(MPI_IN_PLACE, pbuf, pcount, ptype, pop, MyComm->rank,
                 MyComm->prime);
  }
}
void
oh1_init_stats_(int *key, int *ps) {
  if (statsMode) oh1_init_stats(*key, *ps);
}
void
oh1_init_stats(int key, int ps) {
  int i;

  if (!statsMode)  return;
  clear_stats(&Stats.subtotal);
  clear_stats(&Stats.total);
  Stats.curr.time.key = (key<<1) + ps;
  Stats.curr.time.value = MPI_Wtime();
  for (i=0; i<(STATS_TIMINGS<<1)+1; i++) Stats.curr.time.ev[i] = 0;
  MPI_Type_contiguous(sizeof(struct S_statstime), MPI_BYTE, &T_StatsTime);
  MPI_Type_commit(&T_StatsTime);
  MPI_Op_create(stats_reduce_time, 1, &Op_StatsTime);
  MPI_Op_create(stats_reduce_part, 1, &Op_StatsPart);
}
static void
clear_stats(struct S_statstotal *stotal) {
  int i;
  struct S_statstime *st = stotal->time;
  struct S_statspart *sp = stotal->part;

  for (i=0; i<STATS_TIMINGS<<1; i++) {
    st[i].ev = 0;
    st[i].min = DBL_MAX;
    st[i].max = 0.0;
    st[i].total = 0.0;
  }
  for (i=0; i<STATS_PARTS; i++) {
    sp[i].min = INT_MAX;
    sp[i].max = 0;
    sp[i].total = 0;
  }
  sp[STATS_PART_PRIMARY].min = sp[STATS_PART_SECONDARY].min = 0;
}
void
oh1_stats_time_(int *key, int *ps) {
  if (statsMode)  oh1_stats_time(*key, *ps);
}
void
oh1_stats_time(int key, int ps) {
  double t;
  int k=Stats.curr.time.key;

  if (!statsMode)  return;
  t = MPI_Wtime();
  Stats.curr.time.val[k] = t - Stats.curr.time.value;
  Stats.curr.time.ev[k] = 1;
  Stats.curr.time.value = t;
  Stats.curr.time.key = (key<<1) + ps;
}
static void
stats_primary_comm(int currmode) {
  stats_comm(NOfPrimaries, NOfPLocal, Stats.curr.part, nOfSpecies*2);
  Stats.curr.part[STATS_PART_PRIMARY] =
    (currmode==MODE_ANY_PRI) ? 3 : Mode_PS(currmode)+1;
}
static void
stats_secondary_comm(int currmode, int reb) {
  int ns=nOfSpecies, nnns=nOfNodes*ns;

  stats_comm(NOfRecv, NOfSend, Stats.curr.part, ns);
  stats_comm(NOfRecv+nnns, NOfSend+nnns,
             Stats.curr.part+STATS_PART_MOVE_SEC_MIN, ns);
  Stats.curr.part[STATS_PART_SECONDARY] =
    Mode_PS(currmode) ? (reb ? 3 : 2) : 1;
}
static void
stats_comm(int* nrecv, int* nsend, dint* scp, int ns) {
  int i, s, nn=nOfNodes, me=myRank;
  int get=0, put=0, minmove=INT_MAX, maxmove=0, nmove=0;

  for (i=0; i<nn; i++) {
    int g=0, p=0, *npr=nrecv, *nps=nsend;
    for (s=0; s<ns; s++,npr+=nn,nps+=nn) {
      g += nrecv[i];                            /* nrecv[s][i] */
      p += nsend[i];                            /* nsend[s][i] */
    }
    if (i!=me) {
      get += g;  put += p;
      if (g>0) {
        if (g<minmove) minmove = g;
        if (g>maxmove) maxmove = g;
        nmove++;
      }
    }
    if (minmove>maxmove) minmove = 0;
    scp[STATS_PART_MOVE_PRI_MIN] = minmove;
    scp[STATS_PART_MOVE_PRI_MAX] = maxmove;
    scp[STATS_PART_MOVE_PRI_AVE] = nmove;
    scp[STATS_PART_GET_PRI_MIN] = scp[STATS_PART_GET_PRI_MAX]
                                = scp[STATS_PART_PG_PRI_AVE] = get;
    scp[STATS_PART_PUT_PRI_MIN] = scp[STATS_PART_PUT_PRI_MAX] = put;
  }
}
void
oh1_show_stats_(int *step, int *currmode) {
  if (statsMode)  oh1_show_stats(*step, *currmode);
}
void
oh1_show_stats(int step, int currmode) {

  if (!statsMode)  return;
  oh1_stats_time(STATS_TIMINGS,0);
  MPI_Barrier(MCW);
  if (statsMode==2) {
    update_stats(&Stats.subtotal, step, currmode);
    if (step%reportIteration == 0) {
      print_stats(&Stats.subtotal, step, reportIteration);
      clear_stats(&Stats.subtotal);
    }
  }
  update_stats(&Stats.total, step, currmode);
  MPI_Barrier(MCW);
}
#define Round(NUM,DEN) (DEN ? (NUM+(DEN>>1))/DEN : 0)

static void
update_stats(struct S_statstotal *stotal, int step, int currmode) {
  int i, j, k, ev, nn=nOfNodes;
  int evclr = stotal==&Stats.total, reduce = statsMode==1 || !evclr;
  struct S_statstime *st = stotal->time;
  struct S_statspart *sp = stotal->part;
  int pm=Mode_PS(currmode)-1;
  int transkey=pm?STATS_PART_PRIMARY:STATS_PART_SECONDARY;
  dint trans=Stats.curr.part[transkey];

  for (i=0; i<STATS_TIMINGS<<1; i++) {
    if ((ev=Stats.curr.time.ev[i])) {
      double t = Stats.curr.time.val[i];
      st[i].ev++;
      if (t<st[i].min) st[i].min = t;
      if (t>st[i].max) st[i].max = t;
      st[i].total += t;
      if (evclr) Stats.curr.time.ev[i] = 0;
    }
  }
  if (step<=0) return;
  if (myRank!=0) {
    if (reduce)
      MPI_Reduce(st, NULL, STATS_PART_PRIMARY, MPI_LONG_LONG_INT, Op_StatsPart,
                 0, MCW);
    return;
  }
  if (reduce)
    MPI_Reduce(MPI_IN_PLACE, st, STATS_PART_PRIMARY, MPI_LONG_LONG_INT,
               Op_StatsPart, 0, MCW);
  for (i=0,j=0; i<(pm?1:2); i++,j+=STATS_PART_MOVE_SEC_MIN) {
    dint *scp=Stats.curr.part+j;
    struct S_statspart *spps=sp+j;
    if (reduce) {
      scp[STATS_PART_MOVE_PRI_AVE] = Round(scp[STATS_PART_PG_PRI_AVE],
                                           scp[STATS_PART_MOVE_PRI_AVE]);
      scp[STATS_PART_PG_PRI_AVE] = Round(scp[STATS_PART_PG_PRI_AVE], nn);
    }
    for (k=0; k<STATS_PART_MOVE_SEC_MIN; k++) {
      dint n = scp[k];
      if (n<spps[k].min) spps[k].min = n;
      if (n>spps[k].max) spps[k].max = n;
      spps[k].total += n;
    }
  }
  if      (trans==1) sp[transkey].min++;
  else if (trans==2) sp[transkey].max++;
  else if (trans==3) sp[transkey].total++;
}
#define Stats_Reduce_Part_Min(KEY) { if (io[KEY]<in[KEY]) io[KEY] = in[KEY]; }
#define Stats_Reduce_Part_Max(KEY) { if (io[KEY]>in[KEY]) io[KEY] = in[KEY]; }
#define Stats_Reduce_Part_Sum(KEY) { io[KEY] += in[KEY]; }

static void
stats_reduce_part(void* inarg, void* ioarg, int* len, MPI_Datatype* type) {
  dint *in=(dint*)inarg, *io=(dint*)ioarg;
  int ps, statsbase=0;

  for (ps=0; ps<2; ps++,statsbase+=STATS_PART_MOVE_SEC_MIN) {
    Stats_Reduce_Part_Min(statsbase+STATS_PART_MOVE_PRI_MIN);
    Stats_Reduce_Part_Max(statsbase+STATS_PART_MOVE_PRI_MAX);
    Stats_Reduce_Part_Sum(statsbase+STATS_PART_MOVE_PRI_AVE);
    Stats_Reduce_Part_Min(statsbase+STATS_PART_GET_PRI_MIN);
    Stats_Reduce_Part_Max(statsbase+STATS_PART_GET_PRI_MAX);
    Stats_Reduce_Part_Min(statsbase+STATS_PART_PUT_PRI_MIN);
    Stats_Reduce_Part_Max(statsbase+STATS_PART_PUT_PRI_MAX);
    Stats_Reduce_Part_Sum(statsbase+STATS_PART_PG_PRI_AVE);
  }
}
static void
print_stats(struct S_statstotal *stotal, int step, int nstep) {
  int i;
  struct S_statstime *st = stotal->time;
  struct S_statspart *sp = stotal->part;

  if (myRank!=0) {
    MPI_Reduce(st, NULL, STATS_TIMINGS<<1, T_StatsTime, Op_StatsTime, 0, MCW);
    return;
  }
  MPI_Reduce(MPI_IN_PLACE, st, STATS_TIMINGS<<1, T_StatsTime, Op_StatsTime,
             0, MCW);
  printf("\n");
  if (stotal==&Stats.subtotal)
    printf("# Subtotal Statistics for %d..%d\n",  step-nstep+1, step);
  else
    printf("# Total Statistics\n");
  printf("## Execution Times (sec)\n");
  for (i=0; i<STATS_TIMINGS<<1; i++) {
    if (st[i].ev==0)
      printf("  %-29s = -------- / -------- / -------- / ------------\n",
             StatsTimeStrings[i]);
    else
      printf("  %-29s = %8.3f / %8.3f / %8.3f / %12.3f\n",
             StatsTimeStrings[i],
             st[i].min, st[i].max, st[i].total/st[i].ev, st[i].total);
  }
  printf("## Particle Movements\n");
  for (i=0; i<STATS_PARTS; i++) {
    if (i<STATS_PART_PRIMARY && sp[i].min>sp[i].max)
      printf("  %-29s = -------- / -------- / -------- / ------------\n",
             StatsPartStrings[i]);
    else if (i<STATS_PART_PRIMARY)
      printf("  %-29s = %8lld / %8lld / %8lld / %12lld\n",
             StatsPartStrings[i],
             sp[i].min, sp[i].max, Round(sp[i].total,nstep), sp[i].total);
    else
      printf("  %-29s = %8lld / %8lld / %8lld / %12lld\n",
             StatsPartStrings[i],
             sp[i].min, sp[i].max, sp[i].total,
             sp[i].min+sp[i].max+sp[i].total);
  }
}
static void
stats_reduce_time(void* inarg, void* ioarg, int* len, MPI_Datatype* type) {
  struct S_statstime *in=(struct S_statstime*)inarg;
  struct S_statstime *io=(struct S_statstime*)ioarg;
  int n=*len, i;

  for (i=0; i<n; i++) {
    if (in[i].ev>0) {
      io[i].ev += in[i].ev;
      if (in[i].min<io[i].min) io[i].min = in[i].min;
      if (in[i].max>io[i].max) io[i].max = in[i].max;
      io[i].total += in[i].total;
    }
  }
}
void
oh1_print_stats_(int *nstep) {
  if (statsMode)  print_stats(&Stats.total, 0, *nstep);
}
void
oh1_print_stats(int nstep) {
  if (statsMode)  print_stats(&Stats.total, 0, nstep);
}
void
oh1_verbose_(char *message) {
  Verbose(1, vprint(message));
}
void
oh1_verbose(char *message) {
  Verbose(1, vprint(message));
}
#define Vprint(FORMAT, RANKFORMAT) {\
  char buf[1024];\
  va_list v;\
  sprintf(buf, RANKFORMAT, myRank);\
  strcat(buf, FORMAT);\
  strcat(buf, "\n");\
  va_start(v, FORMAT);\
  vprintf(buf, v);\
  fflush(stdout);\
  va_end(v);\
}
#define Vprint_Norank(FORMAT, RANKFORMAT) {\
  char buf[1024];\
  va_list v;\
  sprintf(buf, RANKFORMAT);\
  strcat(buf, FORMAT);\
  strcat(buf, "\n");\
  va_start(v, FORMAT);\
  vprintf(buf, v);\
  fflush(stdout);\
  va_end(v);\
}
void
vprint(char* format, ...) {

  if (verboseMode>=3) { Vprint(format, "*Starting[%d] "); }
  else                { Vprint_Norank(format, "*Starting "); }
}
void
dprint(char* format, ...) {

  if (nOfNodes>=1000)     { Vprint(format, "#Debug[%04d] "); }
  else if (nOfNodes>=100) { Vprint(format, "#Debug[%03d] "); }
  else if (nOfNodes>=10)  { Vprint(format, "#Debug[%02d] "); }
  else                    { Vprint(format, "#Debug[%d] "); }
}
