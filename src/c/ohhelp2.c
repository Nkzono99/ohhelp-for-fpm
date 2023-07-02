/* File: ohhelp2.c
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#define EXTERN extern
#include "ohhelp1.h"
#undef  EXTERN
#define EXTERN
#include "ohhelp2.h"

/* Prototypes for private functions. */
static int   try_primary2(int currmode, int level, int stats);
static int   try_stable2(int currmode, int level, int stats);
static void  rebalance2(int currmode, int level, int stats);
static void  move_to_sendbuf_secondary(int secondary, int stats);
static void  move_to_sendbuf_uw(int ps, int me, int *putmes, int cbase,
                                int *ctp, int nbase, int *ntp,
                                struct S_particle **rbb);
static void  move_to_sendbuf_dw(int ps, int me, int *putmes, int ctail,
                                int *ctp, int ntail, int *ntp);
static void  move_injected_to_sendbuf();
static void  move_injected_from_sendbuf(int *injected, int mysd,
                                        struct S_particle **rbb);
static void  receive_particles(struct S_commlist *rlist, int rlsize, int *req);
static void  send_particles(struct S_commlist *slist, int slsize, int myregion,
                            int parentregion, int *req);

void
oh2_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
          int *totalp, struct S_particle *pbuf, int *pbase, int *maxlocalp,
          struct S_mycommf *mycomm, int *nbor, int *pcoord,
          int *stats, int *repiter, int *verbose) {
  specBase = 1;
  init2(&sdid, *nspec, *maxfrac, &nphgram, &totalp,
        &pbuf, &pbase, *maxlocalp, NULL, mycomm, &nbor, pcoord,
        *stats, *repiter, *verbose);
}
void
oh2_init(int **sdid, int nspec, int maxfrac, int **nphgram,
         int **totalp, struct S_particle **pbuf, int **pbase, int maxlocalp,
         void *mycomm, int **nbor, int *pcoord,
         int stats, int repiter, int verbose) {
  specBase = 0;
  init2(sdid, nspec, maxfrac, nphgram, totalp,
        pbuf, pbase, maxlocalp, (struct S_mycommc*)mycomm, NULL, nbor, pcoord,
        stats, repiter, verbose);
}
void
init2(int **sdid, int nspec, int maxfrac, int **nphgram,
      int **totalp, struct S_particle **pbuf, int **pbase, int maxlocalp,
      struct S_mycommc *mycommc, struct S_mycommf *mycommf,
      int **nbor, int *pcoord, int stats, int repiter, int verbose) {
  int ns, nn, nnns, s;

  init1(sdid, nspec, maxfrac, nphgram, totalp, NULL, NULL,
        mycommc, mycommf, nbor, pcoord, stats, repiter, verbose);

  ns = nOfSpecies;  nn = nOfNodes;  nnns = nn * ns;

  nOfLocalPLimit = totalParts = maxlocalp;
  if (*pbuf)
    Particles = *pbuf;
  else
    Particles = *pbuf =
      (struct S_particle*)mem_alloc(sizeof(struct S_particle),
                                    maxlocalp, "Particles");

  MPI_Type_contiguous(sizeof(struct S_particle), MPI_BYTE, &T_Particle);
  MPI_Type_commit(&T_Particle);

  if (!*pbase)  *pbase = (int*)mem_alloc(sizeof(int), 3, "ParticleBase");
  (*pbase)[0] = (*pbase)[1] = (*pbase)[2] = 0;
  secondaryBase = *pbase + 1;  totalLocalParticles = *pbase + 2;

#ifndef OH_POS_AWARE
  SendBuf = (struct S_particle*)mem_alloc(sizeof(struct S_particle), maxlocalp,
                                          "SendBuf");
#endif
  RecvBufBases = (struct S_particle**)mem_alloc(sizeof(struct S_particle*),
                                                2*ns+1, "RecvBufBases");
  SendBufDisps = (int*)mem_alloc(sizeof(int),  nnns, "SendBufDisps");
  RecvBufDisps = (int*)mem_alloc(sizeof(int),  nn, "RecvBufDisps");
  nOfInjections = 0;

  Requests = (MPI_Request*)mem_alloc(sizeof(MPI_Request),
                                     nnns*4+OH_NEIGHBORS*2, "Requests");
  Statuses = (MPI_Status*) mem_alloc(sizeof(MPI_Status),
                                     nnns*4+OH_NEIGHBORS*2, "Statuses");
}
int
oh2_transbound_(int *currmode, int *stats) {
  return(transbound2(*currmode, *stats, 2));
}
int
oh2_transbound(int currmode, int stats) {
  return(transbound2(currmode, stats, 2));
}
int
transbound2(int currmode, int stats, int level) {
  int ret=MODE_NORM_SEC, nn=nOfNodes, ns=nOfSpecies, nnns2=2*nn*ns;
  int i, s, tp;

  stats = stats && statsMode;
  currmode = transbound1(currmode, stats, level);

  if (try_primary2(currmode, level, stats))  ret = MODE_NORM_PRI;
  else if (!Mode_PS(currmode) || !try_stable2(currmode, level, stats)) {
    rebalance2(currmode, level, stats);  ret = MODE_REB_SEC;
  }
  for (i=0; i<nnns2; i++) NOfPLocal[i] = 0;
  for (s=0,tp=0; s<ns*2; s++) {
    TotalP[s] = TotalPNext[s];  tp += TotalPNext[s];
  }
  for (s=0; s<ns*2; s++)  InjectedParticles[s] = 0;
  totalParts = *totalLocalParticles = tp;  nOfInjections = 0;
  return((currMode=ret));
}
static int
try_primary2(int currmode, int level, int stats) {

  if (!try_primary1(currmode, level, stats))  return(FALSE);
  move_to_sendbuf_primary(Mode_PS(currmode), stats);
  exchange_primary_particles(currmode, stats);
  primaryParts = *secondaryBase = TotalPGlobal[myRank];
  return(TRUE);
}
void
exchange_primary_particles(int currmode, int stats) {
  int i, s, nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns, me=myRank;
  int *np, *rnp, *sbd;

  if (stats) oh1_stats_time(STATS_TB_COMM, 0);
  np = NOfPLocal;                       /* &NOfPLocal[0][0][0] */
  rnp = NOfPrimaries;                   /* &NOfPrimaries[0][0][0] */
  sbd = SendBufDisps;                   /* SendBufDisps[0][0] */
  if (currmode==MODE_NORM_PRI) {
    for (s=0; s<ns; s++,np+=nn,rnp+=nn,sbd+=nn) {
                                        /* np=&NOfPLocal[0][s][0] */
                                        /* rnp=&NOfPrimaries[0][s][0] */
                                        /* sbd=&SendBufDisps[s][0] */
      struct S_particle *rb;
      rb = RecvBufBases[s];             /* RecvBufBases[0][s] */
      for (i=0; i<OH_NEIGHBORS; i++) {
        int dst=DstNeighbors[i];
        int src=SrcNeighbors[i];
        int rc;
        MPI_Status st;
        if (dst==me) continue;
        if (src>=0) {
          rc = rnp[src];                /* NOfPrimaries[0][s][src] */
          if (dst>=0)
            MPI_Sendrecv(SendBuf+sbd[dst], np[dst], T_Particle, dst, 0,
                         rb, rc, T_Particle, src, 0, MCW, &st);
          else
            MPI_Recv(rb, rc, T_Particle, src, 0, MCW, &st);
          rb += rc;
        } else if (dst>=0)
          MPI_Send(SendBuf+sbd[dst], np[dst], T_Particle, dst, 0, MCW);
      }
    }
  } else {
    for (s=0; s<ns; s++,np+=nn,rnp+=nn,sbd+=nn) {
                                        /* np=&NOfPLocal[0][s][0] */
                                        /* sbd=&SendBufDisps[s][0] */
                                        /* rnp=&NOfPrimaries[0][s][0] */
      int rdisp=0;
      rnp[me] = 0;                      /* &NOfPrimaries[0][s][me] */
      for (i=0; i<nn; i++) {
        int rc = rnp[i] + rnp[i+nnns];
                        /* NOfPrimaries[0][s][i]+ NofPrimaries[1][s][i] */
        TempArray[i] = rc;
        RecvBufDisps[i] = rdisp;
        rdisp += rc;
        np[i] += np[i+nnns];            /* += NOfPLocal[1][s][i] */
      }
      MPI_Alltoallv(SendBuf, np, sbd, T_Particle,
                    RecvBufBases[s], TempArray, RecvBufDisps, T_Particle, MCW);
    }
  }
}
static int
try_stable2(int currmode, int level, int stats) {

  if (!try_stable1(currmode, level, stats)) return(FALSE);
  exchange_particles(CommList+SLHeadTail[1], SecSLHeadTail[0],
                     Nodes[myRank].parentid, currmode==MODE_NORM_SEC,
                     currmode, stats);
  return(TRUE);
}
static void
rebalance2(int currmode, int level, int stats) {
  int me=myRank, ns=nOfSpecies, s, oldp, newp;

  rebalance1(currmode, level, stats);
  oldp = NodesNext[me].parentid;  newp=Nodes[me].parentid;
  if (nOfInjections && oldp>=0 && oldp!=newp)
    for (s=0; s<ns; s++)  InjectedParticles[ns+s] = 0;
  if (Mode_Is_Norm(currmode))
    exchange_particles(SecRList, SecRLSize,
                       Mode_PS(currmode) ? oldp : -1,
                       1, currmode, stats);
  else
    exchange_particles(SecRList, SecRLSize, -1, 0, currmode, stats);
}
void
move_to_sendbuf_primary(int secondary, int stats) {
  int me=myRank, ns=nOfSpecies, nn=nOfNodes, nnns=nn*ns;
  int s, i, j, *pp;

  if (stats) oh1_stats_time(STATS_TB_MOVE, 0);
  for (s=0,i=me,pp=NOfPrimaries; s<ns; s++,i+=nn,pp+=nn) {
    int t = 0;
    NOfPLocal[i] = 0;                           /* NOfPLocal[0][s][me] */
    for (j=0; j<nn; j++) t += pp[j] + pp[j+nnns];
                            /* NOfPrimaries[0][s][j] + NOfPrimaries[1][s][j] */
    TotalPNext[s] = t;                          /* TotalPNext[0][s] */
    TotalPNext[ns+s] = 0;                       /* TotalPNext[1][s] */
  }
  set_sendbuf_disps(secondary, -1);
  if (nOfInjections)  move_injected_to_sendbuf();

  move_to_sendbuf_uw(0, me, NOfPLocal+me, 0, TotalP, 0, TotalPNext,
                     RecvBufBases);

  if (secondary) move_to_sendbuf_uw(1, -1, NULL, primaryParts, TotalP+ns,
                                    0, TotalPNext+ns, RecvBufBases+ns);

  move_to_sendbuf_dw(0, me, NOfPLocal+me, primaryParts, TotalP,
                     TotalPGlobal[me], TotalPNext);

  set_sendbuf_disps(secondary, -1);
  if (nOfInjections)
    move_injected_from_sendbuf(InjectedParticles, me, RecvBufBases);
}
static void
move_to_sendbuf_secondary(int secondary, int stats) {
  int me=myRank, ns=nOfSpecies, ns2=ns<<1, nn=nOfNodes;
  struct S_node *node = Nodes+me;
  int put[2]={-node->get.prime, -node->get.sec}, pnext[2];
  int sec=node->parentid;
  int nnns=nn*ns;
  int *mynp[2]={NOfPLocal+me,           /* &NOfPLocal[0][0][me] */
                sec<0 ? NULL : NOfPLocal+nnns+sec};
                                        /* &NOfPLocal[1][0][sec] */
  int *mynps;
  int ps, s, i;

  if (stats) oh1_stats_time(STATS_TB_MOVE, 1);
  for (ps=0,i=0; ps<2; ps++) {
    int putme=put[ps], npnext=0;
    mynps=mynp[ps];
    if (mynps==NULL) {
      pnext[ps] = 0;  break;
    }
    if (putme<0) putme = 0;
    for (s=0; s<ns; s++,i++,mynps+=nn) {        /* i=ps*ns+s */
      int stay=*mynps;                          /* NofPLocal[ps][s][me/sec] */
      int tpni=TotalPNext[i];                   /* TotalPNext[ps][s] */
      int inj=InjectedParticles[i];             /* InjectedParticles[ps][s] */
      if (putme<stay) {
        TotalPNext[i] = tpni = tpni + stay - putme;
        stay = putme;
        putme = 0;
      } else
        putme -= stay;
      if (stay>inj) {
        InjectedParticles[ns2+i] = 0;  *mynps = stay - inj;
      } else {
        InjectedParticles[ns2+i] = inj - stay;  *mynps = 0;
      }
      npnext += tpni;
    }
    pnext[ps] = npnext;
  }
  set_sendbuf_disps(secondary, sec);
  if (nOfInjections)  move_injected_to_sendbuf();

  move_to_sendbuf_uw(0, me, mynp[0],            /* &NOfPLocal[0][0][me] */
                     0, TotalP, 0, TotalPNext, RecvBufBases);

  if (secondary) {
    move_to_sendbuf_uw(1, sec, mynp[1],         /* &NOfPLocal[1][0][sec] */
                       primaryParts, TotalP+ns, pnext[0], TotalPNext+ns,
                       RecvBufBases+ns);
    move_to_sendbuf_dw(1, sec, mynp[1], totalParts, TotalP+ns,
                       pnext[0]+pnext[1], TotalPNext+ns);
  } else {
    struct S_particle *rbb=Particles+pnext[0];
    for (s=0; s<ns; s++) {
      RecvBufBases[ns+s] = rbb;                 /* RecvBufBases[1][s] */
      rbb += TotalPNext[ns+s];                  /* TotalPNext[1][s] */
    }
  }
  move_to_sendbuf_dw(0, me, mynp[0], primaryParts, TotalP, pnext[0],
                     TotalPNext);

  set_sendbuf_disps(secondary, sec);
  if (nOfInjections) {
    move_injected_from_sendbuf(InjectedParticles+ns2, me, RecvBufBases);
    if (sec>=0)
      move_injected_from_sendbuf(InjectedParticles+ns2+ns, sec,
                                 RecvBufBases+ns);
  }
  primaryParts = *secondaryBase = pnext[0];
}
void
set_sendbuf_disps(int secondary, int parent) {
  int nn=nOfNodes, ns=nOfSpecies, me=myRank;
  int i, j, k, s, disp;

  for (s=0,i=0,disp=0; s<ns; s++) {
    for (k=0; k<nn; k++,i++) {
      SendBufDisps[i] = disp;                   /* SendBufDisps[s][k] */
      disp += NOfPLocal[i];                     /* NOfPLocal[0][s][k] */
      if (k==me)  disp += InjectedParticles[s]; /* InjectedParticles[0][s] */
    }
  }
  if (secondary) {
    for (s=0,j=0,disp=0; s<ns; s++) {
      for (k=0; k<nn; k++,i++,j++) {
        SendBufDisps[j] += disp;                /* SendBufDisps[s][k] */
        disp += NOfPLocal[i];                   /* NOfPLocal[1][s][k] */
        if (k==parent)  disp += InjectedParticles[ns+s];
      }                                         /* InjectedParticles[1][s] */
    }
  }
}
void
exchange_particles(struct S_commlist *secrlist, int secrlsize, int oldparent,
                   int neighboring, int currmode, int stats) {
  int me=myRank, nn=nOfNodes, ns=nOfSpecies, nnns=nn*ns;
  int newparent=Nodes[me].parentid;
  int s, i, req;

  move_to_sendbuf_secondary(Mode_PS(currmode), stats);
  if (stats) oh1_stats_time(STATS_TB_COMM, 1);
  if (neighboring) {
    req = 0;
    receive_particles(CommList, SLHeadTail[0], &req);
    if (newparent>=0)
      receive_particles(secrlist, secrlsize, &req);
    if (oldparent!=newparent && oldparent>=0)
      receive_particles(CommList+SLHeadTail[1], SecSLHeadTail[0], &req);
    send_particles(CommList+SLHeadTail[0], SLHeadTail[1]-SLHeadTail[0],
                   newparent, oldparent, &req);
    if (oldparent>=0)
      send_particles(CommList+SLHeadTail[1]+SecSLHeadTail[0],
                     SecSLHeadTail[1]-SecSLHeadTail[0], me, newparent, &req);
    MPI_Waitall(req, Requests, Statuses);
  }
  else {
    int ps;
    int *rcount=NOfRecv;
    int *scount=NOfSend;
    struct S_particle **rbb=RecvBufBases;
    for (ps=0; ps<2; ps++,rbb+=ns) {            /* rbb=&RecvBufBases[p][0] */
      int *sbd0=SendBufDisps, *sbd;
      for (s=0; s<ns; s++,rcount+=nn,scount+=nn,sbd0+=nn) {
                                        /* rcount=&NOfRecv[ps][s][0] */
                                        /* sbd0=&SendBufDisps[s][0] */
        int rdisp=0;
        for (i=0; i<nn; i++) {
          RecvBufDisps[i] = rdisp;  rdisp += rcount[i];
        }
        if (ps==0) sbd = sbd0;                  /* &SendBufDisps[s][0] */
        else {
          sbd = TempArray;
          for (i=0; i<nn; i++) {
            int r=Nodes[i].parentid;
            if (r>=0) {
              sbd[i] = sbd0[r];
              sbd0[r] += scount[i];
            }
            else sbd[i] = 0;            /* not necessary becasuse scount[i]=0
                                           but ... */
          }
        }
        MPI_Alltoallv(SendBuf, scount, sbd, T_Particle,
                      rbb[s], rcount, RecvBufDisps, T_Particle, MCW);
        if (ps==0)
          for (i=0; i<nn; i++) sbd0[i] += scount[i];
      }
    }
  }
}
static void
move_to_sendbuf_uw(int ps, int me, int *putmes, int cbase, int *ctp,
                   int nbase, int *ntp, struct S_particle **rbb) {
  int i, in, j, jn, k, s;
  int ns=nOfSpecies, nn=nOfNodes, *sbd=SendBufDisps;
  Decl_Grid_Info();

  for (s=0,i=cbase,j=nbase,k=0; s<ns; s++,i=in,j=jn,sbd+=nn,k+=nn) {
    int putme = putmes ? putmes[k] : 0; /* NOfPLocal[0/1][s][me/sec] */
    in = i + ctp[s];  jn = j + ntp[s];
    if (j<=i) {                         /* upward move only */
      for (; putme>0; i++) {            /* throw my particles to send buf */
        int dst=Subdomain_Id(Particles[i].nid, ps);
        if (dst<0) continue;
        SendBuf[sbd[dst]++] = Particles[i];
        if (dst==me) putme--;
      }
      for (; i<in; i++) {               /* move upward */
        int dst=Subdomain_Id(Particles[i].nid, ps);
        if (dst<0) continue;
        if (dst==me) Particles[j++] = Particles[i];
        else         SendBuf[sbd[dst]++] = Particles[i];
      }
      rbb[s] = Particles + j;           /* receive to bottom */
    } else if (jn>in) {                 /* downward only and thus skip */
      rbb[s] = Particles + j;           /* receive to top */
    } else {                            /* downward and upward */
      int ib, im, jm;
      for (; putme>0; i++) {            /* throw my particles to send buf */
        int dst=Subdomain_Id(Particles[i].nid, ps);
        if (dst<0) continue;
        SendBuf[sbd[dst]++] = Particles[i];
        if (dst==me) putme--;
      }
      ib = i;
      for (; i<j; i++) {                 /* skip downward movers if any */
        int dst=Subdomain_Id(Particles[i].nid, ps);
        if (dst==me && dst>=0)  j++;
      }
      im = i-1; jm = j-1;
      for (; i<in; i++) {               /* move remainders upward */
        int dst=Subdomain_Id(Particles[i].nid, ps);
        if (dst<0) continue;
        if (dst==me) Particles[j++] = Particles[i];
        else         SendBuf[sbd[dst]++] = Particles[i];
      }
      rbb[s] = Particles + j;           /* receive to bottom */
      for (i=im,j=jm; i>=ib; i--) {     /* move first half downward if any */
        int dst=Subdomain_Id(Particles[i].nid, ps);
        if (dst<0) continue;
        if (dst==me) Particles[j--] = Particles[i];
        else         SendBuf[sbd[dst]++] = Particles[i];
      }
    }
  }
}
static void
move_to_sendbuf_dw(int ps, int me, int *putmes, int ctail, int *ctp, int ntail,
                   int *ntp) {
  int i, in, j, jn, k, s, ns=nOfSpecies, nn=nOfNodes, nnnsm1=nn*(ns-1);
  int *sbd=SendBufDisps+nnnsm1;
  Decl_Grid_Info();

  in = ctail;  jn = ntail;
  for (s=ns-1,i=in-1,j=jn-1,k=nnnsm1; s>=0; s--,i=in-1,j=jn-1,sbd-=nn,k-=nn) {
    int putme = putmes ? putmes[k] : 0; /* NOfPLocal[0/1][s][me/sec] */
    in -= ctp[s];  jn -= ntp[s];
    if (i>=j || in>=jn) continue;       /* not downward only and thus skip */
    for (; putme>0; i--) {              /* throw my particles to send buf */
      int dst=Subdomain_Id(Particles[i].nid, ps);
      if (dst<0) continue;
      SendBuf[sbd[dst]++] = Particles[i];
      if (dst==me) putme--;
    }
    for (; i>=in; i--) {                /* move downward */
      int dst=Subdomain_Id(Particles[i].nid, ps);
      if (dst<0) continue;
      if (dst==me) Particles[j--] = Particles[i];
      else         SendBuf[sbd[dst]++] = Particles[i];
    }
  }
}
static void
move_injected_to_sendbuf() {
  struct S_particle *pbuf=Particles+totalParts;
  int ninj=nOfInjections, nn=nOfNodes, sb=specBase;
  int i;
  Decl_Grid_Info();

  for (i=0; i<ninj; i++) {
    int dst = Subdomain_Id(pbuf[i].nid, 0);
    int s = Particle_Spec(pbuf[i].spec-sb);
    if (dst<0) continue;
#ifdef OH_POS_AWARE
    if (dst>=nn)  Primarize_Id(pbuf+i, dst);
#endif
    SendBuf[SendBufDisps[dst+s*nn]++] = pbuf[i];
  }
}
static void
move_injected_from_sendbuf(int *injected, int mysd, struct S_particle **rbb) {
  int nn=nOfNodes, ns=nOfSpecies;
  int *sdisp=SendBufDisps+mysd;
  int s, i;

  for (s=0; s<ns; s++,sdisp+=nn) {
    struct S_particle *rbuf=rbb[s];
    struct S_particle *sbuf=SendBuf+*sdisp;
    int inj=injected[s];
    for (i=0; i<inj; i++)  rbuf[i] = sbuf[i];
    rbb[s] += inj;  *sdisp += inj;
  }
}
static void
receive_particles(struct S_commlist *rlist, int rlsize, int *req) {
  int me=myRank, i, r=*req, nn=nOfNodes, ns=nOfSpecies, sdisp;
  struct S_particle *rbuf;

  for (i=0; i<rlsize; i++) {
    if (rlist[i].rid==me) {
      int count=rlist[i].count, tag=rlist[i].tag;
      rbuf = RecvBufBases[tag];  RecvBufBases[tag] = rbuf + count;
      MPI_Irecv(rbuf, count, T_Particle, rlist[i].sid, tag, MCW,
                Requests+(r++));
    }
    if (rlist[i].sid==me) {
      int count=rlist[i].count, tag=rlist[i].tag, region=rlist[i].region;
      region += nn * (tag<ns ? tag : tag-ns);
      sdisp = SendBufDisps[region];  SendBufDisps[region] = sdisp + count;
                                                /* SendBufDisps[s][region] */
      MPI_Isend(SendBuf+sdisp, count, T_Particle, rlist[i].rid, tag, MCW,
                Requests+(r++));
    }
  }
  *req = r;
}
static void
send_particles(struct S_commlist *slist, int slsize, int myregion,
               int parentregion, int *req) {
  int me=myRank, i, r=*req, nn=nOfNodes, ns=nOfSpecies, sdisp, region;

  for (i=0; i<slsize; i++) {
    if (slist[i].sid==me && (region=slist[i].region)!=myregion &&
                            region != parentregion) {
      int count=slist[i].count, tag=slist[i].tag;
      region += nn * (tag<ns ? tag : tag-ns);
      sdisp = SendBufDisps[region];  SendBufDisps[region] = sdisp + count;
                                                /* SendBufDisps[s][region] */
      MPI_Isend(SendBuf+sdisp, count, T_Particle, slist[i].rid, tag, MCW,
                Requests+(r++));
    }
  }
  *req = r;
}
void
oh2_inject_particle_(struct S_particle *part) {
  oh2_inject_particle(part);
}
void
oh2_inject_particle(struct S_particle *part) {
  const int ns=nOfSpecies, nn=nOfNodes;
  int inj = totalParts + nOfInjections++;
  int s = Particle_Spec(part->spec - specBase);
  int n=part->nid;

#ifndef OH_HAS_SPEC
  if (ns!=1)
    local_errstop("particles cannot be injected when S_particle does not "
                  "have 'spec' element and you have two or more species");
#endif
  if (inj>=nOfLocalPLimit)
    local_errstop("injection causes local particle buffer overflow");
  Particles[inj] = *part;
  if (n<0)  return;
  if (n==RegionId[1]) {
    NOfPLocal[(ns+s)*nn+n]++;
    InjectedParticles[ns+s]++;
  } else {
    NOfPLocal[nn*s+n]++;
    if (n==myRank)  InjectedParticles[s]++;
  }
}
void
oh2_remap_injected_particle_(struct S_particle *part) {
  oh2_remap_injected_particle(part);
}
void
oh2_remap_injected_particle(struct S_particle *part) {
  const int pidx = part - Particles, ns=nOfSpecies, nn=nOfNodes;
  int s, n;

  if (pidx<totalParts || pidx>=totalParts+nOfInjections)
    local_errstop("'part' argument pointing %c%d%c of the particle buffer is "\
                  "not for injected particles",
                  specBase?'(':'[', pidx+specBase, specBase?')':']');
#ifndef OH_HAS_SPEC
  if (ns!=1)
    local_errstop("particles cannot be injected when S_particle does not "
                  "have 'spec' element and you have two or more species");
#endif
  s = Particle_Spec(part->spec - specBase);
  n = part->nid;
  if (n<0)  return;
  if (n==RegionId[1]) {
    NOfPLocal[(ns+s)*nn+n]++;
    InjectedParticles[ns+s]++;
  } else {
    NOfPLocal[nn*s+n]++;
    if (n==myRank)  InjectedParticles[s]++;
  }
}
void
oh2_remove_injected_particle_(struct S_particle *part) {
  oh2_remove_injected_particle(part);
}
void
oh2_remove_injected_particle(struct S_particle *part) {
  const int pidx = part - Particles, ns=nOfSpecies, nn=nOfNodes;
  int s, n;

  if (pidx<totalParts || pidx>=totalParts+nOfInjections)
    local_errstop("'part' argument pointing %c%d%c of the particle buffer is "\
                  "not for injected particles",
                  specBase?'(':'[', pidx+specBase, specBase?')':']');
#ifndef OH_HAS_SPEC
  if (ns!=1)
    local_errstop("particles cannot be injected when S_particle does not "
                  "have 'spec' element and you have two or more species");
#endif
  s = Particle_Spec(part->spec - specBase);
  n = part->nid;
  if (n<0)  return;
  if (n==RegionId[1]) {
    NOfPLocal[(ns+s)*nn+n]--;
    InjectedParticles[ns+s]--;
  } else {
    NOfPLocal[nn*s+n]--;
    if (n==myRank)  InjectedParticles[s]--;
  }
  part->nid = -1;
}
void
oh2_set_total_particles_() {
  set_total_particles();
}
void
oh2_set_total_particles() {
  set_total_particles();
}
int
oh2_max_local_particles_(dint *npmax, int *maxfrac, int *minmargin) {
  return(oh2_max_local_particles(*npmax, *maxfrac, *minmargin));
}
int
oh2_max_local_particles(dint npmax, int maxfrac, int minmargin) {
  int nn, nplint;
  dint npl, npmargin;

  MPI_Comm_size(MCW, &nn);
  if (npmax<=0) errstop("max # of particles should be greater than 0");
  if (maxfrac<=0 || maxfrac>100)
    errstop("load imbalance factor (%d) should be in the range [1..100]",
            maxfrac);
  npl = (npmax-1)/nn + 1; /* ceil(npmax/nn) */
  npmargin = (npl*maxfrac-1)/100 + 1;
  npl += (npmargin<minmargin) ? minmargin : npmargin;
  if (npl>INT_MAX) mem_alloc_error("Particles", 0);
  nplint = npl;
  return(nplint);
}
