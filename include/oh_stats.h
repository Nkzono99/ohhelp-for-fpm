/* File: oh_stats.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#define STATS_TRANSBOUND  0
#define STATS_TRY_STABLE  (STATS_TRANSBOUND + 1)
#define STATS_REBALANCE   (STATS_TRY_STABLE + 1)
#define STATS_REB_COMM    (STATS_REBALANCE + 1)
#define STATS_TB_MOVE     (STATS_REB_COMM + 1)
#ifdef  OH_LIB_LEVEL_4PS
#define STATS_TB_SORT     (STATS_TB_MOVE + 1)
#define STATS_TB_COMM     (STATS_TB_SORT + 1)
#else
#define STATS_TB_COMM     (STATS_TB_MOVE + 1)
#endif
#define STATS_TIMINGS     (STATS_TB_COMM + 1)

#ifdef OH_DEFINE_STATS
static char *StatsTimeStrings[2*STATS_TIMINGS] = {
  "transbound",         "",
  "try_stable",         "",
  "rebalance",          "",
  "reb comm create",    "",
  "part move[pri]",     "part move[sec]",
#ifdef  OH_LIB_LEVEL_4PS
  "part sort[pri]",     "part sort[sec]",
#endif
  "part comm[pri]",     "part comm[sec]",
};
#endif
