/* File: oh_config.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#ifndef OH_DIMENSION
#define OH_DIMENSION    3
#endif

/* If you want to activate level-4p functions, remove this comment surrounding
   the line below.
#define OH_LIB_LEVEL_4P
*/
/* If you want to activate level-4s functions, remove this comment surrounding
   the line below.
#define OH_LIB_LEVEL_4S
*/
#ifdef  OH_LIB_LEVEL_4P
#define OH_LIB_LEVEL_4PS
#endif
#ifdef  OH_LIB_LEVEL_4S
#define OH_LIB_LEVEL_4PS
#endif
#ifdef  OH_LIB_LEVEL_4PS
#define OH_LIB_LEVEL 4
/* If you want to use level-4p/4s functions with a large simulation space,
   remove this comment surrounding the line below.
#define OH_BIG_SPACE
*/
/* If you want to let level-4p/4s's particle mapping functions run without
   checking the consistency of their arguments, remove this comment
   surrounding the line below.
#define OH_NO_CHECK
*/
#else
#ifndef OH_LIB_LEVEL
#define OH_LIB_LEVEL 3
#endif
#endif
