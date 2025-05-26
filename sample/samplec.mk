#  File: samplec.mk
#  Version 1.1.1 (2015/10/23)
#  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
#                           (ACCMS, Kyoto University)
#  This program can be freely used, redistributed and modified for non-
#  commercial purpose providing that the copyright notice above remains
#  unchanged.
LINKER          = $(CC)

COMMONHDRS      = oh_config.h oh_part.h oh_stats.h
OHHDRS          = ohhelp1.h ohhelp2.h ohhelp3.h ohhelp4p.h
OHOBJS          = ohhelp1.o ohhelp2.o ohhelp3.o ohhelp4p.o

CHDRS           = ohhelp_c.h
COBJS           = simulator.o sample.o

OBJS            = $(COBJS) $(OHOBJS)

simulator:      $(OBJS)
                $(LINKER) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $@

$(COBJS):%.o:   %.c $(COMMONHDRS) $(CHDRS)
                $(CC) $(CFLAGS) -c $< -o $@
$(OHOBJS):%.o:  %.c $(COMMONHDRS) $(OHHDRS)
                $(CC) $(CFLAGS) -c $< -o $@

clean:;
                rm *.o
