#
#FC	 = ifc
#FC	= ifort
FC      = gfortran
CC	= gcc-9
#CC	 = icc
LD	= $(FC)
F77	= $(FC)
SWP	= swplist
RM	= /bin/rm -f
MP	=
ABI	= 
ISA	= 
PROC	= 
ARCH	= $(MP) $(ABI) $(ISA)
OLEVEL  =  -O5
#OLEVEL  =  -O3 -scalar_rep  -vec -prefetch -tpp6
#OLEVEL  =  -O3 -scalar_rep  -vec -prefetch -tpp6  -axK -pad
#FOPTS   = -r8  
FOPTS   = -fdefault-real-8

COPTS	= 
FFLAGS	= $(ARCH) $(OLEVEL) $(FOPTS)
CFLAGS	= $(ARCH) $(COPTS)
#LIBS	 = -lPEPCF90
LIBS	= 
LDFLAGS	= $(ARCH) $(OLEVEL)
PROF	=

FOBJS   = dims.o dimsc.o\
          flo_var.o mesh_var.o in_var.o\
          geo_var.o solv_var.o psm_var.o mg_var.o out_var.o\
          geo_param.o flo_param.o solv_param.o mg_param.o\
          alloc_flo.o alloc_mesh.o alloc_geo.o\
          alloc_solv.o alloc_psm.o alloc_mg.o alloc_out.o\
          addw.o bcfar.o bcwall.o bound.o\
          collc.o coord.o cutime.o\
          delt.o dflux.o dfluxc.o eflux.o euler.o\
          ffmax.o forcf.o geom.o gmesh.o\
          graph.o grid.o zoomgrid.o halo.o idmax.o idmin.o\
          init.o input.o intpl.o mesh.o metric.o\
          monitr.o movco.o movfin.o nsflux.o\
          output.o plot.o prntff.o prntx.o psmoo.o\
          rplot.o setup.o shift.o splif.o step.o\
          tanhds.o turbbl.o viscf.o vmesh.o xpand.o flo103.o\
          turb2.o ismax.o ismin.o cvmgp.o\
          psgs.o lgraph.o pscont.o 

COBJS	= 
OBJS	= $(FOBJS) $(COBJS)
EXEC	= flo103

$(EXEC):	$(OBJS)
	$(LD) -o $@ $(PROF) $(LDFLAGS) $(QV_LOPT) $(OBJS) $(LIBS)

clean:
	$(RM) $(EXEC) $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .F .c .f .swp

.F.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.f.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.c.o:
	$(CC)  -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<

.F.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.f.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.c.swp:
	$(SWP) -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<
