# Makefile created by mkmf 19.3.0


MK_TEMPLATE = Makefile.template
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: ../lib/libensdam.a
anaobs.o: ../src/EnsAnam/anaobs.f90 anaqua.o anatra.o obserror.o storng.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsAnam/anaobs.f90
anaqua.o: ../src/EnsAnam/anaqua.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsAnam/anaqua.f90
anatra.o: ../src/EnsAnam/anatra.f90 storng.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsAnam/anatra.f90
anautil.o: ../src/EnsAnam/anautil.f90 stoutil.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsAnam/anautil.f90
constant.o: ../src/StochTools/constant.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/constant.f90
covariance.o: ../src/EnsStat/covariance.F90 meanstd.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsStat/covariance.F90
ensaugm.o: ../src/EnsAugm/ensaugm.F90 schurprod.o storng.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsAugm/ensaugm.F90
interp.o: ../src/InterpTools/interp.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/InterpTools/interp.f90
meanstd.o: ../src/EnsStat/meanstd.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsStat/meanstd.F90
obserror.o: ../src/ObsError/obserror.F90 storng.o stoutil.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/ObsError/obserror.F90
qxgs.o: ../src/StochTools/qxgs.f90 constant.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/qxgs.f90
schurprod.o: ../src/EnsAugm/schurprod.f90 stoanam.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsAugm/schurprod.f90
score_crps.o: ../src/EnsScores/score_crps.F90 anaqua.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsScores/score_crps.F90
score_entropy.o: ../src/EnsScores/score_entropy.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsScores/score_entropy.F90
score_optimality.o: ../src/EnsScores/score_optimality.F90 meanstd.o stoutil.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsScores/score_optimality.F90
score_rcrv.o: ../src/EnsScores/score_rcrv.F90 anaqua.o anatra.o meanstd.o stoutil.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/EnsScores/score_rcrv.F90
spharea.o: ../src/TranSpHO/spharea.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/TranSpHO/spharea.f90
sphylm.o: ../src/TranSpHO/sphylm.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/TranSpHO/sphylm.F90
stoanam.o: ../src/StochTools/stoanam.f90 stogprod.o stoutil.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/stoanam.f90
stogprod.o: ../src/StochTools/stogprod.f90 constant.o qxgs.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/stogprod.f90
storfg.o: ../src/StochTools/storfg.F90 storng.o sphylm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/storfg.F90
storng.o: ../src/StochTools/storng.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/storng.f90
stotge.o: ../src/StochTools/stotge.f90 storng.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/stotge.f90
stoutil.o: ../src/StochTools/stoutil.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/StochTools/stoutil.f90
./covariance.F90: ../src/EnsStat/covariance.F90
	cp ../src/EnsStat/covariance.F90 .
./storfg.F90: ../src/StochTools/storfg.F90
	cp ../src/StochTools/storfg.F90 .
./spharea.f90: ../src/TranSpHO/spharea.f90
	cp ../src/TranSpHO/spharea.f90 .
./stoutil.f90: ../src/StochTools/stoutil.f90
	cp ../src/StochTools/stoutil.f90 .
./sphylm.F90: ../src/TranSpHO/sphylm.F90
	cp ../src/TranSpHO/sphylm.F90 .
./qxgs.f90: ../src/StochTools/qxgs.f90
	cp ../src/StochTools/qxgs.f90 .
./interp.f90: ../src/InterpTools/interp.f90
	cp ../src/InterpTools/interp.f90 .
./schurprod.f90: ../src/EnsAugm/schurprod.f90
	cp ../src/EnsAugm/schurprod.f90 .
./anautil.f90: ../src/EnsAnam/anautil.f90
	cp ../src/EnsAnam/anautil.f90 .
./score_entropy.F90: ../src/EnsScores/score_entropy.F90
	cp ../src/EnsScores/score_entropy.F90 .
./score_optimality.F90: ../src/EnsScores/score_optimality.F90
	cp ../src/EnsScores/score_optimality.F90 .
./meanstd.F90: ../src/EnsStat/meanstd.F90
	cp ../src/EnsStat/meanstd.F90 .
./stoanam.f90: ../src/StochTools/stoanam.f90
	cp ../src/StochTools/stoanam.f90 .
./score_crps.F90: ../src/EnsScores/score_crps.F90
	cp ../src/EnsScores/score_crps.F90 .
./anaobs.f90: ../src/EnsAnam/anaobs.f90
	cp ../src/EnsAnam/anaobs.f90 .
./storng.f90: ../src/StochTools/storng.f90
	cp ../src/StochTools/storng.f90 .
./ensaugm.F90: ../src/EnsAugm/ensaugm.F90
	cp ../src/EnsAugm/ensaugm.F90 .
./stogprod.f90: ../src/StochTools/stogprod.f90
	cp ../src/StochTools/stogprod.f90 .
./obserror.F90: ../src/ObsError/obserror.F90
	cp ../src/ObsError/obserror.F90 .
./score_rcrv.F90: ../src/EnsScores/score_rcrv.F90
	cp ../src/EnsScores/score_rcrv.F90 .
./anatra.f90: ../src/EnsAnam/anatra.f90
	cp ../src/EnsAnam/anatra.f90 .
./stotge.f90: ../src/StochTools/stotge.f90
	cp ../src/StochTools/stotge.f90 .
./constant.f90: ../src/StochTools/constant.f90
	cp ../src/StochTools/constant.f90 .
./anaqua.f90: ../src/EnsAnam/anaqua.f90
	cp ../src/EnsAnam/anaqua.f90 .
SRC = ../src/ObsError/obserror.F90 ../src/StochTools/stogprod.f90 ../src/StochTools/storfg.F90 ../src/StochTools/storng.f90 ../src/EnsStat/meanstd.F90 ../src/StochTools/stotge.f90 ../src/EnsScores/score_crps.F90 ../src/EnsAnam/anautil.f90 ../src/EnsAugm/schurprod.f90 ../src/EnsScores/score_optimality.F90 ../src/StochTools/stoanam.f90 ../src/EnsAnam/anatra.f90 ../src/TranSpHO/spharea.f90 ../src/InterpTools/interp.f90 ../src/EnsScores/score_rcrv.F90 ../src/EnsAugm/ensaugm.F90 ../src/EnsAnam/anaqua.f90 ../src/EnsScores/score_entropy.F90 ../src/StochTools/constant.f90 ../src/StochTools/stoutil.f90 ../src/EnsStat/covariance.F90 ../src/TranSpHO/sphylm.F90 ../src/StochTools/qxgs.f90 ../src/EnsAnam/anaobs.f90
OBJ = obserror.o stogprod.o storfg.o storng.o meanstd.o stotge.o score_crps.o anautil.o schurprod.o score_optimality.o stoanam.o anatra.o spharea.o interp.o score_rcrv.o ensaugm.o anaqua.o score_entropy.o constant.o stoutil.o covariance.o sphylm.o qxgs.o anaobs.o
OFF = ../src/EnsStat/covariance.F90 ../src/StochTools/storfg.F90 ../src/TranSpHO/spharea.f90 ../src/StochTools/stoutil.f90 ../src/TranSpHO/sphylm.F90 ../src/StochTools/qxgs.f90 ../src/InterpTools/interp.f90 ../src/EnsAugm/schurprod.f90 ../src/EnsAnam/anautil.f90 ../src/EnsScores/score_entropy.F90 ../src/EnsScores/score_optimality.F90 ../src/EnsStat/meanstd.F90 ../src/StochTools/stoanam.f90 ../src/EnsScores/score_crps.F90 ../src/EnsAnam/anaobs.f90 ../src/StochTools/storng.f90 ../src/EnsAugm/ensaugm.F90 ../src/StochTools/stogprod.f90 ../src/ObsError/obserror.F90 ../src/EnsScores/score_rcrv.F90 ../src/EnsAnam/anatra.f90 ../src/StochTools/stotge.f90 ../src/StochTools/constant.f90 ../src/EnsAnam/anaqua.f90
clean: neat
	-rm -f .../lib/libensdam.a.cppdefs $(OBJ) ../lib/libensdam.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

../lib/libensdam.a: $(OBJ)
	$(AR) $(ARFLAGS) ../lib/libensdam.a $(OBJ)

