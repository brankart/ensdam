# Makefile to generate proposed pyf files from FORTRAN modules
# ------------------------------------------------------------
# J.-M. Brankart, January 2019
# 
include Makefile.macro
#
SRCDIR = ../src
BLDDIR = ../build
#
#OBJ = anaobs.so anaqua.so anatra.so anautil.so covariance.so ensaugm.so interp.so meanstd.so obserror.so schurprod.so score_crps.so score_entropy.so score_optimality.so score_rcrv.so spharea.so sphylm.so stoanam.so stogprod.so storfg.so storng.so stotge.so stoutil.so
OBJ = anaobs.so anaqua.so anatra.so anautil.so covariance.so ensaugm.so interp.so meanstd.so obserror.so schurprod.so score_crps.so score_entropy.so score_rcrv.so spharea.so sphylm.so stoanam.so stogprod.so storng.so stotge.so stoutil.so
#
all: $(OBJ)
	@echo "Building all shared libraries: $(OBJ)"
	-mv *.so modules
#
clean:
	-rm -f $(OBJ)
	-rm -fr build
	-mkdir build
#
anaobs.so: $(SIGDIR)/anaobs.pyf $(SRCDIR)/EnsAnam/anaobs.f90 $(BLDDIR)/anaqua.o $(BLDDIR)/anatra.o $(BLDDIR)/obserror.o $(BLDDIR)/storng.o $(BLDDIR)/stoutil.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/anaobs.pyf $(SRCDIR)/EnsAnam/anaobs.f90 $(BLDDIR)/anaqua.o $(BLDDIR)/anatra.o $(BLDDIR)/obserror.o $(BLDDIR)/storng.o $(BLDDIR)/stoutil.o
anaqua.so: $(SIGDIR)/anaqua.pyf $(SRCDIR)/EnsAnam/anaqua.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/anaqua.pyf $(SRCDIR)/EnsAnam/anaqua.f90
anatra.so: $(SIGDIR)/anatra.pyf $(SRCDIR)/EnsAnam/anatra.f90 $(BLDDIR)/storng.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/anatra.pyf $(SRCDIR)/EnsAnam/anatra.f90 $(BLDDIR)/storng.o
anautil.so: $(SIGDIR)/anautil.pyf $(SRCDIR)/EnsAnam/anautil.f90 $(BLDDIR)/stoutil.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/anautil.pyf $(SRCDIR)/EnsAnam/anautil.f90 $(BLDDIR)/stoutil.o
covariance.so: $(SIGDIR)/covariance.pyf $(SRCDIR)/EnsStat/covariance.F90 $(BLDDIR)/meanstd.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/covariance.pyf $(SRCDIR)/EnsStat/covariance.F90 $(BLDDIR)/meanstd.o
ensaugm.so: $(SIGDIR)/ensaugm.pyf $(SRCDIR)/EnsAugm/ensaugm.F90 $(BLDDIR)/schurprod.o $(BLDDIR)/storng.o $(BLDDIR)/stoanam.o $(BLDDIR)/stoutil.o $(BLDDIR)/stogprod.o $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/ensaugm.pyf $(SRCDIR)/EnsAugm/ensaugm.F90 $(BLDDIR)/schurprod.o $(BLDDIR)/storng.o $(BLDDIR)/stoanam.o $(BLDDIR)/stoutil.o $(BLDDIR)/stogprod.o $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
interp.so: $(SIGDIR)/interp.pyf $(SRCDIR)/InterpTools/interp.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/interp.pyf $(SRCDIR)/InterpTools/interp.f90
meanstd.so: $(SIGDIR)/meanstd.pyf $(SRCDIR)/EnsStat/meanstd.F90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/meanstd.pyf $(SRCDIR)/EnsStat/meanstd.F90
obserror.so: $(SIGDIR)/obserror.pyf $(SRCDIR)/ObsError/obserror.F90 $(BLDDIR)/storng.o $(BLDDIR)/stoutil.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/obserror.pyf $(SRCDIR)/ObsError/obserror.F90 $(BLDDIR)/storng.o $(BLDDIR)/stoutil.o
schurprod.so: $(SIGDIR)/schurprod.pyf $(SRCDIR)/EnsAugm/schurprod.f90 $(BLDDIR)/stoanam.o $(BLDDIR)/stoutil.o $(BLDDIR)/stogprod.o $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/schurprod.pyf $(SRCDIR)/EnsAugm/schurprod.f90 $(BLDDIR)/stoanam.o $(BLDDIR)/stoutil.o $(BLDDIR)/stogprod.o $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
score_crps.so: $(SIGDIR)/score_crps.pyf $(SRCDIR)/EnsScores/score_crps.F90 $(BLDDIR)/anaqua.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/score_crps.pyf $(SRCDIR)/EnsScores/score_crps.F90 $(BLDDIR)/anaqua.o
score_entropy.so: $(SIGDIR)/score_entropy.pyf $(SRCDIR)/EnsScores/score_entropy.F90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/score_entropy.pyf $(SRCDIR)/EnsScores/score_entropy.F90
#score_optimality.so: $(SIGDIR)/score_optimality.pyf $(SRCDIR)/EnsScores/score_optimality.F90 $(BLDDIR)/meanstd.o $(BLDDIR)/stoutil.o
	#-cp $(WRAPDIR)/score_optimality-f2pywrappers2.f90 build/src.linux-x86_64-2.7/$(SIGDIR)
	#$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/score_optimality.pyf $(SRCDIR)/EnsScores/score_optimality.F90 $(BLDDIR)/meanstd.o $(BLDDIR)/stoutil.o
score_rcrv.so: $(SIGDIR)/score_rcrv.pyf $(SRCDIR)/EnsScores/score_rcrv.F90 $(BLDDIR)/anaqua.o $(BLDDIR)/anatra.o $(BLDDIR)/meanstd.o $(BLDDIR)/stoutil.o $(BLDDIR)/storng.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/score_rcrv.pyf $(SRCDIR)/EnsScores/score_rcrv.F90 $(BLDDIR)/anaqua.o $(BLDDIR)/anatra.o $(BLDDIR)/meanstd.o $(BLDDIR)/stoutil.o $(BLDDIR)/storng.o
spharea.so: $(SIGDIR)/spharea.pyf $(SRCDIR)/TranSpHO/spharea.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/spharea.pyf $(SRCDIR)/TranSpHO/spharea.f90
sphylm.so: $(SIGDIR)/sphylm.pyf $(SRCDIR)/TranSpHO/sphylm.F90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/sphylm.pyf $(SRCDIR)/TranSpHO/sphylm.F90
stoanam.so: $(SIGDIR)/stoanam.pyf $(SRCDIR)/StochTools/stoanam.f90 $(BLDDIR)/stogprod.o $(BLDDIR)/stoutil.o $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/stoanam.pyf $(SRCDIR)/StochTools/stoanam.f90 $(BLDDIR)/stogprod.o $(BLDDIR)/stoutil.o $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
stogprod.so: $(SIGDIR)/stogprod.pyf $(SRCDIR)/StochTools/stogprod.f90 $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/stogprod.pyf $(SRCDIR)/StochTools/stogprod.f90 $(BLDDIR)/constant.o $(BLDDIR)/qxgs.o
#storfg.so: $(SIGDIR)/storfg.pyf $(SRCDIR)/StochTools/storfg.F90 $(BLDDIR)/storng.o $(BLDDIR)/sphylm.o
	#-cp $(WRAPDIR)/storfg-f2pywrappers2.f90 build/src.linux-x86_64-2.7/$(SIGDIR)
	#$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/storfg.pyf $(SRCDIR)/StochTools/storfg.F90 $(BLDDIR)/storng.o $(BLDDIR)/sphylm.o
storng.so: $(SIGDIR)/storng.pyf $(SRCDIR)/StochTools/storng.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/storng.pyf $(SRCDIR)/StochTools/storng.f90
stotge.so: $(SIGDIR)/stotge.pyf $(SRCDIR)/StochTools/stotge.f90 $(BLDDIR)/storng.o
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/stotge.pyf $(SRCDIR)/StochTools/stotge.f90 $(BLDDIR)/storng.o
stoutil.so: $(SIGDIR)/stoutil.pyf $(SRCDIR)/StochTools/stoutil.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c $(SIGDIR)/stoutil.pyf $(SRCDIR)/StochTools/stoutil.f90
