# Makefile to generate proposed pyf files from FORTRAN modules
# ------------------------------------------------------------
# J.-M. Brankart, January 2019
# 
include ../macro/make.ada_ifort
#
OBJ = anaobs.so anaqua.so anatra.so anautil.so covariance.so ensaugm.so interp.so meanstd.so obserror.so schurprod.so score_crps.so score_entropy.so score_optimality.so score_rcrv.so spharea.so sphylm.so stoanam.so stogprod.so storfg.so storng.so stotge.so stoutil.so
#
all: $(OBJ)
	@echo "Building all shared libraries: $(OBJ)"
	-mv $(OBJ) lib
#
clean:
	-rm -f $(OBJ)
	-rm -fr build
	-mkdir build
#
anaobs.so: signatures/anaobs.pyf ../src/EnsAnam/anaobs.f90 anaqua.so anatra.so obserror.so storng.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/anaobs.pyf ../src/EnsAnam/anaobs.f90
anaqua.so: signatures/anaqua.pyf ../src/EnsAnam/anaqua.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/anaqua.pyf ../src/EnsAnam/anaqua.f90
anatra.so: signatures/anatra.pyf ../src/EnsAnam/anatra.f90 storng.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/anatra.pyf ../src/EnsAnam/anatra.f90
anautil.so: signatures/anautil.pyf ../src/EnsAnam/anautil.f90 stoutil.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/anautil.pyf ../src/EnsAnam/anautil.f90
covariance.so: signatures/covariance.pyf ../src/EnsStat/covariance.F90 meanstd.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/covariance.pyf ../src/EnsStat/covariance.F90
ensaugm.so: signatures/ensaugm.pyf ../src/EnsAugm/ensaugm.F90 schurprod.so storng.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/ensaugm.pyf ../src/EnsAugm/ensaugm.F90
interp.so: signatures/interp.pyf ../src/InterpTools/interp.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/interp.pyf ../src/InterpTools/interp.f90
meanstd.so: signatures/meanstd.pyf ../src/EnsStat/meanstd.F90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/meanstd.pyf ../src/EnsStat/meanstd.F90
obserror.so: signatures/obserror.pyf ../src/ObsError/obserror.F90 storng.so stoutil.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/obserror.pyf ../src/ObsError/obserror.F90
schurprod.so: signatures/schurprod.pyf ../src/EnsAugm/schurprod.f90 stoanam.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/schurprod.pyf ../src/EnsAugm/schurprod.f90
score_crps.so: signatures/score_crps.pyf ../src/EnsScores/score_crps.F90 anaqua.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/score_crps.pyf ../src/EnsScores/score_crps.F90
score_entropy.so: signatures/score_entropy.pyf ../src/EnsScores/score_entropy.F90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/score_entropy.pyf ../src/EnsScores/score_entropy.F90
score_optimality.so: signatures/score_optimality.pyf ../src/EnsScores/score_optimality.F90 meanstd.so stoutil.so
	-cp wrappers/score_optimality-f2pywrappers2.f90 build/src.linux-x86_64-2.7/signatures
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/score_optimality.pyf ../src/EnsScores/score_optimality.F90
score_rcrv.so: signatures/score_rcrv.pyf ../src/EnsScores/score_rcrv.F90 anaqua.so anatra.so meanstd.so stoutil.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/score_rcrv.pyf ../src/EnsScores/score_rcrv.F90
spharea.so: signatures/spharea.pyf ../src/TranSpHO/spharea.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/spharea.pyf ../src/TranSpHO/spharea.f90
sphylm.so: signatures/sphylm.pyf ../src/TranSpHO/sphylm.F90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/sphylm.pyf ../src/TranSpHO/sphylm.F90
stoanam.so: signatures/stoanam.pyf ../src/StochTools/stoanam.f90 stogprod.so stoutil.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/stoanam.pyf ../src/StochTools/stoanam.f90
stogprod.so: signatures/stogprod.pyf ../src/StochTools/stogprod.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/stogprod.pyf ../src/StochTools/stogprod.f90
storfg.so: signatures/storfg.pyf ../src/StochTools/storfg.F90 storng.so sphylm.so
	-cp wrappers/storfg-f2pywrappers2.f90 build/src.linux-x86_64-2.7/signatures
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/storfg.pyf ../src/StochTools/storfg.F90
storng.so: signatures/storng.pyf ../src/StochTools/storng.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/storng.pyf ../src/StochTools/storng.f90
stotge.so: signatures/stotge.pyf ../src/StochTools/stotge.f90 storng.so
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/stotge.pyf ../src/StochTools/stotge.f90
stoutil.so: signatures/stoutil.pyf ../src/StochTools/stoutil.f90
	$(F2PY) $(F2PYFLAGS) --build-dir build --f90flags="$(FPPFLAGS) $(FFLAGS)" -c signatures/stoutil.pyf ../src/StochTools/stoutil.f90
