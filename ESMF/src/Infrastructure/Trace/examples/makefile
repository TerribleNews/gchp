# $Id$

ALL:      build_examples
run:      run_examples
run_uni:  run_examples_uni

LOCDIR	  = src/Infrastructure/Trace/examples

.NOTPARALLEL:
EXAMPLES_BUILD    = $(ESMF_EXDIR)/ESMF_TraceEx \
                    $(ESMF_EXDIR)/ESMF_TraceUserEx

EXAMPLES_RUN      = run_ESMF_TraceEx \
                    run_ESMF_TraceUserEx \

EXAMPLES_RUN_UNI  = run_ESMF_TraceEx_uni \
                    run_ESMF_TraceUserEx_uni

include $(ESMF_DIR)/makefile

DIRS	  = 

CLEANDIRS   = 
CLEANFILES  = $(EXAMPLES_BUILD) $(ESMF_EXDIR)/traceout
CLOBBERDIRS = 

#
#  ESMF_Trace
#

# --- ESMF_TraceEx

run_ESMF_TraceEx:
	env ESMF_RUNTIME_TRACE=ON $(MAKE) EXNAME=Trace NP=4 exfrun

run_ESMF_TraceEx_uni:
	env ESMF_RUNTIME_TRACE=ON $(MAKE) EXNAME=Trace NP=1 exfrun

# --- ESMF_TraceUserEx

run_ESMF_TraceUserEx:
	env ESMF_RUNTIME_TRACE=ON $(MAKE) EXNAME=TraceUser NP=4 exfrun

run_ESMF_TraceUserEx_uni:
	env ESMF_RUNTIME_TRACE=ON $(MAKE) EXNAME=TraceUser NP=1 exfrun

