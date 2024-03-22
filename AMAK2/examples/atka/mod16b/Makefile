EXEC = amak
PROJ = main
ifeq ($(SAFE),TRUE)
	DIST = ~/_mymods/afsc-assessments/AMAK/src/
else
	DIST = ~/_mymods/afsc-assessments/AMAK/src/
endif
ARGS = -nox -iprint 50 

ifdef ComSpec
    RM=del /F /Q
else
    RM=rm -rf
endif

all: mpd
$(EXEC): $(DIST)$(EXEC).tpl
	ln -sf $(DIST)$(EXEC) $@
	ln -sf $(DIST)$(PROJ) $(PROJ)

$(DIST)$(EXEC).tpl:
	$(MAKE) --directory=~/_mymods/AMAK/src

mpd: $(EXEC)
	./$(EXEC) $(ARGS)

mcmc: $(EXEC)
	./$(EXEC) $(ARGS) -mcmc 1000000 -mcsave 200
	./$(EXEC) $(ARGS) -mceval

retro: 
	./runretro.sh    

proj: $(PROJ)
	./$(PROJ) 

debug: $(EXEC)
	./$(EXEC) $(ARGS)

plot:
	R CMD BATCH plot.R

clean:
	@$(RM) $(EXEC)
	@$(RM) $(EXEC)  $(EXEC).[brces]* $(EXEC).*[0123456789] *.rpt *.log variance gradient.* *tmp
	@$(RM) admodel.*
	@$(RM) checkfile.rep
	@$(RM) mcout.rep
	@$(RM) plot.Rout
	@$(RM) Rplots.pdf
	@$(RM) For_R.rep
	@$(RM) Fprof.yld
	@$(RM) amak.prj
	@$(RM) amak.par
	@$(RM) SIS_out.rep
	@$(RM) mceval.dat

