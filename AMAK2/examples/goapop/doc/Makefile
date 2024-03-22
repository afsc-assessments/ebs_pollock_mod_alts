ifdef ComSpec
    RM=del /F /Q
else
    RM=rm -rf
endif

all: pdf 
# app clean

R_OPTS=--no-save --no-restore --no-init-file --no-site-file

pdf: goapop_amak.Rmd 
#	cp smbkc.R smbkc.Rmd
	R ${R_OPTS} -e 'library(knitr); knit("goapop_amak.Rmd") '
#knit("app.Rmd")
	pandoc -V geometry:margin=1in goapop_amak.md -o goapop_amak.pdf 
	open goapop_amak.pdf

word: goapop_amak.Rmd 
#	cp smbkc.R smbkc.Rmd
	# R ${R_OPTS} -e 'library(knitr); knit("smbkc.Rmd"); knit("app.Rmd")'
	pandoc -V geometry:margin=1in goapop_amak.md -o goapop_amak.docx --bibliography=references/Gmacs.bib
	open goapop_amak.docx


clean:
	$(RM) goapop_amak.md 
