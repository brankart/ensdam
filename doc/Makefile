# Makefile to compile documentation
# ---------------------------------
# J.-M. Brankart, November 2018
# 
# Name of latex compiler
LATEX = pdflatex

# Define source and executables
src = $(wildcard *.tex)
pdf = $(src:.tex=.pdf)

# Compilation rules
%.pdf : %.tex
	$(LATEX) $<
	$(LATEX) $<

doc: $(pdf) clean

clean:
	rm -f *.aux *.log
