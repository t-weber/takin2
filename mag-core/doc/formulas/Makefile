#
# makefile for takin formula notes
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2021
# @license see 'LICENSE' file
#



# -----------------------------------------------------------------------------
# tools
# -----------------------------------------------------------------------------
PDFLATEX = pdflatex
BIBTEX = bibtex
CONVERT = convert
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# meta rules
# -----------------------------------------------------------------------------
.PHONY: all clean figures_svg2pdf

all: figures_svg2pdf formulas.pdf

clean:
	find . \( -name "*.aux" -o -name "*.out" \
		-o -name "*.bbl" -o -name "*.blg" \
		-o -name "*.toc" -o -name "*.log" \) \
			-exec rm -fv {} \;
	rm -fv formulas.pdf

figures_svg2pdf: $(patsubst %.svg, %.pdf, $(wildcard *.svg))
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# general rules
# -----------------------------------------------------------------------------
%.pdf: %.tex
	@echo "Building $< -> $@..."
	$(PDFLATEX) $<
	$(BIBTEX) $(patsubst %.tex, %, $<)
	$(PDFLATEX) $<
	$(PDFLATEX) $<

%.pdf: %.svg
	@echo "Converting $< -> $@..."
	$(CONVERT) $< $@
# -----------------------------------------------------------------------------
