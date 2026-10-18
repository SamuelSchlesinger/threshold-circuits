.PHONY: all check clean

PAPER = threshold_circuit_representations

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(PAPER).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(PAPER).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(PAPER).tex

check:
	python3 weighted_pivot_validation.py
	python3 interval_entropy_validation.py
	python3 zonotope_validation.py
	python3 scalar_coding_validation.py

clean:
	rm -f $(PAPER).aux $(PAPER).log $(PAPER).out $(PAPER).toc
