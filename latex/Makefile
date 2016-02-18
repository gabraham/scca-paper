PAPER=paper
SUPP=supplementary

full: ${PAPER}.pdf ${SUPP}.pdf

scca_timing-crop.pdf: ../code/scca_timing.pdf
	pdfcrop ../code/scca_timing.pdf scca_timing-crop.pdf

scca_timing_full-crop.pdf: ../code/scca_timing_full.pdf
	pdfcrop ../code/scca_timing_full.pdf scca_timing_full-crop.pdf

${PAPER}.pdf: ${PAPER}.tex ${PAPER}.bib scca_timing-crop.pdf 
	pdflatex ${PAPER}
	bibtex ${PAPER}
	pdflatex ${PAPER}
	pdflatex ${PAPER}	

${SUPP}.pdf: ${SUPP}.tex ${PAPER}.bib scca_timing_full-crop.pdf
	pdflatex ${SUPP}
	bibtex ${SUPP}
	pdflatex ${SUPP}
	pdflatex ${SUPP}	

clean:
	rm -f *.out *.aux *.bbl *.blg *.dvi *.log *.toc ${PAPER}.ps ${PAPER}.pdf
