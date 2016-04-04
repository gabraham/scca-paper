PAPER=paper
SUPP=supplementary
PAPER_PLAIN=paper_plain

full: ${PAPER}.pdf ${SUPP}.pdf

scca_timing-crop.pdf: ../code/scca_timing.pdf
	pdfcrop ../code/scca_timing.pdf scca_timing-crop.pdf

scca_timing-crop.eps: scca_timing-crop.pdf
	pdftops scca_timing-crop.pdf scca_timing-crop.eps

scca_timing_full-crop.pdf: ../code/scca_timing_full.pdf
	pdfcrop ../code/scca_timing_full.pdf scca_timing_full-crop.pdf

${PAPER}.pdf: ${PAPER}.tex ${PAPER}.bib scca_timing-crop.eps
	latex ${PAPER}
	bibtex ${PAPER}
	latex ${PAPER}
	latex ${PAPER}	
	dvips ${PAPER}.dvi
	ps2pdf ${PAPER}.ps

${PAPER_PLAIN}.pdf: ${PAPER_PLAIN}.tex ${PAPER}.bib scca_timing-crop.pdf 
	latex ${PAPER_PLAIN}
	bibtex ${PAPER_PLAIN}
	latex ${PAPER_PLAIN}
	latex ${PAPER_PLAIN}	
	dvips ${PAPER_PLAIN}.dvi
	ps2pdf ${PAPER_PLAIN}.ps

${SUPP}.pdf: ${SUPP}.tex ${PAPER}.bib scca_timing_full-crop.pdf
	pdflatex ${SUPP}
	bibtex ${SUPP}
	pdflatex ${SUPP}
	pdflatex ${SUPP}	

clean:
	rm -f *.out *.aux *.bbl *.blg *.dvi *.log *.toc ${PAPER}.ps ${PAPER}.pdf
