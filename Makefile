qtlX.pdf: qtlX.tex qtlX.bib genetics.bst Makefile Figs/fig_sim.eps Figs/fig_prob305.eps
	pdflatex qtlX
	bibtex qtlX
	pdflatex qtlX
	pdflatex qtlX
	pdflatex qtlX

pdf: qtlX.pdf

all: pdf qtlX_source.tgz

Figs/fig_sim.eps: R/plot_sim.R
	cd R; R CMD BATCH plot_sim.R plot_sim.Rout

Figs/fig_prob305.eps: R/plot_analysis.R
	cd R; R CMD BATCH plot_analysis.R plot_analysis.Rout

clean:
	\rm -f *.aux *.bbl *.blg *.log *.bak *~ *.Rout */*~ */*.Rout */*.aux */*.log

cleanall:
	\rm -f *.aux *.bbl *.blg *.log *.pdf *.bak *~ *.Rout */*~ */*.Rout Figs/*.eps */*.aux */*.log

qtlX_source.tgz: qtlX.tex genetics.bst qtlX.bib Figs/fig1.eps Figs/fig2.eps Figs/fig3.eps Figs/fig4.eps
	tar czvf qtlX_source.tgz qtlX.tex genetics.bst phyloqtl.bib Figs/fig?.eps
