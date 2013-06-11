all:
	latex -shell-escape ester.tex
	bibtex ester.aux
	latex -shell-escape ester.tex
	latex -shell-escape ester.tex # To set the correct page numbers in ToC
	dvipdfm -o ../ester.pdf ester.dvi
fast:
	latex -shell-escape ester.tex
	dvipdfm -o ../ester.pdf ester.dvi
clean:
	rm -f *.aux
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.fcn
	rm -f *.dvi
