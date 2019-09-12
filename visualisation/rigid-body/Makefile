all: tikz-RKMK4.pdf

%.pdf: %.tex
	latexmk -f -pdf -pdflatex="pdflatex -shell-escape -interaction=nonstopmode" -use-make $<

clean:
	latexmk -c
