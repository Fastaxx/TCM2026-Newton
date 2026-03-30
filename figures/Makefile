# directory containing your fig1.tex, fig2.tex
FIGS = fig1 fig2

.PHONY: all clean
all: $(FIGS:=.png)

%.pdf: %.tex
	pdflatex -shell-escape -interaction=nonstopmode $< >/dev/null

%.png: %.pdf
	pdftoppm -r 900 -png -singlefile $< $*

clean:
	rm -f *.aux *.log $(FIGS:=.pdf) $(FIGS:=.png)
