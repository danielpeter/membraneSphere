#!/bin/sh

/bin/rm -rf *.dvi >  /dev/null
/bin/rm -rf *.log >  /dev/null
/bin/rm -rf *.out >  /dev/null
/bin/rm -rf *.aux >  /dev/null
/bin/rm -rf *.toc >  /dev/null
/bin/rm -rf *.blg >  /dev/null
/bin/rm -rf *.bbl >  /dev/null
/bin/rm -rf *.lof >  /dev/null
/bin/rm -rf *.lot >  /dev/null
/bin/rm -rf *.plt >  /dev/null
/bin/rm -rf *.fff >  /dev/null
/bin/rm -rf *.ttt >  /dev/null
/bin/rm -rf *.tit >  /dev/null
/bin/rm -rf *.spl >  /dev/null
/bin/rm -rf *.brf >  /dev/null
/bin/rm -rf *.mtc* >  /dev/null
/bin/rm -rf *.maf >  /dev/null
/bin/rm -rf *.idx >  /dev/null

pdflatex --shell-escape user_manual
bibtex user_manual
pdflatex --shell-escape user_manual
pdflatex --shell-escape user_manual
pdflatex --shell-escape user_manual
pdflatex --shell-escape user_manual

/bin/rm -rf *.dvi >  /dev/null
/bin/rm -rf *.log >  /dev/null
/bin/rm -rf *.out >  /dev/null
/bin/rm -rf *.aux >  /dev/null
/bin/rm -rf *.toc >  /dev/null
/bin/rm -rf *.blg >  /dev/null
/bin/rm -rf *.bbl >  /dev/null
/bin/rm -rf *.lof >  /dev/null
/bin/rm -rf *.lot >  /dev/null
/bin/rm -rf *.plt >  /dev/null
/bin/rm -rf *.fff >  /dev/null
/bin/rm -rf *.ttt >  /dev/null
/bin/rm -rf *.tit >  /dev/null
/bin/rm -rf *.spl >  /dev/null
/bin/rm -rf *.brf >  /dev/null
/bin/rm -rf *.mtc* >  /dev/null
/bin/rm -rf *.maf >  /dev/null
/bin/rm -rf *.idx >  /dev/null

if [ -e user_manual.pdf ]; then
  echo
  echo
  mv -v user_manual.pdf ../membraneSphere-user_manual.pdf
fi

