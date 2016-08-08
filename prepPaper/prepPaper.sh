#!/bin/bash

rm -rf build/
mkdir -p build/
cp -r ./* ./build
cd build

#make bbl file
pdflatex main; bibtex main;

#do inputs
cp main.tex main.step0.tex

while
	cat main.tex | grep -E "\\input{" | grep -v "^\s*%" > inputs.txt
	python inputs.py
	chmod +x inputs.sh
	./inputs.sh
	[[ -s inputs.txt ]]
do
	continue
done

#do lhcb logo and bibliography
vim -E -s main.tex <<-EOF
	:g/Logo format choice/d
	:g/lhcb-logo.eps/d
	:g/lhcb-logo.pdf/s/^{\\(.*\\)}%$/\\1
	:g/setboolean{inbibliography}{true}/d
	:g/bibliographystyle{LHCb}/d
	:g/bibliography{.*}/r main.bbl
	:g/bibliography{.*}/d
	:wq
EOF
mv figs/lhcb-logo.pdf .

cp main.tex main.step1.tex
#do comments and blank lines
# - delete all comment lines (% preceeded by on whitespace)
# - delete all trailing comments (characters after % unless preceeded by \)
# - remove multiple blank lines
vim -E -s main.tex <<-EOF
	:g/^\\s*%/d
	:g/%/s/\\(.*[^\\\\]\\)%.*/\\1/g
	:%s/\\s\\+$//e
	:%s/\n\\{3,}/\r\r/e
	:wq
EOF

cp main.tex main.step2.tex
#do figures
cat main.tex | grep -E "begin{figure|includegraphics" | grep -v "^\s*%" > figs.txt
vim -E -s figs.txt <<-EOF
	:g/lhcb-logo/d
	:wq
EOF
python figs.py
chmod +x figs.sh
./figs.sh

# tar it up
mkdir -p tar
cp main.tex fig*.pdf lhcb-logo.pdf *.xml *.sty tar/
tar -zcf ../PAPER.tar.gz tar/*

# make the pdf
cd tar
pdflatex main
pdflatex main
cp main.pdf ../../PAPER.pdf
