# report/.latexmkrc

# Directory for auxiliary files (aux, log, toc, etc.)
$aux_dir = 'build';

# Directory for output files
$out_dir = 'build';

# Ensure PDF goes to the project root instead of build/
$pdf_mode = 1;

# Customize pdflatex command to output to build/ but generate PDF in root
$pdflatex = 'pdflatex -interaction=nonstopmode -synctex=1 -output-directory=build %O %S';
$bibtex   = 'bibtex %B';
$makeindex= 'makeindex %O -o %D %S';
$dvipdf   = 'dvipdf %O %S';