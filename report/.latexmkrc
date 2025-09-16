$aux_dir = 'build';
$out_dir = 'build';

$pdf_mode = 1;

$pdflatex = 'pdflatex -interaction=nonstopmode -output-directory=build %O %S';
$bibtex   = 'bibtex %B';
$makeindex= 'makeindex %O -o %D %S';
$dvipdf   = 'dvipdf %O %S';
