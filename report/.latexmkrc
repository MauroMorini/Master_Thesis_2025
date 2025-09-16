$pdf_mode = 1;

# Everything (auxiliary files, PDF, synctex) goes into build/
$aux_dir = 'build';
$out_dir = 'build';

# Use default pdflatex (latexmk passes -synctex automatically)
$pdflatex = 'pdflatex -interaction=nonstopmode -synctex=1 -output-directory=build %O %S';
