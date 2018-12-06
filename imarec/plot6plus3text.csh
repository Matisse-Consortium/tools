#!/bin/tcsh
# create reconstruction overview page

set out = plot6plus3text.tex

cat <<EOF >$out
\documentclass[a4paper]{report}
\usepackage{graphicx}
\usepackage{color}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\usepackage{txfonts}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\topmargin-3.0cm
\textwidth200mm
\textheight290mm
\oddsidemargin-1.5cm
\begin{document}
\unitlength1mm
\pagestyle{empty}
%
EOF
set indent = `awk 'NF>1{printf "l|"}'        param.txt`
set head   = `awk 'NF>1{printf "%s & ",$1}' param.txt`
set line   = `awk 'NF>1{for (i=2;i<=NF;i++)printf "%s ",$i; printf "& "}' param.txt`
set superr = `tail -2 param.txt | head -1`
set targ   = `tail -1 param.txt`

cat <<EOF >>$out
\begin{tabular}{|${indent}l|}
\hline
$head:q \\\\
\hline
$line:q \\\\
\hline
\end{tabular}

\begin{picture}(200,70)
  \put(0,0) {\includegraphics[height=70mm,angle=0]{t5.eps}} 
  \put(75,70){\includegraphics[height=100mm,angle=270]{t66.eps}}
\end{picture}
EOF

set t1 = `head -1 results.txt | tail -1`
set t2 = `head -2 results.txt | tail -1`
set t3 = `head -3 results.txt | tail -1`
set t4 = `head -4 results.txt | tail -1`
set indent = "|ll|ll|ll|lll|l|"
#set head   = "$t1[1] & ResRatio & $t2[1] & ResRatio & $t3[1] & ResRatio & $t4[1] & $t4[2] & $t4[3] & \\"
set head   = ($t1[1] $t2[1] $t3[1] $t4[1] $t4[2] $t4[3])
set line1  = "$t1[2] & $t1[3]   & $t2[2] & $t2[3]   & $t3[2] & $t3[3]   & $t4[4] & $t4[5] & $t4[6] & with uv weight\\"
set line2  = "       &          & $t2[4] & $t2[5]   & $t3[4] & $t3[5]   & $t4[7] & $t4[8] & $t4[9] & without uv weight\\"

cat <<EOF >>$out

\vspace*{4mm}
\begin{tabular}{${indent}}
\hline
\$\chi^{2}\$ $head[1] & ResRatio & \$\chi^{2}\$ $head[2] & ResRatio & \$\chi^{2}\$ $head[3] & ResRatio & $head[4] & $head[5] & $head[6] & \\\\
\hline
$line1:q
$line2:q 
\hline
\end{tabular}

\begin{picture}(200,85)
  \put( 0,0){\includegraphics[height=80mm,angle=0]{t1.eps}} \put( 2,75){\large\bf\textcolor{white} {startima (unconv.) / model (unconv.)}} 
  \put(85,0){\includegraphics[height=80mm,angle=0]{t2.eps}} \put(87,75){\large\bf\textcolor{white} {priorima (unconv.) / model (conv.)}}   \put(87,2){\large\bf\textcolor{white} {$superr}}
  % NE-Grafik
  \put(160,3) {\textcolor{white} {\linethickness{0.3mm} \vector(-1,0){10}}}
  \put(160,3) {\textcolor{white} {\linethickness{0.3mm} \vector(0,1){10}}}
  \put(147,2) {\large\bf\textcolor{white}{E}}
  \put(159,14) {\large\bf\textcolor{white}{N}}

\end{picture}

\begin{picture}(200,85)
  \put( 0,0){\includegraphics[height=80mm,angle=0]{t3.eps}} \put( 2,75){\large\bf\textcolor{white} {reconstruction (unconvolved)}}
  \put( 2,2){{\large\bf\textcolor{white} {$targ}}}
  \put(85,0){\includegraphics[height=80mm,angle=0]{t4.eps}} \put(87,75){\large\bf\textcolor{white} {reconstruction (convolved)}}   \put(87,2){\large\bf\textcolor{white} {$superr}}
\end{picture}

\newpage

  \vspace*{-10mm}
  \includegraphics[height=125mm,angle=270]{tt1.eps} \\\\

  \vspace{-5mm}
  \includegraphics[height=125mm,angle=270]{tt2.eps} \\\\

  \vspace{-5mm}
  \includegraphics[height=125mm,angle=270]{tt3.eps}

\end{document}
EOF

latex $out;
dvips $out:r.dvi

exit
