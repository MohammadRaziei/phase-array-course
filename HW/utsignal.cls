\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{utsignal}

\LoadClass[a4paper, dvipsnames]{article}
\RequirePackage{kvoptions}
\RequirePackage{float}
\RequirePackage[inline]{enumitem}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{pstricks}
\RequirePackage{graphicx}
\RequirePackage{geometry}
\RequirePackage{tikz}
\RequirePackage{tikzscale}
\RequirePackage{tocbibind}
\RequirePackage{listings}
\RequirePackage[exponent-to-prefix]{siunitx}
\RequirePackage{fontspec}
\RequirePackage[normalem]{ulem}
\RequirePackage{array}
\RequirePackage[colorlinks, citecolor=Green, linkcolor=Red, filecolor=Cyan, urlcolor=Magenta]{hyperref}
% \RequirePackage{xepersian}
\usepackage[ugly]{nicefrac}


\DeclareDefaultOption{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessKeyvalOptions*

\makeatletter
\AtBeginDocument{
    \tolerance=5000
    \setlength{\parskip}{3ex plus1mm minus1mm}
    \settextfont{BNazanin}
    \setlatintextfont[Scale=0.8]{DejaVu Serif}
    \setmonofont[Scale=0.8]{DejaVu Sans Mono}
    \renewcommand{\labelitemi}{$\circ$}
    \renewcommand{\labelitemii}{$\diamond$}
    \definecolor{codegreen}{rgb}{0,0.6,0}
    \definecolor{codegray}{rgb}{0.5,0.5,0.5}
    \definecolor{codepurple}{rgb}{0.58,0,0.82}
    \definecolor{backcolour}{rgb}{0.95,0.95,0.92}

    \definecolor{identifiercolor}{rgb}{.4,.6,.56}
    \definecolor{stringcolor}{gray}{0.5}
    \definecolor{inactivecolor}{rgb}{0.15,0.15,0.5}
	\definecolor{bluekeywords}{rgb}{0.13,0.13,1}
	\definecolor{greencomments}{rgb}{0,0.5,0}
	\definecolor{redstrings}{rgb}{0.9,0,0}
	
    \lstdefinestyle{matlabstyle}{
    	backgroundcolor=\color{backcolour},   
    	commentstyle=\color{codegreen},
    	keywordstyle=\color{bluekeywords},
    	numberstyle=\tiny\color{codegray},
    	stringstyle=\color{redstrings},
    	basicstyle=\footnotesize,
    	breakatwhitespace=false,         
    	breaklines=true,                 
    	captionpos=b,                    
    	keepspaces=true,                 
    	numbers=left,                    
    	numbersep=5pt,                  
    	showspaces=false,                
    	showstringspaces=false,
    	showtabs=false,                  
    	tabsize=2,
    	inputpath=codes/
    }
    
    \lstset{style=matlabstyle}
    
    \renewcommand{\maketitle}[0]{
        \clearpage
        {
            \begin{tikzpicture}[remember picture, overlay]
            \def\logosize{2cm}
            \def\widelogosize{2.9cm}
            \def\pcenter{\textwidth / 2 - \widelogosize / 4}
            \draw[thick] (-\widelogosize / 2,-0.5) -- (\textwidth,-0.5);
            \node [anchor=center, inner sep=0]  at (\pcenter,\logosize / 2 + \baselineskip) 
            {\Large \rl{\textbf{\@course}}};
            \node [anchor=center, inner sep=0]  at (\pcenter,\logosize / 3.1)
            {\large \rl{\@title}};
            \node [anchor=east, inner sep=0]  at (\textwidth,-0.8) 
            {\small\rl{استاد: \textbf{\@lecturer}}};
            \if\@deadline
            \else
            \node [anchor=center, inner sep=0]  at (\pcenter,-0.1)
            {\footnotesize\rl{مهلت تحویل: \@deadline}};
            \fi
            \end{tikzpicture}
        }

    }

    \hypersetup{pdftitle=\@course\ :: \@title}
    \hypersetup{pdfauthor=\@author}
}

\definecolor{light-gray}{gray}{0.8}
\newcommand{\course}[1]{\def\@course{#1}}
\newcommand{\lecturer}[1]{\def\@lecturer{#1}}
\newcommand{\deadline}[1]{\def\@deadline{#1}}
\deadline{}
\newenvironment{inlinelist}[1][,\ ]{\begin{itemize*}[label={}, afterlabel={}, itemjoin={#1}, itemjoin*={\ and\ }]}{\end{itemize*}}

\makeatother
