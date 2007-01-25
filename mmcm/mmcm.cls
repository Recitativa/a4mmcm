\NeedsTeXFormat{LaTeX2e}[1994/12/01]
\ProvidesClass{mmcm}[2005/06/14 v050601, by Mao Ziyang]

\RequirePackage{graphicx} 
\RequirePackage{amsmath}
\RequirePackage{amsfonts} 
\RequirePackage{amssymb}
\RequirePackage{amsbsy}
%\usepackage{amsthm}
\RequirePackage{mathrsfs}
\RequirePackage{color}
\RequirePackage{pifont}

\RequirePackage{theorem}
\RequirePackage{CJK}

\@ifundefined{exebook}{
\LoadClass[11pt,oneside, a4paper]{article}
\newcommand{\docbegin}{\begin{document}\begin{CJK*}{GBK}{song}}
\newcommand{\docend}{\end{CJK*}\end{document}}
{
\theorembodyfont{\normalfont}
\newtheorem{Exec}{{Exercise}}[section]
}
\def\setsection#1{{\setcounter{section}{#1}}}
}{
\LoadClass[11pt,oneside,a4paper]{book}
\newcommand{\docbegin}{}
\newcommand{\docend}{}
\def\section#1{\chapter{#1}}
{
\theorembodyfont{\normalfont}
\newtheorem{Exec}{{Exercise}}[chapter]
}
\def\setsection#1{{\setcounter{chapter}{#1}}}
\def\setsection#1{}
}

\newcommand{\setexenum}[1]{{\setcounter{Exec}{#1}}}
\newcommand{\stepexenum}{{\stepcounter{Exec}}}


%\endlinechar -1

%\topmargin=-15.4mm
\oddsidemargin=-0.4mm
%\evensidemargin=-0.4mm
\textwidth=160mm
%\textheight=247mm

%%\newcommand\sectionformat{\centering}

\def\ae#1#2{\text{a.e.} {#1} \in {#2}}

\def\brn{\mathbb{R}^N}
\def\br{\mathbb{R}}
\def\bn{\mathbb{N}}
\def\bc{\mathbb{C}}
\def\bq{\mathbb{Q}}
\def\bz{\mathbb{Z}}
\def\bbf{\mathbb{F}}
\def\ch{\mathcal{H}}
\def\cb{\mathscr{B}}
\def\cw{\mathcal{W}}
\def\cd{\mathcal{D}}
\def\cc{\mathcal{C}}
\def\ck{\mathcal{K}}
\def\cl{\mathcal{L}}
\def\cs{\mathcal{S}}
\def\cm{\mathcal{M}}
\def\scl{\mathscr{L}}
\def\lpr{L^2_{\br}[0,\pi]}
\def\lpc{L^2_{\bc}[0,\pi]}
\def\lppr{L^2_{\br}[-\pi,\pi]}
\def\lppc{L^2_{\bc}[-\pi,\pi]}
\def\lpf{L^2_{\bbf}[-\pi,\pi]}
\def\cpr{C_{\br}[0,\pi]}
\def\brk{\mathbb{R}^k}
\def\brm{\mathbb{R}^m}
\def\bpk{\mathbb{P}^{k-1}}
\def\ixy{\langle x,y\rangle}
\def\ixa{\langle x,a\rangle}
\def\ixx{\langle x,x\rangle}
\def\iyx{\langle y,x\rangle}
\def\ixz{\langle x,z\rangle}
\def\iyz{\langle y,z\rangle}
\def\cat{\text{\rm cat}}
\def\ve{\varepsilon}
\def\bve{\overline{\ve}}
\def\ox{\overline{X}}
\def\tf{\tilde{f}}
\def\tg{\tilde{g}}
\def\th{\tilde{h}}
\def\tx{\tilde{X}}
\def\ta{\tilde{A}}
\def\vp{\varphi}
\def\tvp{\tilde{\vp}}
\def\wc{\rightharpoonup}
\def\hd{\hat{\delta}}
\def\hve{\hat{\ve}}
\def\unk{u_{n_k}}

\def\b0{{\boldsymbol{0}}}
\def\bu{{\boldsymbol{u}}}
\def\bv{{\boldsymbol{v}}}
\def\bw{{\boldsymbol{w}}}
\def\bx{{\boldsymbol{x}}}
\def\bxd{{\boldsymbol{\dot{x}}}}
\def\bxdd{{\boldsymbol{\ddot{x}}}}


\def\Re{\mathfrak{Re}}
\def\Im{\mathfrak{Im}}

\def\esssup{\text{ess}\,\sup}
\def\sspan{{\text{\rm span}}}

%\newtheorem{theorem}{{Theorem}}[section]
%\newtheorem{proposition}[theorem]{{Proposition}}
%\newtheorem{lemma}[theorem]{{Lemma}}
%\newtheorem{corollary}[theorem]{{Corollary}}
%\newtheorem{definition}{{Definition}}[section]
%\newtheorem{example}{{Example}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Exercises counter and environment
%

%\newcounter{exercises}[section]

%\newenvironment{exercises}[1][\empty]{\par
%  \pushQED{\qed}%
%  \normalfont \topsep6\p@\@plus6\p@\relax
%  \trivlist
%  \stepcounter{exenum}
%  \item[\hskip\labelsep
%        {\textbf{Exercises}}
%        \arabic{section}.\arabic{exenum}
%        %\@ifnotempty{#1}{({#1})}
%        \ifx#1\empty{}\else{({\itshape #1})}\fi: 
%        ]
%        \hspace{0.5em}
%        \ignorespaces
%}{%
%  \popQED\endtrivlist\@endpefalse
%}


\newenvironment{solution}[1][\empty]{\par
  \par%
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \item[\hskip\labelsep
        {\textbf{Solution}} 
        %\@ifnotempty{#1}{({#1})}
        \ifx#1\empty{}\else{({\itshape #1})}\fi: 
        ]
        \hspace{0.5em}
        \ignorespaces
}{%
  \endtrivlist\@endpefalse
}

\newenvironment{proof}[1][\empty]{\par
  \par%
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \item[\hskip\labelsep
        {\textbf{Proof}} 
        %\@ifnotempty{#1}{({#1})}
        \ifx#1\empty{}\else{({\itshape #1})}\fi: 
        ]
        \hspace{0.5em}
        \ignorespaces
}{%
  \endtrivlist\@endpefalse
}

\newcommand{\dinf} {\displaystyle\inf}
\newcommand{\dint} {\displaystyle\int}
\newcommand{\dlim} {\displaystyle\lim}
\newcommand{\dsum} {\displaystyle\sum}
\newcommand{\dsup} {\displaystyle\sup}


\newcommand{\inp}[2]{{\left\langle {#1} ,\;{#2}  \right\rangle}}

\newcommand{\rproof}{{($\Rightarrow$)\quad}}
\newcommand{\lproof}{{($\Leftarrow$)\quad}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Mark Set
%
\newcommand{\red}[1]{{\color{red}{#1}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% LIST ENVIRONMENTS
%
% Alphlist:  A., B., ...
%
\newcounter{Alphnum}
\newenvironment{Alphlist}{\begin{list}{\Alph{Alphnum}.}
        {\usecounter{Alphnum}\setlength{\topsep}{1ex}
        \setlength{\itemsep}{1pt}\setlength{\leftmargin}{4em}}
        \rm}{\end{list}}
%
% alphlist:  (a), (b), ...
%
\newcounter{alphnum}
\newenvironment{alphlist}{\begin{list}{(\alph{alphnum})}
        {\usecounter{alphnum}\setlength{\leftmargin}{1.5em}}
        \rm}{\end{list}}
%
% arablist:  (1), (2), ...
%
\newcounter{arabicnum}
\newenvironment{arablist}{\begin{list}{(\arabic{arabicnum})}
        {\usecounter{arabicnum}}
         \rm}{\end{list}}
%
% romlist: (i), (ii), ...
%
\newcounter{romnum}
\newenvironment{romlist}{\begin{list}{(\roman{romnum})}
    {\usecounter{romnum}\setlength{\topsep}{1pt}
    \setlength{\itemsep}{1pt}}
    \rm}{\end{list}}
%
% Romlist:  I., II., ...
%
\newcounter{Romnum}
\newenvironment{Romlist}{\begin{list}{\Roman{Romnum}.}
    {\usecounter{Romnum}\setlength{\topsep}{1pt}
    \setlength{\itemsep}{1pt}}
    \rm}{\end{list}}




%%\DeclareGraphicsExtensions{.eps,.eps.gz,{}}
%%\DeclareGraphicsRule{*}{eps}{*}{}

\endlinechar `\^^M

%%\RequirePackage[dvips,CJKbookmarks]{hyperref}

\endinput

