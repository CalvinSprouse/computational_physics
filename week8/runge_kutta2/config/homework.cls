% a custom document class used for homework

% custom class header
\NeedsTeXFormat{LaTeX2e}

% make sure the name matches the file hierarchy used
\ProvidesClass{config/homework}[2023/10/04 Custom Homework Class]

% Define a base class (article in this case) and pass options to it
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax

% Load the base class with the passed options
\LoadClass{article}

% Page configurations (you can customize these)
\linespread{1.3}

\setlength{\topmargin}{0pt}

\addtolength{\topmargin}{-\headheight}
\addtolength{\topmargin}{-\headsep}

\setlength{\textheight}{8.9in}
\setlength{\oddsidemargin}{0pt}
\setlength{\evensidemargin}{\oddsidemargin}
\setlength{\marginparwidth}{0.5in}

\pagestyle{empty}

\setlength{\textwidth}{6.25in}

\endinput
