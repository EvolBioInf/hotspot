(TeX-add-style-hook
 "asproDoc"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "graphics"
    "eurosym"
    "latexsym"
    "listings"
    "times")
   (TeX-add-symbols
    "be"
    "ee"
    "bi"
    "ei"
    "I"
    "ty"
    "kr"
    "version")))

