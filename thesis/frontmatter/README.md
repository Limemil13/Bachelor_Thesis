# University front matter

`university_frontpage_template.pdf` is the official two-page LMU/TUM template
used by `main.tex`.

The source PDF is not a fillable form. main.tex preserves it as a background and
places the thesis metadata over the template placeholders when LaTeX compiles.
Edit the metadata commands near the top of main.tex; do not edit the PDF itself.

If your program distinguishes a formal Supervisor from an Advisor, put the names
in \thesissupervisor and \thesisadvisor respectively. Otherwise leave
\thesisadvisor empty and keep both names in \thesissupervisor.
