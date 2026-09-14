# Bachelor thesis LaTeX project

This is the active manuscript source. Open or upload this directory as one
LaTeX project; `main.tex` is the root document.

The manuscript now uses a findings-first order modelled on the recommended
chair example: Introduction, Results, Discussion, Conclusion, and Methods and
Materials. Results are organized by biological question and evidence layer;
the exact implementation details remain together in the final Methods chapter.

## Where to edit

- `main.tex`: thesis title, author, supervisor/advisor, submission date,
  document layout, and section order.
- `chapters/`: manuscript text, split by thesis section.
- `bibliography.bib`: BibTeX references.
- `figures/`: final figures actually used in the manuscript.
- `tables/`: frozen table exports actually used in the manuscript.
- `supplementary/`: final supplementary files.
- `frontmatter/university_frontpage_template.pdf`: official two-page LMU/TUM
  template required during compilation.

The authoritative computational results stay under `../analyses/`. Do not move
the analysis folders into the LaTeX project. Copy only the final figure or table
chosen for the manuscript into `figures/` or `tables/`, and record its source in
the caption or manuscript notes.

## University front-page integration

The template is already integrated. `main.tex` uses both PDF pages as
backgrounds and replaces the example title, author, supervisor/advisor, and date
with the metadata declared near the top of the file. Do not edit the template
PDF itself.

Before the final submission, replace `\today` in `\thesisdate` with the fixed
official submission date and confirm whether the program wants both names on
the Supervisor line or one on the Advisor line.

## Overleaf

Upload the contents of this directory together, including the `frontmatter/`
folder and its PDF. Set `main.tex` as the main document. The project uses
standard pdfLaTeX-compatible packages, so pdfLaTeX is the safest compiler
choice.

Generated PDF and auxiliary build files are outputs, not source files; keep
them out of this directory or place local build output under `build/`.

The ZIP snapshot supplied on 7 September 2026 is an older manuscript export.
It is useful as a front-page and asset reference, but it must not overwrite this
newer source tree.
