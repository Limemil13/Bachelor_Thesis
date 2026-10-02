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

Upload this directory as one project, including `frontmatter/` and its PDF,
or import the prepared ZIP with `main.tex` at its root. Set `main.tex` as the
main document and use pdfLaTeX. This repository has no configured Overleaf
Git remote, so edits made here do not appear automatically in an existing
Overleaf project. To update an existing project without replacing it, upload
the changed `main.tex`, `chapters/*.tex`, and updated expression and synteny
figures into the same relative paths in that project. A project-specific
Overleaf Git remote would enable direct synchronization later.

The current manuscript describes the 98-row SynVoy locus audit as an
exact-assembly GFF/accession assessment with detailed review of 25 priority
exceptions. It does not claim that all 98 rows were independently inspected
visually in NCBI Genome Data Viewer.

Generated PDF and auxiliary build files are outputs, not source files; keep
them out of this directory or place local build output under `build/`.

The ZIP snapshot supplied on 7 September 2026 is an older manuscript export.
It is useful as a front-page and asset reference, but it must not overwrite this
newer source tree.
