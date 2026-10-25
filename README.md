# NYCU Solid State Physics and Theory

English LaTeX reconstructions of eleven Solid State Theory and Solid State Physics exam sets, with a separate book of independently derived, detailed solutions.

- **[Question collection (PDF)](questions/questions.pdf)** — for attempting the problems.
- **[Detailed solutions (PDF)](solutions/solutions.pdf)** — derivations, assumptions, diagrams, and corrections to unreliable reference answers.
- Editable entry points: [questions/main.tex](questions/main.tex) and [solutions/main.tex](solutions/main.tex).

This is an independent study archive, not an official university publication or answer key. The exams span different courses and years; their printed titles and dates are preserved. All public prose, diagram labels, code comments, and document metadata are in English.

## Exam index

Each directory contains `questions.tex` and `solutions.tex`. These are chapters included by the two complete books.

| Exam | Surviving questions | Source status |
| --- | --- | --- |
| [Theory: May 6, 2014](exams/theory-2014-05-midterm) | 1–3 | Only the first printed page is available; visible weights total 60%. |
| [Theory: December 5, 2014](exams/theory-2014-12-midterm) | 1–7 | Both printed pages; printed weights total 110%. |
| [Physics: April 30, 2015](exams/physics-2015-04-midterm) | 1–8 | Both printed pages. |
| [Physics: June 25, 2015](exams/physics-2015-06-final) | 1–9 | Complete supplied printed page; duplicate photograph consolidated. |
| [Theory: April 29, 2016](exams/theory-2016-04-midterm) | 1–12 | All twelve questions on the supplied printed page. |
| [Theory: April 26, 2017](exams/theory-2017-04-midterm) | 1–8 | Both printed pages. |
| [Physics: November 8, 2017](exams/physics-2017-11-midterm) | 1–9 | Both printed pages. |
| [Theory: June 24, 2020](exams/theory-2020-06-final) | 1–12 | Complete supplied printed page; printed weights total 85%. |
| [Theory: April 23, 2021](exams/theory-2021-04-midterm) | 1–4 | Both printed question pages from the supplied PDF. |
| [Theory: Midterm Exam II, May 11, 2026](exams/theory-2026-05-midterm-ii) | 1–7 | Complete questions in a six-page PDF with supplied answers; printed weights total 100%. |
| [Theory: undated fragment](exams/theory-undated-midterm-fragment) | 2(a)–(c), 3–7 | Question 1 and the setup of Question 2 are cropped out; date unavailable. |

The missing May 2014 second page has **not** been invented. Its solution chapter includes clearly marked supplementary exercises inferred from handwritten reference topics. These are not presented as recovered original questions. The undated Question 2 cannot have unique numerical answers without its missing Hamiltonian; its solution explains the limitation and gives a separately labelled conditional method.

## Coverage and editorial choices

The local audit covered **51 supplied files: 49 photographs and two PDFs containing eight and six pages**. One final-exam photograph is an exact duplicate, leaving 50 distinct files by content hash. [source-coverage.json](source-coverage.json) accounts for every input using anonymous source identifiers. It records question pages, reference-only pages, the duplicate, the excluded identity/score cover in one PDF, and the later PDF containing questions with supplied answers. Original local paths and identity mappings are not needed to use this repository.

Printed problems were transcribed and typeset. Supplied solutions, whether handwritten or typeset, were consulted as fallible references and were not treated as authoritative. Figures needed to understand the problems or solutions are reconstructed in editable LaTeX/TikZ.

**Original source PDFs are never uploaded**, including PDFs that contain only the original question sheets. The only PDF files permitted in the public repository are the generated books `questions/questions.pdf` and `solutions/solutions.pdf`. Original photographs, personal identifiers, handwriting scans, and grading covers are also excluded. The archive checker enforces the PDF allowlist in the public working tree and, when Git metadata is available, in tracked paths and available Git history. Local sources and private review material are kept outside the public archive.

Corrections and qualifications include the distinction between Dirac contacts and differentiable stationary points, graphene curvature and density-of-states factors, occupied FCC plane spacing, silicon primitive-cell band counts, mass-weighted phonon eigenvectors, and the difference between exchange-correlation energy and its functional derivative. The 2026 solutions also distinguish bound ionic charge-transfer excitation energies from separated-pair gaps, check the sign of Thomas–Fermi induced charge, and account for the zero-frequency Drude conductivity term in the Kramers–Kronig relations. Printed approximate formulas and unusual premises remain visible in the question book; the solutions explain their implications rather than silently changing the questions.

## Build from a fresh download

Reading the two PDFs requires no build tools. Rebuilding requires **Python 3.10 or newer** and **XeLaTeX**, with TeX Gyre and Latin Modern fonts and the standard LaTeX mathematics/TikZ packages. The build has no dependency on the original source folders.

On Debian or Ubuntu, install:

```sh
sudo apt-get update
sudo apt-get install texlive-xetex texlive-latex-extra texlive-science fonts-texgyre
```

From the repository root:

```sh
python tools/verify_archive.py
python tools/run_checks.py
python tools/build_documents.py --publish
```

The last command runs two LaTeX passes for each book and copies the completed PDFs to `questions/questions.pdf` and `solutions/solutions.pdf`. Intermediate files and compiler logs remain in the ignored `.build/` directory. With MiKTeX, install missing packages beforehand; the script disables automatic installation during compilation. If XeLaTeX is not on `PATH`, add `--engine /path/to/xelatex`.

To build one chapter for review:

```sh
python tools/build_documents.py --exam theory-2020-06-final
```

The preview stays under `.build/preview/`; `--publish` is reserved for the complete books. You can also open either book's `main.tex` in a LaTeX editor and select XeLaTeX.

## Verification

[GitHub Actions](../../actions) runs the archive checks, the numerical checks, and both full-book builds on Linux. To inspect PDF text and metadata locally as well:

```sh
python -m pip install pypdf
python tools/verify_archive.py --pdf
```

The five numerical scripts use the Python standard library. They check selected derivations through independent finite differences, state counting, reciprocal-lattice enumeration, direct Hamiltonian action, density-of-states integration, mass-weighted phonon identities, and the additional 2026 exam calculations. They are reproducible checks of those formulas, not a proof that every conceptual answer is beyond revision. The rendered books are also reviewed for missing glyphs, clipped equations, overlapping objects, and source omissions.

Primary scientific references are linked beside the relevant derivations. The editable sources make assumptions and corrections reviewable without relying on an unseen answer sheet.
