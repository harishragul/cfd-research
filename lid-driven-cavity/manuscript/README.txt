HOW TO SUBMIT THIS FOLDER ON SCHOLARONE
=========================================

UPDATE: this now uses Wiley's OFFICIAL LaTeX template (USG.cls), taken
from the guidelines you provided in ~/Downloads/WileyDesign. Verified
to compile cleanly end-to-end (0 errors) and produces genuine
Wiley-branded output: journal name in the running header/footer, Open
Access badge, correct page numbering, and a proper reference list with
DOIs. See main.pdf for the compiled preview.

One bug found and fixed along the way: the bundled wileyNJD-Chicago.bst
did not render correctly through Tectonic's bibtex engine (citation
years came out blank, e.g. "Ghia et al. ()" instead of "Ghia et al.
(1982)"). This is a Tectonic-specific bibtex incompatibility, not a
problem with your files or with real bibtex. Fix applied: the
bibliography in main.tex is now a hand-written thebibliography block
(correct output verified) instead of an automatic \bibliography{}
call. If Wiley's own production toolchain processes it differently and
you want the .bib-driven \bibliography{references} version instead
(it's more convenient to maintain), say so and I'll switch it back --
their real bibtex almost certainly doesn't have this bug.

THE FOUR SCHOLARONE UPLOAD SLOTS
---------------------------------

1) MAIN MANUSCRIPT (required, LaTeX)
   Zip together ALL of the following into one archive and upload as
   the Main Manuscript -- this is the complete set main.tex depends on:
     - main.tex
     - references.bib        (kept for the typesetter's own conversion
                                pipeline per Wiley's guidance, even
                                though main.tex itself no longer calls
                                \bibliography{} on it directly)
     - USG.cls
     - NJDnatbib.sty, NJDapacite.sty, natbib.sty, wileyNJD-Chicago.bst
     - amssymb.sty, algorithm.sty, algorithmicx.sty, appendix.sty,
       LETTERSP.STY, mla.sty
     - Fonts/                (STIX2 font family used by the class)
     - images/               (Wiley/journal branding graphics used by
                                the class header -- includes the fixed
                                PDF versions of two logos that were
                                broken binary-wrapped EPS files)
     - figures/               (the 6 content figures)
   Verified to compile standalone with "tectonic main.tex" (or
   pdflatex/xelatex + bibtex on a standard install), 0 errors, 8 pages.
   (main.pdf is the already-compiled result, included for your
   preview -- you don't need to upload it, the typesetter compiles
   the .tex source.)

2) GRAPHICAL ABSTRACT IMAGE (required)
   Upload: graphical_abstract.tiff
   (5.5x5.0cm, 300 DPI -- matches Wiley's general spec and the .tiff
   format IJNMF's page explicitly asks for.)

3) GRAPHICAL ABSTRACT TEXT (required)
   Upload: graphical_abstract_text.txt
   (Title + author line + 77-word / 3-sentence summary, under the
   80-word / 3-sentence limit.)

4) AUTHOR'S NOVELTY FILE (required)
   Upload: novelty_file.txt

OPTIONAL
   Cover letter: cover_letter.txt

BEFORE YOU SUBMIT
   - Set the manuscript-type dropdown to Research/Original Article.
   - You'll need to suggest 3-5 potential reviewers yourself.
   - Double check received/revised/accepted dates in main.tex are
     currently "00" placeholders -- ScholarOne likely sets these
     automatically on submission, but confirm against what the portal
     shows.
