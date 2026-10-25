"""Build the two XeLaTeX books, or an individual exam for review.

All intermediate files are retained under .build. Use --publish only after
review to copy complete-book PDFs to their public locations.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', default='xelatex')
    parser.add_argument('--kind', choices=['questions', 'solutions', 'both'], default='both')
    parser.add_argument('--exam', help='Build only this exam module for review')
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    engine = shutil.which(args.engine)
    if not engine:
        parser.error('XeLaTeX was not found; install a TeX distribution or supply --engine')
    metadata = json.loads((ROOT / 'exams.json').read_text(encoding='utf-8'))
    ids = {item['id'] for item in metadata}
    if args.exam and args.exam not in ids:
        parser.error('Unknown exam id')
    if args.exam and args.publish:
        parser.error('--publish is reserved for the two complete books')
    version = subprocess.run([engine, '--version'], capture_output=True, text=True,
                             errors='replace', timeout=30, check=True)
    miktex = 'miktex' in (version.stdout + version.stderr).lower()
    results = []
    for kind in (['questions', 'solutions'] if args.kind == 'both' else [args.kind]):
        if args.exam:
            output = ROOT / '.build/preview' / args.exam / kind
            output.mkdir(parents=True, exist_ok=True)
            source = output / 'main.tex'
            prefix = ROOT / 'tex/preamble.tex'
            chapter = ROOT / 'exams' / args.exam / f'{kind}.tex'
            source.write_text('\\documentclass[11pt,a4paper,oneside]{book}\n'
                              f'\\input{{{prefix.as_posix()}}}\n'
                              '\\begin{document}\n'
                              f'\\input{{{chapter.as_posix()}}}\n'
                              '\\end{document}\n', encoding='utf-8')
        else:
            output = ROOT / '.build' / kind
            output.mkdir(parents=True, exist_ok=True)
            source = ROOT / kind / 'main.tex'
        command = [engine] + (['--disable-installer'] if miktex else [])
        command += ['-interaction=nonstopmode', '-halt-on-error', '-file-line-error',
                    f'-output-directory={output}', source.name]
        for number in (1, 2):
            run = subprocess.run(command, cwd=source.parent, capture_output=True,
                                 text=True, errors='replace', timeout=240)
            (output / f'pass-{number}.txt').write_text(run.stdout + run.stderr, encoding='utf-8')
            if run.returncode:
                print('\n'.join((run.stdout + run.stderr).splitlines()[-35:]))
                raise SystemExit(f'Build failed; see {output}')
        pdf = output / 'main.pdf'
        if not pdf.is_file():
            raise SystemExit('No PDF was generated')
        log = (output / 'main.log').read_text(encoding='utf-8', errors='replace')
        warnings = re.findall(r'^(?:Overfull|Underfull|LaTeX Warning|Package .* Warning).*$', log, re.M)
        results.append({'kind': kind, 'exam': args.exam, 'warnings': warnings})
        print(f'PASS {kind}: {args.exam or "complete book"}; {len(warnings)} layout/reference warnings')
        if args.publish:
            shutil.copyfile(pdf, ROOT / kind / f'{kind}.pdf')
    report = ROOT / '.build' / ('last-preview.json' if args.exam else 'books-build.json')
    report.write_text(json.dumps(results, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
