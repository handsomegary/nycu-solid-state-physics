"""Check chapter coverage, English sources, the PDF publication policy, and optional PDFs."""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent.parent
NON_ENGLISH = re.compile('[' + chr(0x3400) + '-' + chr(0x4dbf)
                         + chr(0x4e00) + '-' + chr(0x9fff)
                         + chr(0xf900) + '-' + chr(0xfaff) + ']')
LOCAL_PATH = re.compile(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]'
                        r'(?=[^{}\s<>]*[\\/]|[^{}\s<>]*\.[A-Za-z0-9]{1,8}\b)'
                        r'|/h[o]me/|/U[s]ers/')
TEXT_EXTENSIONS = {'.tex', '.md', '.json', '.py', '.yml', '.yaml', '.txt'}
ALLOWED_PDFS = {'questions/questions.pdf', 'solutions/solutions.pdf'}


def check_pdf_paths(paths, location):
    """Reject original-source PDFs even if they have been staged or later deleted."""
    pdfs = {name for name in paths if Path(name).suffix.lower() == '.pdf'}
    unexpected = pdfs - ALLOWED_PDFS
    assert not unexpected, f'Unapproved PDF paths in {location}: {sorted(unexpected)}'
    return pdfs


def check_git_pdf_paths():
    """Inspect available repository history; a ZIP download needs no Git installation."""
    if not (ROOT / '.git').exists():
        print('NOTE: no Git metadata; PDF policy checked in the public working tree only')
        return
    git = shutil.which('git')
    if not git:
        print('NOTE: Git is unavailable; PDF policy checked in the public working tree only')
        return

    def run(*arguments):
        return subprocess.run([git, *arguments], cwd=ROOT, check=True,
                              capture_output=True, text=True, encoding='utf-8').stdout

    tracked = [path for path in run('ls-files', '-z').split('\0') if path]
    check_pdf_paths(tracked, 'tracked files')
    history = [path.strip('\r\n') for path in
               run('log', '--all', '--format=', '--name-only', '--no-renames', '-z').split('\0')
               if path.strip('\r\n')]
    pdfs = check_pdf_paths(history, 'available Git history')
    shallow = run('rev-parse', '--is-shallow-repository').strip() == 'true'
    qualifier = 'available shallow' if shallow else 'all available'
    print(f'PASS: PDF allowlist in tracked files and {qualifier} Git history: {sorted(pdfs)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', action='store_true', help='Also inspect the built PDFs; requires pypdf')
    args = parser.parse_args()
    exams = json.loads((ROOT / 'exams.json').read_text(encoding='utf-8'))
    ids = [item['id'] for item in exams]
    assert len(ids) == len(set(ids)) == 11
    for kind, prefix in [('questions','q'), ('solutions','s')]:
        master = (ROOT / kind / 'main.tex').read_text(encoding='utf-8')
        included_ids = re.findall(r'\\input\{\.\./exams/([^/]+)/' + kind + r'\.tex\}', master)
        assert included_ids == ids, f'{kind} chapter order differs from exams.json'
        for exam_id in ids:
            relative = f'exams/{exam_id}/{kind}.tex'
            content = (ROOT / relative).read_text(encoding='utf-8')
            assert f'../{relative}' in master, relative
            assert f'\\label{{{prefix}:{exam_id}}}' in content, relative
            assert '\\chapter{' in content, relative
    sources = json.loads((ROOT / 'source-coverage.json').read_text(encoding='utf-8'))
    records = sources['files']
    assert sources['source_file_count'] == len(records) == 51
    assert sources['unique_file_hash_count'] == 50
    assert len({row['id'] for row in records}) == 51
    assert sum(row['format'] == 'photograph' for row in records) == 49
    assert sum(row['format'] == 'pdf' for row in records) == 2
    assert sum('duplicate_of' in row for row in records) == 1
    for row in records:
        assert row['role'] and row['exam_ids']
        assert set(row['exam_ids']) <= set(ids)
        if 'duplicate_of' in row:
            assert row['duplicate_of'] in {other['id'] for other in records}
    addition = next(row for row in records if row['id'] == 'S3-051')
    assert addition['format'] == 'pdf' and addition['pages'] == 6
    assert addition['exam_ids'] == ['theory-2026-05-midterm-ii']
    public_paths = []
    for path in ROOT.rglob('*'):
        relative = path.relative_to(ROOT)
        if any(part in {'.git', '.build', '__pycache__'} for part in relative.parts):
            continue
        assert not NON_ENGLISH.search(relative.as_posix()), str(relative)
        if path.is_file():
            public_paths.append(relative.as_posix())
        if not path.is_file() or path.suffix not in TEXT_EXTENSIONS:
            continue
        content = path.read_text(encoding='utf-8')
        assert not NON_ENGLISH.search(content), f'Non-English script in {relative}'
        assert not LOCAL_PATH.search(content), f'Private absolute path in {relative}'
        if path.suffix == '.tex':
            for target in re.findall(r'\\(?:input|include)\{([^}]+)\}', content):
                dependency = path.parent / target
                if not dependency.suffix:
                    dependency = dependency.with_suffix('.tex')
                assert dependency.is_file(), f'Missing input {target} in {relative}'
                assert dependency.resolve().is_relative_to(ROOT.resolve()), str(dependency)
    check_pdf_paths(public_paths, 'public working tree')
    check_git_pdf_paths()
    print('PASS: 11 paired chapters, all 51 inputs accounted for, portable English sources')
    if args.pdf:
        from pypdf import PdfReader
        for kind in ('questions', 'solutions'):
            reader = PdfReader(ROOT / kind / f'{kind}.pdf')
            assert len(reader.pages) >= 10
            content = '\n'.join(page.extract_text() or '' for page in reader.pages)
            metadata = str(reader.metadata)
            assert not NON_ENGLISH.search(content + metadata), kind
            assert not LOCAL_PATH.search(content + metadata), kind
            assert 'NYCU Solid State Physics' in content, kind
            assert reader.metadata.title and kind[:-1].lower() in reader.metadata.title.lower(), kind
            print(f'PASS: {kind}.pdf, {len(reader.pages)} pages, English text and metadata')


if __name__ == '__main__':
    main()
