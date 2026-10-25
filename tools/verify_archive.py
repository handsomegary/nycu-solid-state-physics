"""Check chapter completeness, portable inputs, English text, and optional PDFs."""
import argparse
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent.parent
NON_ENGLISH = re.compile('[' + chr(0x3400) + '-' + chr(0x4dbf)
                         + chr(0x4e00) + '-' + chr(0x9fff)
                         + chr(0xf900) + '-' + chr(0xfaff) + ']')
LOCAL_PATH = re.compile(r'[A-Za-z]:[\\/]Users[\\/]|/h[o]me/|/U[s]ers/')
TEXT_EXTENSIONS = {'.tex', '.md', '.json', '.py', '.yml', '.yaml', '.txt'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', action='store_true', help='Also inspect the built PDFs; requires pypdf')
    args = parser.parse_args()
    exams = json.loads((ROOT / 'exams.json').read_text(encoding='utf-8'))
    ids = [item['id'] for item in exams]
    assert len(ids) == len(set(ids)) == 10
    for kind, prefix in [('questions','q'), ('solutions','s')]:
        master = (ROOT / kind / 'main.tex').read_text(encoding='utf-8')
        for exam_id in ids:
            relative = f'exams/{exam_id}/{kind}.tex'
            content = (ROOT / relative).read_text(encoding='utf-8')
            assert f'../{relative}' in master, relative
            assert f'\\label{{{prefix}:{exam_id}}}' in content, relative
            assert '\\chapter{' in content, relative
    sources = json.loads((ROOT / 'source-coverage.json').read_text(encoding='utf-8'))
    records = sources['files']
    assert len(records) == 50
    assert len({row['id'] for row in records}) == 50
    assert sum(row['format'] == 'photograph' for row in records) == 49
    assert sum(row['format'] == 'pdf' for row in records) == 1
    assert sum('duplicate_of' in row for row in records) == 1
    for row in records:
        assert row['role'] and row['exam_ids']
        assert set(row['exam_ids']) <= set(ids)
    for path in ROOT.rglob('*'):
        relative = path.relative_to(ROOT)
        if any(part in {'.git', '.build', '__pycache__'} for part in relative.parts):
            continue
        assert not NON_ENGLISH.search(relative.as_posix()), str(relative)
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
    print('PASS: 10 paired chapters, all 50 inputs accounted for, portable English sources')
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
