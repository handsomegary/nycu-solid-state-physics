"""Run every published mathematical verification script with this interpreter."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent


def main():
    scripts = sorted(path for folder in ('tools', 'exams')
                     for path in (ROOT / folder).rglob('verify_*.py')
                     if path.name != 'verify_archive.py')
    if len(scripts) != 5 or ROOT / 'tools/verify_2026.py' not in scripts:
        raise SystemExit('Expected five mathematical verification scripts, including verify_2026.py')
    for script in scripts:
        print(f'Checking {script.relative_to(ROOT).as_posix()}', flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print(f'PASS: {len(scripts)} mathematical verification scripts')


if __name__ == '__main__':
    main()
