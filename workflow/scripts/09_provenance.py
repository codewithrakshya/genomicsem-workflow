#!/usr/bin/env python3
"""Record stage invocations and content hashes; verify archived input manifests."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import uuid
from paths import results_dir


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def fingerprint(path, root):
    try:
        name = str(path.relative_to(root))
    except ValueError:
        name = str(path)
    if not path.is_file():
        return {'path': name, 'missing': True}
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    stat = path.stat()
    return {'path': name, 'sha256': digest.hexdigest(), 'bytes': stat.st_size,
            'mtime_ns': stat.st_mtime_ns}


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def start(root, settings_file, stage):
    settings = {r['key']: r['value'] for r in rows(root / settings_file)}
    traits = rows(root / settings['traits_config'])
    models = rows(root / settings['models_config'])
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    directory = root / '.work/runs' / run_id
    directory.mkdir(parents=True)
    code = [root / 'run_pipeline.sh', root / 'Snakefile', root / settings_file,
            root / settings['traits_config'], root / settings['models_config'],
            root / 'config/genomicsem_ref.txt', root / 'config/cohort_metadata.tsv', root / 'config/cohort_overlap.md', root / 'workflow/snakemake/environment.yml']
    code += [root / m['model_file'] for m in models]
    code += list((root / 'workflow/scripts').glob('*.py')) + list((root / 'workflow/scripts').glob('*.R'))
    code += list((root / 'workflow/bin').glob('*'))
    code = sorted(set(code))
    snapshots = directory / 'snapshot'
    snapshot_map = []
    for i, path in enumerate(code):
        if path.is_file():
            target = snapshots / f'{i:03d}_{path.name}'
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            snapshot_map.append({'source': str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
                                 'snapshot': str(target.relative_to(directory))})
    inputs = [root / t['raw_file'] for t in traits]
    for key in ('hm3_reference', 'gwas_reference', 'ld_reference', 'weight_reference'):
        path = root / settings[key]
        inputs.extend(sorted(p for p in path.rglob('*') if p.is_file()) if path.is_dir() else [path])
    # Reporting alone does not read raw GWAS/reference files; do not imply it validates them.
    if stage in ('report', 'html-report', 'setup', 'check', 'validate'):
        inputs = []
    manifest = {'schema_version': 1, 'run_id': run_id, 'stage': stage, 'status': 'running',
                'started_utc': now(), 'project_root': str(root),
                'command': ['./run_pipeline.sh', stage], 'settings': str(settings_file), 'results_directory': str(results_dir(root).relative_to(root)),
                'python': sys.version, 'platform': platform.platform(),
                'environment': {k: os.environ[k] for k in (
                    'GENOMICSEM_RUN_DATE', 'GENOMICSEM_SETTINGS', 'GENOMICSEM_CORES', 'GENOMICSEM_PILOT_SNPS',
                    'GENOMICSEM_REF', 'GENOMICSEM_REINSTALL', 'R_LIBS_USER', 'NSLOTS', 'SLURM_CPUS_PER_TASK') if k in os.environ},
                'code_and_configuration': [fingerprint(p, root) for p in code],
                'input_files': [fingerprint(p, root) for p in sorted(set(inputs))],
                'snapshot_map': snapshot_map,
                'existing_results': [fingerprint(p, root) for p in sorted(results_dir(root).glob('*')) if p.is_file() and p.suffix != '.html']}
    r_code = r'''
.libPaths(c(Sys.getenv("R_LIBS_USER", unset="workflow/R/library"), .libPaths()))
cat(capture.output(sessionInfo()), sep="\n")
cat("\nINSTALLED PACKAGES\n")
x <- installed.packages(fields=c("RemoteSha", "RemoteRepo", "RemoteUsername"))
write.table(x[, c("Package", "Version", "RemoteSha", "RemoteRepo", "RemoteUsername")], stdout(), sep="\t", quote=FALSE, row.names=FALSE)
'''
    try:
        completed = subprocess.run(['Rscript', '--vanilla', '-e', r_code], cwd=root,
                                   capture_output=True, text=True, timeout=60)
        (directory / 'r_environment.txt').write_text(completed.stdout + completed.stderr)
        manifest['r_environment_exit_code'] = completed.returncode
    except (OSError, subprocess.TimeoutExpired) as error:
        manifest['r_environment_error'] = str(error)
    write(directory / 'manifest.json', manifest)
    return directory / 'manifest.json'


def finish(root, manifest_path, code):
    manifest = json.loads(manifest_path.read_text())
    current = [fingerprint(p, root) for p in sorted((root / manifest.get('results_directory', str(results_dir(root).relative_to(root)))).glob('*')) if p.is_file() and p.suffix != '.html']
    before = {r['path']: r for r in manifest['existing_results']}
    manifest.update(status='succeeded' if code == 0 else 'failed', exit_code=code, finished_utc=now())
    manifest['changed_results'] = [r for r in current if r != before.get(r['path'])]
    manifest['note'] = 'Status describes this invocation only. Unchanged or pre-existing results are not attributed to this run. HTML is excluded to avoid a circular report hash.'
    write(manifest_path, manifest)


def verify(root, manifest_path):
    manifest = json.loads(manifest_path.read_text())
    failures = []
    for saved in manifest['code_and_configuration'] + manifest['input_files']:
        actual = fingerprint(root / saved['path'], root)
        if saved.get('missing') or actual.get('sha256') != saved.get('sha256'):
            failures.append(saved['path'])
    if failures:
        raise SystemExit('Missing or changed files:\n' + '\n'.join(failures))
    print('All recorded code, configuration, and input hashes match. R package versions must also match the archived environment.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['start', 'finish', 'verify'])
    parser.add_argument('--project-root', default='.')
    parser.add_argument('--settings', default='config/workflow.tsv')
    parser.add_argument('--stage', default='all')
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--exit-code', type=int, default=0)
    args = parser.parse_args()
    root = Path(args.project_root).resolve()
    if args.action == 'start':
        print(start(root, Path(args.settings), args.stage))
    elif args.action == 'finish':
        finish(root, args.manifest, args.exit_code)
    else:
        verify(root, args.manifest)
