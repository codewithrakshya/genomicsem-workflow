"""Shared dated results paths for Python stages and the workflow launcher."""
from datetime import date, datetime
import os
from pathlib import Path
from zoneinfo import ZoneInfo


def run_date():
    value = os.environ.get('GENOMICSEM_RUN_DATE') or datetime.now(ZoneInfo('America/Los_Angeles')).date().isoformat()
    if date.fromisoformat(value).isoformat() != value:
        raise ValueError('GENOMICSEM_RUN_DATE must be YYYY-MM-DD')
    return value


def results_dir(root):
    return Path(root) / 'results' / run_date()


if __name__ == '__main__':
    print(run_date())
