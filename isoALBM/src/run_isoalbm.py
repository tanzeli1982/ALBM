"""Command-line entry point for configurable isoALBM analyses."""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Sequence

from analysis_runner import run_config_file


for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run isoALBM analyses from a TOML configuration file."
    )
    parser.add_argument(
        "--config",
        default="configs/paper.toml",
        help="Path to a TOML config file.",
    )
    parser.add_argument(
        "--condition",
        action="append",
        help="Run only this condition name. Can be supplied more than once.",
    )
    parser.add_argument(
        "--output-dir",
        help="Override the configured output directory. Best for one-run configs.",
    )
    parser.add_argument("--maxiter", type=int, help="Override optimizer maxiter.")
    parser.add_argument("--workers", type=int, help="Override optimizer workers.")
    parser.add_argument("--popsize", type=int, help="Override optimizer popsize.")
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip all configured plot generation.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    run_config_file(
        args.config,
        condition_names=args.condition,
        output_dir_override=args.output_dir,
        maxiter_override=args.maxiter,
        workers_override=args.workers,
        popsize_override=args.popsize,
        no_plots=args.no_plots,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
