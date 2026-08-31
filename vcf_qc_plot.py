#!/usr/bin/env python3
"""
vcf_qc_plot.py

Quick multi-panel QC diagnostic plots for VCF hard-filtering decisions.

It shells out to `bcftools query` to pull QUAL and a set of INFO fields
(one row per variant site), then plots histograms with:
  - the site-level distribution
  - a robust x-axis (auto-trimmed at the 99th percentile so one huge
    outlier doesn't crush the histogram)
  - dashed reference lines for conventional GATK-style hard-filter
    thresholds (edit CUTOFFS below to match your own choices)

It also writes a summary table (min/max/mean/median/quantiles + fraction
of sites that would be flagged by the reference cutoff) to <out>_summary.tsv,
so you have numbers to cite in your methods, not just a picture.

Requirements: bcftools on PATH, python3 with pandas / numpy / matplotlib.

Usage
-----
    python3 vcf_qc_plot.py --vcf variants.vcf.gz --out qc

    # restrict to a subset of fields, e.g. if some INFO tags aren't present
    python3 vcf_qc_plot.py --vcf variants.vcf.gz --out qc \
        --fields QUAL,INFO/DP,INFO/QD,INFO/MQ

    # add extra INFO fields not in the default set
    python3 vcf_qc_plot.py --vcf variants.vcf.gz --out qc \
        --fields QUAL,INFO/DP,INFO/QD,INFO/MQ,INFO/FS,INFO/SOR,INFO/MQRankSum,INFO/ReadPosRankSum

    # only look at a region / BED (e.g. your BUSCO BED) rather than the whole VCF
    python3 vcf_qc_plot.py --vcf variants.vcf.gz --out qc --regions-file busco.bed
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Conventional GATK-style hard-filter reference cutoffs (SNP-oriented).
# These are starting points, not gospel -- adjust to match what you justify
# in methods. `direction` says which side of the line is the "bad" side.
# ---------------------------------------------------------------------------
CUTOFFS = {
    "QUAL": {"value": 30, "direction": "below", "label": "QUAL < 30"},
    "INFO/DP": None,  # depth cutoffs are dataset-specific (coverage-dependent); see --dp-low/--dp-high
    "INFO/QD": {"value": 2, "direction": "below", "label": "QD < 2"},
    "INFO/MQ": {"value": 40, "direction": "below", "label": "MQ < 40"},
    "INFO/FS": {"value": 60, "direction": "above", "label": "FS > 60"},
    "INFO/SOR": {"value": 3, "direction": "above", "label": "SOR > 3"},
    "INFO/MQRankSum": {"value": -12.5, "direction": "below", "label": "MQRankSum < -12.5"},
    "INFO/ReadPosRankSum": {"value": -8, "direction": "below", "label": "ReadPosRankSum < -8"},
}

DEFAULT_FIELDS = ["QUAL", "INFO/DP", "INFO/QD", "INFO/MQ", "INFO/FS", "INFO/SOR"]


def check_bcftools():
    if shutil.which("bcftools") is None:
        sys.exit("Error: bcftools not found on PATH. Load/install bcftools and retry.")


def build_query_format(fields):
    # bcftools query wants %QUAL, %INFO/DP, etc., tab-separated, newline-terminated
    return "\t".join(f"%{f}" for f in fields) + "\n"


def run_bcftools_query(vcf, fields, regions_file=None, targets_file=None):
    fmt = build_query_format(fields)
    cmd = ["bcftools", "query", "-f", fmt]
    if regions_file:
        cmd += ["-R", regions_file]
    if targets_file:
        cmd += ["-T", targets_file]
    cmd += [vcf]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"bcftools query failed:\n{proc.stderr}")
    if not proc.stdout.strip():
        sys.exit("bcftools query returned no rows -- check your VCF/region filters.")
    return proc.stdout


def parse_to_dataframe(raw_tsv, fields):
    from io import StringIO

    df = pd.read_csv(
        StringIO(raw_tsv),
        sep="\t",
        header=None,
        names=fields,
        na_values=[".", "nan", "NA", ""],
    )
    for f in fields:
        df[f] = pd.to_numeric(df[f], errors="coerce")
    return df


def summarize(df, fields, cutoffs, out_prefix):
    rows = []
    for f in fields:
        s = df[f].dropna()
        n_total = len(df[f])
        n_missing = df[f].isna().sum()
        row = {
            "field": f,
            "n_sites": n_total,
            "n_missing": n_missing,
            "min": s.min() if len(s) else np.nan,
            "p01": s.quantile(0.01) if len(s) else np.nan,
            "median": s.median() if len(s) else np.nan,
            "mean": s.mean() if len(s) else np.nan,
            "p99": s.quantile(0.99) if len(s) else np.nan,
            "max": s.max() if len(s) else np.nan,
        }
        cutoff = cutoffs.get(f)
        if cutoff and len(s):
            if cutoff["direction"] == "below":
                flagged = (s < cutoff["value"]).sum()
            else:
                flagged = (s > cutoff["value"]).sum()
            row["cutoff_label"] = cutoff["label"]
            row["n_flagged_by_cutoff"] = flagged
            row["frac_flagged_by_cutoff"] = round(flagged / len(s), 4)
        else:
            row["cutoff_label"] = ""
            row["n_flagged_by_cutoff"] = ""
            row["frac_flagged_by_cutoff"] = ""
        rows.append(row)

    summary = pd.DataFrame(rows)
    out_path = f"{out_prefix}_summary.tsv"
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote summary table: {out_path}")
    print(summary.to_string(index=False))
    return summary


def plot_panels(df, fields, cutoffs, out_prefix, dp_low=None, dp_high=None, log_fields=None):
    log_fields = set(log_fields or [])
    n = len(fields)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.array(axes).reshape(-1)

    for ax, f in zip(axes, fields):
        s = df[f].dropna()
        if len(s) == 0:
            ax.set_title(f"{f} (no data)")
            ax.axis("off")
            continue

        use_log = f in log_fields
        if use_log:
            n_nonpositive = (s <= 0).sum()
            s_plot = s[s > 0]
            if len(s_plot) == 0:
                ax.set_title(f"{f} (no positive values; cannot log-scale)")
                ax.axis("off")
                continue
            if n_nonpositive:
                print(f"Note: {f} has {n_nonpositive} value(s) <= 0 excluded from the "
                      f"log-scaled panel (log undefined at/below zero).")
            # trim at 1st/99th percentile of the positive values, then use
            # log-spaced bins so bin widths are even on a log axis
            lo, hi = s_plot.quantile(0.01), s_plot.quantile(0.99)
            if lo <= 0:
                lo = s_plot.min()
            if lo == hi:
                lo, hi = s_plot.min(), s_plot.max()
            plot_data = s_plot[(s_plot >= lo) & (s_plot <= hi)]
            bins = np.logspace(np.log10(lo), np.log10(hi), 80)
            ax.hist(plot_data, bins=bins, color="#4C72B0", edgecolor="none")
            ax.set_xscale("log")
            ax.set_xlabel(f"{f} (log scale)")
        else:
            # trim x-axis at 1st/99th percentile so outliers don't dominate the view
            lo, hi = s.quantile(0.01), s.quantile(0.99)
            if lo == hi:
                lo, hi = s.min(), s.max()
            plot_data = s[(s >= lo) & (s <= hi)]
            ax.hist(plot_data, bins=80, color="#4C72B0", edgecolor="none")
            ax.set_xlabel(f)

        ax.set_title(f)
        ax.set_ylabel("n sites")

        cutoff = cutoffs.get(f)
        if f == "INFO/DP":
            if dp_low is not None:
                ax.axvline(dp_low, color="red", linestyle="--", linewidth=1.2,
                           label=f"DP < {dp_low}")
            if dp_high is not None:
                ax.axvline(dp_high, color="red", linestyle="--", linewidth=1.2,
                           label=f"DP > {dp_high}")
            if dp_low is not None or dp_high is not None:
                ax.legend(fontsize=8)
        elif cutoff:
            ax.axvline(cutoff["value"], color="red", linestyle="--", linewidth=1.2,
                       label=cutoff["label"])
            ax.legend(fontsize=8)

    # turn off any unused axes
    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle("VCF site-level QC distributions (x-axes trimmed to 1st-99th percentile)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = f"{out_prefix}_qc_panels.png"
    fig.savefig(out_path, dpi=150)
    print(f"Wrote plot: {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="Input VCF/BCF (indexed if using -R/-T).")
    ap.add_argument("--out", required=True, help="Output prefix for plot + summary table.")
    ap.add_argument("--fields", default=",".join(DEFAULT_FIELDS),
                     help=f"Comma-separated bcftools query fields. Default: {','.join(DEFAULT_FIELDS)}")
    ap.add_argument("--regions-file", default=None,
                     help="Optional BED/region file passed to bcftools query -R "
                          "(e.g. your BUSCO BED, to QC only those sites).")
    ap.add_argument("--targets-file", default=None,
                     help="Optional file passed to bcftools query -T (streaming, no index needed).")
    ap.add_argument("--dp-low", type=float, default=None,
                     help="Optional low-DP reference line (dataset/coverage-specific).")
    ap.add_argument("--dp-high", type=float, default=None,
                     help="Optional high-DP reference line (dataset/coverage-specific).")
    ap.add_argument("--log-fields", default="",
                     help="Comma-separated subset of --fields to plot on a log10 x-axis "
                          "(uses log-spaced bins). Good candidates: QUAL, INFO/DP, INFO/FS, "
                          "INFO/SOR -- these are typically right-skewed and strictly positive. "
                          "Values <= 0 are excluded from that panel (log undefined there) and "
                          "reported to stdout. Rank-sum fields (can be negative) generally "
                          "should NOT be log-scaled. Example: --log-fields QUAL,INFO/DP,INFO/FS")
    args = ap.parse_args()

    check_bcftools()

    fields = [f.strip() for f in args.fields.split(",") if f.strip()]

    print(f"Querying {len(fields)} field(s) from {args.vcf} ...")
    raw = run_bcftools_query(args.vcf, fields, args.regions_file, args.targets_file)
    df = parse_to_dataframe(raw, fields)
    print(f"Parsed {len(df)} site(s).")

    log_fields = [f.strip() for f in args.log_fields.split(",") if f.strip()]
    unknown_log_fields = set(log_fields) - set(fields)
    if unknown_log_fields:
        sys.exit(f"--log-fields contains field(s) not in --fields: {sorted(unknown_log_fields)}")

    summarize(df, fields, CUTOFFS, args.out)
    plot_panels(df, fields, CUTOFFS, args.out, dp_low=args.dp_low, dp_high=args.dp_high,
                log_fields=log_fields)


if __name__ == "__main__":
    main()
