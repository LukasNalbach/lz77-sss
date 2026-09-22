#!/usr/bin/env bash
#
# Regenerates the figures and the table of the paper from the measurement data and builds
# them into a single PDF.
#
#   ./make-plots.sh --paper     # from results-paper/, the data the paper was written from
#   ./make-plots.sh             # from results/, your own measurements
#   ./make-plots.sh --no-pdf    # run sqlplot-tools only, skip pdflatex
#
# Everything is staged in build-plots-paper/ resp. build-plots/; the files in charts/ and
# tables/ are never modified, so you can diff your numbers against the paper's:
#
#   diff -u charts/lz.tex build-plots/charts/lz.tex
#
# Needs sqlplot-tools with the SQLite backend (on the PATH or in $SQLPLOT_TOOLS) and a
# pdflatex with pgfplots and subcaption. See replicate.md section 5.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

use_paper=0
make_pdf=1

for arg in "$@"; do
    case "$arg" in
        --paper)  use_paper=1 ;;
        --no-pdf) make_pdf=0 ;;
        -h|--help) awk 'NR>2 && /^#/ { sub(/^# ?/, ""); print; next } NR>2 { exit }' \
                       "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown option: $arg (try --help)" >&2; exit 1 ;;
    esac
done

if [ "$use_paper" = 1 ]; then
    RESULTS=results-paper
    STAGE=build-plots-paper
else
    RESULTS=results
    STAGE=build-plots
fi

if [ ! -d "$RESULTS" ]; then
    echo "error: $RESULTS/ does not exist" >&2
    [ "$use_paper" = 1 ] || echo "       run ./measure-all.sh first, or use --paper" >&2
    exit 1
fi

# ---------- locate sqlplot-tools ----------
SQLPLOT=${SQLPLOT_TOOLS:-}
if [ -z "$SQLPLOT" ]; then
    for c in sqlplot-tools sqlplot-tools-sqlite3; do
        if command -v "$c" > /dev/null 2>&1; then SQLPLOT=$(command -v "$c"); break; fi
    done
fi
if [ -z "$SQLPLOT" ] || [ ! -x "$SQLPLOT" ]; then
    echo "error: sqlplot-tools not found" >&2
    echo "       put it on your PATH or set SQLPLOT_TOOLS to the binary (see replicate.md)" >&2
    exit 1
fi

# ---------- stage ----------
rm -rf "$STAGE"
mkdir -p "$STAGE/results"
cp -r charts tables plot-styles.tex plots.tex "$STAGE/"
cp "$RESULTS"/*.txt "$STAGE/results/"

# ---------- run sqlplot-tools over every chart and table ----------
cd "$STAGE"
stale=()

for f in charts/*.tex tables/*.tex; do
    # files without an IMPORT-DATA block hold no generated numbers (the texts table is
    # written by hand), so there is nothing for sqlplot-tools to do
    grep -q '^% IMPORT-DATA' "$f" || continue

    if ! "$SQLPLOT" "$f" > "${f%.tex}.sqlplot.log" 2>&1; then
        echo "warning: sqlplot-tools failed on $f, keeping the paper's version" >&2
        echo "         (see $STAGE/${f%.tex}.sqlplot.log)" >&2
        stale+=("$f")
    fi
done

if [ ${#stale[@]} -gt 0 ]; then
    echo >&2
    echo "the following files still hold the paper's numbers:" >&2
    printf '  %s\n' "${stale[@]}" >&2
    echo "Every chart selects the paper's three texts by name (WHERE text_name='sars2.50Gi'" >&2
    echo "and so on). Measuring a text of your own leaves those queries without rows. To plot" >&2
    echo "your own text, replace the text names in the %% query block of the file." >&2
    echo >&2
fi

# ---------- build the PDF ----------
if [ "$make_pdf" = 1 ]; then
    if ! command -v pdflatex > /dev/null 2>&1; then
        echo "warning: pdflatex not found, skipping the PDF" >&2
    else
        for i in 1 2; do
            pdflatex -interaction=nonstopmode -halt-on-error plots.tex > plots.build.log 2>&1 || {
                echo "error: pdflatex failed -- see $STAGE/plots.build.log" >&2
                exit 1
            }
        done
        echo "wrote $STAGE/plots.pdf"
    fi
fi

echo "regenerated charts and tables are in $STAGE/"
