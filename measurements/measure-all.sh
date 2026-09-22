#!/usr/bin/env bash
#
# Runs the measurements of the paper and writes results/results-lz.txt and
# results/results-zip.txt, the two files the charts read.
#
#   ./measure-all.sh                      # all three texts, 1..32 threads
#   ./measure-all.sh -p 16                # use at most 16 threads
#   ./measure-all.sh -t sars2.50Gi        # a single text
#   ./measure-all.sh -z                   # skip the factorization benchmark
#   ./measure-all.sh -l                   # skip the (de)compression benchmark
#
# The texts are expected in texts/ under the names the charts select by
# (sars2.50Gi, chr19.50Gi, dewiki.50Gi); see replicate.md section 2.
#
# This truncates the two result files before it starts. Each of the three tools appends its
# own RESULT lines:
#   lz77-sss-bench  the factorization runs (Figure 3)
#   zip-bench       lz4, 7z, gzip, bzip2, xz, zstd, bsc and alz (Figure 4)
#   ssszip          its own compression and decompression runs (Figure 4)

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

BUILD=${BUILD:-../build}
TEXTS=texts
RESULTS=results

MAX_THREADS=32
ONLY_TEXT=""
run_lz=1
run_zip=1

while getopts "p:t:zlh" opt; do
    case "$opt" in
        p) MAX_THREADS=$OPTARG ;;
        t) ONLY_TEXT=$OPTARG ;;
        z) run_lz=0 ;;
        l) run_zip=0 ;;
        h) awk 'NR>2 && /^#/ { sub(/^# ?/, ""); print; next } NR>2 { exit }' \
               "${BASH_SOURCE[0]}"; exit 0 ;;
        *) exit 1 ;;
    esac
done

if [ -n "$ONLY_TEXT" ]; then
    TEXT_LIST=("$ONLY_TEXT")
else
    TEXT_LIST=(sars2.50Gi chr19.50Gi dewiki.50Gi)
fi

# CMake puts the benchmark tools in build/bench/ and the CLI tools in build/cli/
LZ_BENCH="$BUILD/bench/lz77-sss-bench"
ZIP_BENCH="$BUILD/bench/zip-bench"
SSSZIP="$BUILD/cli/ssszip"

for t in "$LZ_BENCH" "$ZIP_BENCH" "$SSSZIP"; do
    if [ ! -x "$t" ]; then
        echo "error: $t not found -- build the project first (see replicate.md section 1)" >&2
        exit 1
    fi
done

mkdir -p "$RESULTS"
: > "$RESULTS/results-lz.txt"
: > "$RESULTS/results-zip.txt"

for text in "${TEXT_LIST[@]}"; do
    file="$TEXTS/$text"

    if [ ! -s "$file" ]; then
        echo "warning: $file missing, skipping $text" >&2
        continue
    fi

    if [ "$run_lz" = 1 ]; then
        echo ">>> [$text] LZ77 factorization, up to $MAX_THREADS threads" >&2
        "$LZ_BENCH" "$file" "$MAX_THREADS" "$RESULTS/results-lz.txt" \
            || echo "warning: lz77-sss-bench failed on $text" >&2
    fi

    if [ "$run_zip" = 1 ]; then
        echo ">>> [$text] competitors, 1..$MAX_THREADS threads" >&2
        # zip-bench walks 1, 2, 4, ... up to MAX_THREADS itself; the sequential encoders
        # (lz4, gzip, bzip2) are run once. It needs /usr/bin/time, taskset and the encoder
        # binaries on the PATH -- see replicate.md section 1.
        "$ZIP_BENCH" "$file" 1 "$MAX_THREADS" "$RESULTS/results-zip.txt" \
            || echo "warning: zip-bench failed on $text" >&2

        # ssszip writes its own RESULT lines, so it is run separately from zip-bench.
        # Figure 4 uses ssszip_bsc; ssszip_zstd is measured as well because the paper's
        # results-zip.txt contains it.
        for enc in bsc zstd; do
            for p in 1 "$MAX_THREADS"; do
                echo ">>> [$text] ssszip -e $enc, $p threads" >&2
                "$SSSZIP" -t "$p" -e "$enc" -r "$RESULTS/results-zip.txt" -k "$file" \
                    || { echo "warning: ssszip -e $enc failed on $text" >&2; continue; }

                # decompression is single-threaded; the output name carries the
                # .decompressed suffix the charts' LIKE 'text%' queries expect
                "$SSSZIP" -d -o "$file.decompressed" -r "$RESULTS/results-zip.txt" \
                    "$file.ssszip.$enc" \
                    || echo "warning: ssszip -d -e $enc failed on $text" >&2
                rm -f "$file.ssszip.$enc" "$file.decompressed"
            done
        done
    fi
done

echo
echo "wrote $RESULTS/results-lz.txt  ($(grep -c '^RESULT' "$RESULTS/results-lz.txt" 2>/dev/null || echo 0) rows)"
echo "wrote $RESULTS/results-zip.txt ($(grep -c '^RESULT' "$RESULTS/results-zip.txt" 2>/dev/null || echo 0) rows)"
echo "run ./make-plots.sh to turn them into the figures"
