# Reproducing the LZ77-SSS measurements

This directory reproduces the LZ77 factorization and the (de)compression measurements of the
[paper](https://arxiv.org/abs/2609.30193) (accepted at ALENEX 2027) and regenerates the
figures and the table from the results.

You provide the input texts in [`texts/`](texts/); the scripts run every benchmark and write
the measurements to [`results/`](results/).

## 1. Prerequisites

**Build the project.** Compile it once, following the top-level `README.md`:

```shell
# from the repository root
git clone https://github.com/LukasNalbach/lz77-sss.git
mkdir build && cd build
cmake ..
make
```

CMake checks out the submodules if they are missing, copies the patched submodule sources
into them on every configure run, and builds `alz`.

The tools land in `build/cli/` and `build/bench/`, not in `build/` itself. The three that
matter here are `build/bench/lz77-sss-bench`, `build/bench/zip-bench` and `build/cli/ssszip`.
`measure-all.sh` looks for them at `../build` relative to this folder.

**Compressors on the PATH.** `zip-bench` runs the competitors as external programs, so each
one has to be installed and callable by name:

| Program            | Where it comes from                                              |
| ------------------ | ---------------------------------------------------------------- |
| `lz4`, `gzip`, `bzip2`, `xz`, `zstd` | your distribution (`lz4`, `gzip`, `bzip2`, `xz-utils`, `zstd`) |
| `7z`               | your distribution (`p7zip-full`)                                 |
| `bsc`              | [libbsc](https://github.com/IlyaGrebnov/libbsc), built and installed by hand |
| `alz`              | the `external/alz` submodule, built by this project into `external/alz/alz` |

`alz` is the one exception: it is built from the submodule whenever the benchmark tools are
built, and `zip-bench` runs it from `external/alz/alz`, so it does not have to be installed.
Its own CMake project downloads a few dependencies while configuring, so the first `cmake ..`
needs network access.

`zip-bench` prints `skipping <encoder>: <program> not found` and carries on for any encoder it
cannot find, so a missing one costs you its series in Figure 4 and nothing else.

**GNU time and taskset.** `zip-bench` wraps every competitor in `/usr/bin/time -v` and pins it
to a fixed set of cores with `taskset`. Both are available by default on Ubuntu.

**Memory and disk.** The paper uses 50 GB texts. `zip-bench` writes the compressed file and,
for decompression, the decompressed copy next to the input, so plan for about twice the text
size in free disk space per run.

## 2. Provide the texts

Place each input text in [`texts/`](texts/) under the name the charts select by. All three are
on Zenodo ([10.5281/zenodo.22892236](https://doi.org/10.5281/zenodo.22892236)), compressed
with bsc:

| Text   | File name     | Download | md5 of the download |
| ------ | ------------- | -------- | ------------------- |
| sars2  | `sars2.50Gi`  | [`sars2.50Gi.bsc`](https://zenodo.org/records/22892236/files/sars2.50Gi.bsc?download=1) (533 MB) | `ffc9778eab631d22a365d139359b173d` |
| chr19  | `chr19.50Gi`  | [`chr19.50Gi.bsc`](https://zenodo.org/records/22892236/files/chr19.50Gi.bsc?download=1) (10.6 GB) | `e2432c90a77f4bb2a8d1661bd5099d0c` |
| dewiki | `dewiki.50Gi` | [`dewiki.50Gi.bsc`](https://zenodo.org/records/22892236/files/dewiki.50Gi.bsc?download=1) (214 MB) | `4210641285be377ca26bdb940e85f243` |

Each file decompresses to a text of 50 GB. Decompress it with
[bsc](https://github.com/IlyaGrebnov/libbsc):

```shell
bsc d sars2.50Gi.bsc texts/sars2.50Gi
```

To measure a text of your own, put it in `texts/` and pass its name with `-t` (see below).

## 3. Run

**Everything, for all three texts:**

```shell
./measure-all.sh
```

This truncates `results/results-lz.txt` and `results/results-zip.txt` and then, for each text,
runs three tools that append their own `RESULT` lines:

| Tool             | Invocation                                            | Produces                        |
| ---------------- | ----------------------------------------------------- | ------------------------------- |
| `lz77-sss-bench` | `<text> <max_threads> <result_file>`                   | the factorization runs (Figure 3) |
| `zip-bench`      | `<text> 1 <max_threads> <result_file>`                 | `lz4`, `7z`, `gzip`, `bzip2`, `xz`, `zstd`, `bsc_2047`, `alz_6`, `alz_8` (Figure 4) |
| `ssszip`         | `-t <p> -e <encoder> -r <result_file> -k <text>`       | `ssszip_bsc` and `ssszip_zstd` (Figure 4) |

`ssszip` writes its own `RESULT` lines, which is why it is run separately from `zip-bench`.
`measure-all.sh` runs it with both encoders, at one thread and at `-p` threads, and
decompresses each result again.

Options:

| Flag | Meaning                                       | Default              |
| ---- | --------------------------------------------- | -------------------- |
| `-p` | highest thread count to measure               | 32, as in the paper  |
| `-t` | measure only this text in `texts/`            | all three            |
| `-z` | skip the factorization benchmark              | off                  |
| `-l` | skip the (de)compression benchmark            | off                  |

`zip-bench` walks the thread counts 1, 2, 4, ... up to `-p` on its own. `lz4`, `gzip` and
`bzip2` have no parallel mode and are measured once. Decompression is single-threaded
throughout.

## 4. Output

| File                      | Contents                                                   |
| ------------------------- | ---------------------------------------------------------- |
| `results/results-lz.txt`  | one `RESULT` line per (text, algorithm, thread count)       |
| `results/results-zip.txt` | one `RESULT` line per (text, encoder, thread count, compress/decompress) |

A factorization line looks like this:

```
RESULT text_name=sars2.50Gi n=50000000000 alg=lz77_sss num_threads=1 tau=512 phr_mode=2 fact_mode=1 transf_mode=0 ... num_factors=25448659 comp_ratio=1964.74 time=1613177668077 throughput=30.9947 mem_peak=8722412122
```

`time` is in nanoseconds, `throughput` in MB/s and `mem_peak` in bytes. The charts plot
`(mem_peak + n) / n` on the x-axis, so the text buffer is counted even where `mem_peak` does
not include it. `alg` is `lz77_sss`, `lpf` or `kkp2`; which variant of `lz77_sss` a line
belongs to follows from `phr_mode`, `fact_mode` and `transf_mode`, and the chart selects on
those three.

A compression line looks like this:

```
RESULT text_name=sars2.50Gi type=compress num_threads=1 n=50000000000 encoder=lz4 time=100256801381 throughput=498.719 mem_peak=7200000 bytes_comp=16784613861 comp_ratio=2.97892
```

On decompression `ssszip` names the row after its output file, so those lines carry
`text_name=sars2.50Gi.decompressed`. The charts match with `LIKE 'sars2.50Gi%'` and pick up
both.

[`results-paper/`](results-paper/) holds the measurement data the paper was written from.

## 5. Regenerating the figures and the table

The two figures are generated from the results by
[sqlplot-tools](https://github.com/bingmann/sqlplot-tools). Each file in
[`charts/`](charts/) begins with a `% IMPORT-DATA` line and carries a block of `%% SELECT`
directives. sqlplot-tools loads the named results file into an SQLite database, runs each
query and writes the resulting `\addplot` coordinates into the same file, below its query.

| File                             | Paper float | Reads              |
| -------------------------------- | ----------- | ------------------ |
| `tables/texts.tex`               | Table 2     | nothing, see below |
| `charts/lz.tex`                  | Figure 3    | `results-lz.txt`   |
| `charts/compress-decompress.tex` | Figure 4    | `results-zip.txt`  |

Both charts already contain the numbers of the paper, so they compile without running
sqlplot-tools first.

**Building them.**

```shell
./make-plots.sh --paper     # use results-paper/, the data of the paper
./make-plots.sh             # use results/, your own measurements
```

`make-plots.sh` copies the charts, the table and the results into the directory
`build-plots-paper/` or `build-plots/`, runs sqlplot-tools on every chart there and compiles
all of them into one PDF `plots.pdf`. Each float gets the number it has in the paper.

The files in `charts/` and `tables/` are not modified, so you can compare your numbers with
the numbers of the paper:

```shell
diff -u charts/lz.tex build-plots/charts/lz.tex
```

Running `./make-plots.sh --paper` produces files that are identical to the ones in `charts/`.

**What you need.** sqlplot-tools with the SQLite backend, either on your `PATH` or in the
variable `SQLPLOT_TOOLS`, and pdflatex with the packages `pgfplots` and `subcaption`. On
Ubuntu:

```shell
sudo apt install cmake libboost-all-dev libsqlite3-dev libpq-dev \
                 texlive-latex-recommended texlive-pictures texlive-science
git clone https://github.com/bingmann/sqlplot-tools.git
cd sqlplot-tools && mkdir build && cd build && cmake .. && make
export SQLPLOT_TOOLS=$PWD/src/sqlplot-tools
```

`libpq-dev` is in this list because the CMake file of sqlplot-tools searches for PostgreSQL
even if you only build the SQLite backend. If your sqlplot-tools has both backends, it tries
PostgreSQL first and then uses SQLite. It prints the line `Connection to PostgreSQL failed`
when it does this. That line can be ignored.

[`plot-styles.tex`](plot-styles.tex) contains the colors, the markers, the axis settings and
the annotation macros the charts use, taken from the paper's preamble. Include this file if
you want to use one of the charts in a different document.

**If you measure your own text.** Every chart selects the three texts of the paper by name,
for example `WHERE text_name='sars2.50Gi'`. If you measured a different text, these queries
find no rows. `make-plots.sh` then prints the name of the file, keeps the version with the
numbers of the paper, and continues. To plot your own text, replace the text names in the `%%`
query block of that file.
