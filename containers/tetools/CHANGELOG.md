# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## 1.5 - 2022-02-09
### Updated
- RepeatModeler: `2.0.2a` -> `2.0.3`
### Changed
- Updated and expanded upon several sections of the README

## 1.4 - 2021-05-03
### Added
- Installed `ps` in the container, for compatibility with Nextflow tracing/metrics
### Updated
- RepeatModeler: `2.0.1` -> `2.0.2`
- UCSC utilities: `v407` -> `v413`

## 1.3.1 - 2021-04-01
### Updated
- RepeatMasker: `4.1.2` -> `4.1.2-p1`

## 1.3 - 2021-03-23
### Added
- Installed `strace` in the container
- Added several of the UCSC utilities (`faToTwoBit`, `twoBitInfo`, `twoBitToFa`).
  These tools are used by some of the scripts in `RepeatModeler/util`.
### Updated
- HMMER: `3.3` -> `3.3.2`
- RepeatMasker: `4.1.1` -> `4.1.2`
- RMBlast: `2.10.0` -> `2.11.0`
    * Upstream BLAST+ 2.11.0 introduced opt-out reporting of usage statistics,
      which is also included in RMBlast 2.11.0. See
      <http://www.repeatmasker.org/RMBlast.html> for more details. To opt out
      of this, set the environment variable `BLAST_USAGE_REPORT=false` (e.g.
      with the `-e` option to `docker run`, or the `--env` option for
      `singularity exec`)
### Fixed
- coseg: set a working path to a perl interpreter in `.pl` scripts

## 1.2 - 2020-09-09
### Added
- TRF: `4.09.1`
### Updated
- LTR_retriever: `2.8` -> `2.9.0`
- MAFFT: `7.453` -> `7.471`
- RepeatMasker: `4.1.0` -> `4.1.1`
### Changed
- dfam-tetools.sh: `--trf_prgm` is no longer necessary and is ignored
### Fixed
- dfam-tetools.sh: Fix compatibility with macOS (fixes #3)
- dfam-tetools.sh: Use the '--init' option when running docker to avoid a zombie apocalypse (fixes #2)

## 1.1 - 2020-01-13
### Updated
- LTR_retriever: `2.7` -> `2.8`
- genometools: `1.5.10` -> `1.6.0`
- HMMER: `3.2.1` -> `3.3`
- RepeatModeler: `2.0` -> `2.0.1`
- RMBlast: `2.9.0-p2` -> `2.10.0`
- MAFFT: `7.429` -> `7.453`

### Changed
- Use the C locale for singularity when running via dfam-tetools.sh,
  silencing some warnings from perl programs.

### Fixed
- Fixed the HMMER path for RepeatMasker, so that -engine hmmer will work properly.

## 1.0 - 2019-11-15
### Added
- First release.
