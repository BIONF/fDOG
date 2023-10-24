# Changelog

## [Dev version]

### Added
-
### Changed
-
### Fixed
-

## [0.1.23] - 2023.08.16

### Added
- Option to NOT adding all search taxa (`--notAddingTaxa`); OFF by default,
i.e. all search taxa will be present in phyloprofile output
- Check invalid min-/max rank for referece species (specified by --minDist and --maxDist).
If the specified ranks (by default, `--minDist genus` `--maxDist kingdom`) are not available,
the next valid ranks will be suggested (or automatically applied if default ranks are used)
- Check if seed sequence cannot be retrieved by blast. Return with the blast command

### Changed
-
### Fixed
- Fixed issue with long directory path for FASTA36 v36.3.8g Dec 2017

## [0.1.12] - 2023.03.10

### Added
- Option to not check annotations for fdog.checkData (option `--ignoreAnno`)
- Check compatibility between blastp and blast DBs

### Changed
- Work with MUSCLE v5.1
- Replace MuscleCommandline and MafftCommandline by subprocess.run

### Fixed
-
