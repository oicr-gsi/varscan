# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-04-21
[GRD-795](https://jira.oicr.on.ca/browse/GRD-795) - expanded built-in documentation

## [2.5.2] - 2025-03-10
## Added
- Added vcfCombine task to create a final bgzipped vcf with both indels and cnv

## [2.5.1] - 2024-06-27
## Added
- Added vidarr file labels

## [2.5.0] - 2024-06-25
### Added
[GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [2.4.1] - 2024-02-12
### Added
- Added minMemory parameter to tasks using RAM scaling (to handle small chromosomes better)

## [2.4.0] - 2023-12-16
### Added
- Added scaling RAM assignment by chromosome for scattered tasks

## [2.3.0] - 2023-06-21
### Added
- Assembly-specific settings go into wdl

## [2.2.5] - 2023-06-08
### Changed
- Increment version to install varscan_matched in Vidarr

## [Unreleased] - 2023-05-16
### Changed
- updated README, explained passing interval file better

## [2.2.4] - 2022-05-25
### Changed
- Change most output types to non-optional. Classified as a bug fix. Filtered outputs are still optional!

## [2.2.3] - 2022-02-02
### Added
- Add timeout and modules param to interval-processing task so that we can use modularized data

## [2.2.2] - 2021-06-02
### Fixed
- This version increment was to resolve installation issues in vidarr

## [Unreleased] - 2021-11-10
### Fixed
[GP-2891](https://jira.oicr.on.ca/browse/GP-2891) Making RT tests more robust

## [2.2.1] - 2021-02-01
### Fixed
- Increment version to avoid overlap with a compromized installation

## [2.2]   - 2021-01-15
### Changed
- Re-implemented to cut on merging pileups and scatter variant-calling tasks

## [2.1.1] - 2020-04-06
### Added
- Added timeout setting to mpileup concatenation step. 5 hr is not enough for some data

## [2.1] - 2020-04-02
### Changed
- Introduced compression of mpileup data to conserve disk space

## [2.0] - 2020-01-31
### Changed
- Converting to wdl workflow, adding more output types to mesh well with sequenza

## [1.0] - 2017-08-16
### Added
- Initial implementation as a stand-alone workflow


