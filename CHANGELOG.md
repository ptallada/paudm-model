# Change Log (http://keepachangelog.com/)
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [3.1.0] - 2020-07-20
### Added
- In table image_zp, new column n_stars

### Changed
- In table star_zp, column zp_weight -> chi2

### Removed
- Tables template, template_band, template_fit, star_template_zp

## [3.0.0] - 2019-04-17
### Changed
- Compatible with environment refactoring (v3.0)
- Left only model code - removed common and grid sections

### Removed
- Grid Section
- Common Objects module

## [2.2.0] - 2019-03-29
### Changed
- Python 3.X compatible
- Added install db from environment

## [2.1.0] - 2019-03-21
### Changed
- Minor changes to calibration tables

## [2.0.0] - 2018-01-15
### Added
- New Forced Photometry model (non-sparse)
- Added crosstalk tables
### Removed
- Major changes and cleanup of deprecated tables

## [1.2.0] - 2017-07-26
### Added
- Photometric Calibration model
- Quality control model

## [1.1.0] - 2016-04-28
### Added
- Reference catalogues & tables for MEMBA

## [1.0.0] - 2015-07-27
### Added
- brownthrower tables
- reference tables for Nightly
- fixed with commissioning data

## [0.2.0] - 2013-07-29
### Added
- event jobs tables 
- tables for simulations

## [0.1.0] - 2013-06-13
### Added
- sqlite compatibility
- first set of tables