# Change Log
All notable changes to this project will be documented in this file.


## 0.4.1 - 2018-02-08

### Fixed
- Fix incorrect arrangement of node weights in symmetric parameterization.
- Fix the position numbering to use one-based index


## 0.4.0 - 2016-06-03

### Added
- Calculate fitness potential for a sequence. Use `pmrf-infer` command.
- Calculate coevolution scores for an MRF model. Use `pmrf-stat` command.
- Show the MRF model parameters. Use `pmrf-show` command.
- Asymmetric parameterization is added and used as a default.
- `pmrf-build` uses sequence weighting by sequence identity. Use `--seqwt clstr` option. Use `--clstr-maxidt` option to specify the maximum sequence identity between sequence clusters.
- Option to customize L-BFGS parameters.

### Changed
- MRF model file is written to a binary file.
- Use Eigen3 library replacing Blitz++ library. This significantly reduces the computation time.
- The option names and argument values for `pmrf-build` are changed.

### Compatibility Warning
- Not support PMRF v0.2.x MRF file format.


## 0.2.0 - 2016-02-04

### Added
- Build MRF model from an MSA. Use `pmrf-build` command.
- Specify the MRF architecture by using `--edge` option.
- Use L2-regularization for node and edge weights.
- Henikoff's position-based weight for reducing the redundancy in MSA.
