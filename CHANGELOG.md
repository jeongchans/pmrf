# Change Log
All notable changes to this project will be documented in this file.


## Unreleased

### Added
- Calculate pseudolikelihood value for an aligned sequence. Use `pmrf-infer` command.
- Calculate coevolution scores for an MRF model. Use `pmrf-stat` command.
- Show the MRF model parameters. Use `pmrf-show` command.

### Changes
- MRF model file is written to a binary file.
- Use Eigen3 library replacing Blitz++ library. This significantly reduces the computation time.


## 0.2.0 - 2016-02-04

### Added
- Build MRF model from an MSA. Use `pmrf-build` command.
- Specify the MRF architecture by using `--edge` option.
- Use L2-regularization for node and edge weights.
- Henikoff's position-based weight for reducing the redundancy in MSA.
