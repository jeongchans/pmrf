# PMRF - Protein Markov random field
This is an evolutionary analysis tool using Markov random field.


## Installation

### Prerequisites
The following libraries are essential for running this software.

- [GSL - GNU Scientific Library] (http://www.gnu.org/software/gsl/)
- [Blitz++] (http://blitz.sourceforge.net)

The following library is required only for test.

- [Google Test] (https://github.com/google/googletest)

### Download
TBA

### Building the software

1. Compile the software.

```
  $ make
```

2. The executable will be placed in the directory where you untarred PMRF.


## Usage
TBA

### Help messages
`pmrf --help` shows the general help page.
`pmrf <COMMAND> --help` shows the help page for the specific command.


## Examples

### Generating MRF model from an MSA
TBA

### Generating MRF model by specifying the network architecture
`pmrf-build` can specify the network architecture of MRF model with `--edge` option.

### Calculating the pseudolikelihood for aligned sequeces
TBA


## Acknowledgements

This software includes and uses [libLBFGS] (http://www.chokkan.org/software/liblbfgs/).
