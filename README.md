# Surface hopping documentation

## Requirements
- `lapack`. If not available with -llapack, update the makefile to your actual location.

## Compilation 
There are two compilation modes:
- `make debug` (default): No parallelization, many check flags, and if the code breaks, it will show at what line.
- `make fast`: Compiles with -fopenmp (parallelization) and -O2 flags.

## Run
`hopping.x input.in` 
where `input.in` is a file with input arguments.

Normally, the outputs are 
- `popd1.out` Diabatic Pi populations
- `popd2.out` Diabatic P populations
- `popa1.out` Adiabatic Pi populations
- `popa2.out` Adiabatic P populations

(For Tully models with the input argument `branching`, the outputs
are instead `branch1.out` and `branch2.out` for the Pi and P branching
probabilities)

## Available models
- tul = Tully models 
- sbm = Spin-boson model
- fmo = Fenna Matthews Olson complex
- pyr = pyrazine
- dma = DMABN
- ful = fulvene
- so2 = sulfur dioxide
- rhe = rhenium complex

## Input arguments
General arguments:
- model : three-character code for the model (see above)
- ntraj : number of trajectories
- n     : number of modes per bath (sbm/fmo), number of modes (pyr, choose 3/24), model number (tul, choose 1/2/3)
- dt    : time step (in fs or reduced units, depending on model)
- nt    : number of time steps
- method: f for FSSH/Ehrenfest, m for MASH
- hop_opt: Hopping option (used if method==f). 
  - d = FSSH with hopping probability from NACV
  - t = FSSH with hopping probability from TDC with finite differences 
  - e = Ehrenfest dynamics
- res_opt: Rescaling option in FSSH: 
  - d = NACV with mass-scaled quantities
  - h = NACV with Hammes-Schiffer & Tully's 1994 formula
  - v,p = Entire velocity/momentum
- dec_opt: Decoherence option in FSSH:
  - e,c = Energe-based decoherence correction
  - i = instaneous decoherence
  - o = Gaussian overlap decoherence rate
  - v = Gaussian overlaps including only velocity difference part
- dec_param : Decoherence parameter (E0 in a.u. for dec_opt=e/c or r0 for dec_opt=i/o/v)
- nuc_opt: Nuclear sampling
  - wig = Wigner (default)
  - cla = Classical Boltzmann (recommended for sbm/fmo)

System-specific arguments:
- temp : Temperature in kelvin (fmo) or reduced units (sbm)
- eps : epsilon (sbm)
- lambda : reorganization energy (sbm)
- omegac : characteristic frequency (sbm)
- Delta : diabatic coupling (sbm)

NOTE: the initial diabatic state is hard-coded for each model in main.f90

## Examples
- Fulvene with INT and r0=0.01
- SBM (column 3) with GONT and r0=0.01