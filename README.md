# Structurally Informed Temporal Multilayer Functional Connectivity (SiTMFC)

## Summary

This repository provides the study-specific MATLAB implementation accompanying:

> S. M. Razavi et al., "Structurally Constrained Brain Network Dynamics Reveal Reduced Functional Flexibility in Cocaine Use Disorder" [1].

SiTMFC applies a subject-specific structural-connectivity (SC) filter to every temporal functional-connectivity (FC) layer before temporal multilayer community detection. The structurally informed functional connectivity (SiFC) transformation is

$$
L_{\mathrm{SC}}=I-D^{-1/2}SD^{-1/2},
$$

$$
G_{\tau}=(I+\tau L_{\mathrm{SC}})^{-1},
$$

$$
\widetilde{A}^{(s)}=G_{\tau}A^{(s)}G_{\tau}^{\mathsf T}.
$$

Here, $S$ is a subject-specific SC matrix, $D$ is its degree matrix, $\tau$ is the structural smoothing parameter, and $A^{(s)}$ is the FC matrix for temporal layer $s$. The transformed matrix $\widetilde{A}^{(s)}$ replaces the original FC matrix when constructing the temporal multilayer modularity matrix.

The normalized graph Laplacian and graph spectral filtering are established graph-signal-processing methods [5]. The study-specific contribution is their subject-specific bilateral application to every temporal FC layer and integration into the SiTMFC framework.

The repository contains:

- `apply_sc_prior.m`: transforms either one `N x N` FC matrix or an `L x N x N` tensor of ordered temporal FC layers using the matching `N x N` subject-specific SC matrix.

## Instruction

### Installation

The required dependencies are:

1. [MATLAB](https://www.mathworks.com/products/matlab.html) [2]. The `apply_sc_prior.m` function uses standard MATLAB functions and does not require an additional MATLAB toolbox.
2. [GenLouvain v2.2](https://github.com/GenLouvain/GenLouvain) [3], including its `multiord` helper, for temporal multilayer community detection based on the Mucha et al. framework [4].

Download or clone this repository:

```bash
git clone https://github.com/3sigmalab/SiTMFC.git
```

Download GenLouvain separately from its original repository. In MATLAB, add both directories and all GenLouvain subdirectories to the MATLAB path:

```matlab
addpath('/path/to/SiTMFC');
addpath(genpath('/path/to/GenLouvain'));
```

GenLouvain includes precompiled MEX files. If they are incompatible with the local MATLAB release or operating system, follow the GenLouvain instructions and run `compile_mex.m` from its `MEX_SRC` directory.

The input data must satisfy the following requirements:

- `FC` is either one `N x N` FC matrix or an `L x N x N` tensor of ordered temporal FC layers.
- `SC` is the matching participant's `N x N` structural-connectivity matrix.
- FC and SC use the same parcellation and identical node ordering.
- FC and SC contain finite values, and every SC node has positive structural degree.
- `tau` is a nonnegative structural smoothing parameter.

### Sample Run

The following example uses the parameter settings from the accompanying study:

| Parameter | Value |
|---|---:|
| Nodes, $N$ | 200 |
| Temporal layers, $L$ | 19 |
| SC smoothing, $\tau$ | 0.3 |
| Modularity resolution, $\gamma$ | 1 |
| Ordinal coupling, $\omega$ | 0.5 |
| Optimization repetitions, $O$ | 100 |

```matlab
% Load one participant's temporal FC tensor and matching SC matrix.
load('subject_fc.mat','FC');  % FC: 19 x 200 x 200
load('subject_sc.mat','SC');  % SC: 200 x 200

% Construct structurally informed FC for every temporal layer.
tau = 0.3;
[SiFC,G] = apply_sc_prior(FC,SC,tau);

% Construct the ordered temporal multilayer modularity matrix.
N = size(SiFC,2);
L = size(SiFC,1);

A = cell(L,1);
for s = 1:L
    A{s} = squeeze(SiFC(s,:,:));
end

gamma = 1;
omega = 0.5;
reps  = 100;

[B,twomu] = multiord(A,gamma,omega); %#ok<ASGLU>

% Perform repeated GenLouvain optimizations.
Ci_all = zeros(N*L,reps);
Q_all  = zeros(reps,1);

for r = 1:reps
    [Ci,Q] = genlouvain(B,10000,0);
    Ci_all(:,r) = Ci;
    Q_all(r) = Q;
end
```

`SiFC` has the same dimensions as `FC`. Each `SiFC(s,:,:)` matrix is the structurally informed FC matrix for temporal layer $s$. `Ci_all` contains the multilayer community assignments from the repeated optimizations.

## Acknowledgement

The subject-specific SC-informed FC transformation in `apply_sc_prior.m` was developed for the accompanying SiTMFC study [1].

## References

1. Razavi, S. M., et al. "Structurally Constrained Brain Network Dynamics Reveal Reduced Functional Flexibility in Cocaine Use Disorder."

2. The MathWorks, Inc. *MATLAB*. Natick, Massachusetts, United States. [https://www.mathworks.com/products/matlab.html](https://www.mathworks.com/products/matlab.html)

3. Jeub, L. G. S., Bazzi, M., Jutla, I. S., and Mucha, P. J. "A generalized Louvain method for community detection implemented in MATLAB." [https://github.com/GenLouvain/GenLouvain](https://github.com/GenLouvain/GenLouvain) (2011--2019).

4. Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., and Onnela, J.-P. "Community structure in time-dependent, multiscale, and multiplex networks." *Science* 328, 876--878 (2010). [https://doi.org/10.1126/science.1184819](https://doi.org/10.1126/science.1184819)

5. Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and Vandergheynst, P. "The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains." *IEEE Signal Processing Magazine* 30(3), 83--98 (2013). [https://doi.org/10.1109/MSP.2012.2235192](https://doi.org/10.1109/MSP.2012.2235192)

6. Angeles-Valdez, D., Rasgado-Toledo, J., Issa-Garcia, V., et al. "The Mexican magnetic resonance imaging dataset of patients with cocaine use disorder: SUDMEX CONN." *Scientific Data* 9, 133 (2022). [https://doi.org/10.1038/s41597-022-01251-3](https://doi.org/10.1038/s41597-022-01251-3)
