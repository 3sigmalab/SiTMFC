# Structurally Informed Temporal Multilayer Functional Connectivity (SiTMFC)

Study-specific MATLAB implementation accompanying:

> S. M. Razavi et al., "Structurally Constrained Brain Network Dynamics Reveal Reduced Functional Flexibility in Cocaine Use Disorder."

This repository implements a subject-specific structural-connectivity (SC) filtering step that transforms temporal functional-connectivity (FC) layers before standard multilayer community detection.

## Scope and attribution

The study-specific contribution is the construction of structurally informed functional connectivity (SiFC):

$$
L_{\mathrm{SC}} = I-D^{-1/2}SD^{-1/2},
\qquad
G_{\tau}=(I+\tau L_{\mathrm{SC}})^{-1},
\qquad
\widetilde A^{(s)}=G_{\tau}A^{(s)}G_{\tau}^{\mathsf T}.
$$

Here, $\(S\)$ is a subject-specific SC matrix, $\(D\)$ is its degree matrix, $\(\tau\)$ is the structural smoothing parameter, and $\(A^{(s)}\)$ is the FC matrix for temporal layer $\(s\).$ The transformed matrix $\(\widetilde A^{(s)}\)$ replaces the original FC matrix when constructing the multilayer modularity matrix.

The normalized graph Laplacian and graph spectral filtering are established graph-signal-processing methods. The study-specific contribution is their subject-specific bilateral application to every temporal FC layer and integration into the SiTMFC analysis framework.

**GenLouvain is an external dependency. Its source code and optimization algorithm were not written or modified by the authors of this repository.** GenLouvain remains the work of its original authors and is distributed under its own license.

## Repository file

- `apply_sc_prior.m`: converts an FC matrix or an `L x N x N` FC tensor into structurally informed FC using the corresponding `N x N` subject-specific SC matrix.

## Requirements

- MATLAB with support for local functions in function files.
- FC and SC matrices constructed using the same parcellation and identical node ordering.
- [GenLouvain v2.2](https://github.com/GenLouvain/GenLouvain) for multilayer community detection.

The study configuration used:

| Parameter | Value |
|---|---:|
| Nodes, \(N\) | 200 |
| Temporal layers, \(L\) | 19 |
| SC smoothing, $\(\tau\)$ | 0.3 |
| Modularity resolution, $\(\gamma\)$ | 1 |
| Ordinal coupling, $\(\omega\)$ | 0.5 |
| Stochastic repetitions | 100 |

Sensitivity analyses in the study also considered other values of \(\tau\) and temporal-window length.

## Apply the SC-informed transformation

Expected variables:

- `FC`: either one `N x N` FC matrix or an `L x N x N` tensor of ordered FC layers.
- `SC`: the matching participant's `N x N` structural-connectivity matrix.
- `tau`: a nonnegative smoothing parameter.

```matlab
load('subject_fc.mat','FC');  % L x N x N
load('subject_sc.mat','SC');  % N x N, same node order as FC

tau = 0.3;
[SiFC,G] = apply_sc_prior(FC,SC,tau);
```

The function symmetrizes SC and FC, removes self-connections, retains nonnegative weights, constructs the normalized structural Laplacian, and applies the filter to both endpoints of each FC edge. It assumes quality-controlled inputs: SC and FC must be finite, and every SC node must have positive structural degree.

Setting `\tau = 0` gives `G = I` and provides the positive-only FC baseline after the same common preprocessing.

## Construct the multilayer modularity matrix and run GenLouvain

The transformed layers can be used with GenLouvain's `multiord` helper:

```matlab
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

Ci_all = zeros(N*L,reps);
Q_all  = zeros(reps,1);

for r = 1:reps
    [Ci,Q] = genlouvain(B,10000,0);
    Ci_all(:,r) = Ci;
    Q_all(r) = Q;
end
```

## How to cite

If this SC-informed transformation is used, cite the accompanying study:

> Razavi, S. M., et al. "Structurally Constrained Brain Network Dynamics Reveal Reduced Functional Flexibility in Cocaine Use Disorder." 

If GenLouvain is used, retain its attribution and cite the software as requested by its authors:

> Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha. "A generalized Louvain method for community detection implemented in MATLAB." https://github.com/GenLouvain/GenLouvain (2011-2019).

Also cite the multilayer modularity framework:

> Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., and Onnela, J.-P. "Community structure in time-dependent, multiscale, and multiplex networks." *Science* 328, 876-878 (2010). https://doi.org/10.1126/science.1184819

For graph spectral filtering and graph signal processing, a foundational reference is:

> Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and Vandergheynst, P. "The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains." *IEEE Signal Processing Magazine* 30(3), 83-98 (2013). https://doi.org/10.1109/MSP.2012.2235192

## Code and data availability wording

If the repository contains the complete study-specific analysis pipeline:

> All study-specific code used to generate the reported results is publicly available at https://github.com/3sigmalab/SiTMFC. GenLouvain is an external dependency and is available separately from its original authors at https://github.com/GenLouvain/GenLouvain.

If the repository contains only `apply_sc_prior.m` and selected examples:

> The study-specific MATLAB implementation of the SC-informed functional-connectivity transformation is publicly available at https://github.com/3sigmalab/SiTMFC. GenLouvain is an external dependency and is available separately from its original authors at https://github.com/GenLouvain/GenLouvain.

The study uses the open SUDMEX-CONN dataset cited in the manuscript (Angeles-Valdez, Diego, et al. "The Mexican magnetic resonance imaging dataset of patients with cocaine use disorder: SUDMEX CONN." Scientific data 9.1 (2022): 133.)

