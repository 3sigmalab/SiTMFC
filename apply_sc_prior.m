function [SiFC, G] = apply_sc_prior(FC, SC, tau)
%APPLY_SC_PRIOR Construct structurally informed functional connectivity.
%
%   [SiFC,G] = APPLY_SC_PRIOR(FC,SC,TAU) applies the subject-specific
%   structural graph filter used in:
%
%       "Structurally Constrained Brain Network Dynamics Reveal Reduced
%       Functional Flexibility in Cocaine Use Disorder."
%
%
%   INPUTS
%   FC    N-by-N FC matrix or L-by-N-by-N tensor of temporal FC layers.
%   SC    N-by-N subject-specific structural-connectivity matrix.
%   tau   structural smoothing parameter.
%
%   OUTPUTS
%   SiFC  Structurally informed FC with the same dimensions as FC.
%   G     Structural graph filter G = (I + tau*L_SC)^(-1).
%
%   METHOD
%       L_SC   = I - D^(-1/2)*SC*D^(-1/2)
%       A_SiFC = G*A*G'
%
%   Because G is symmetric, G*A*G' is equivalent to G*A*G.
%
%   Example:
%       tau = 0.3;
%       [SiFC,G] = apply_sc_prior(FC,SC,tau);
%
%   GenLouvain is an external dependency available from:
%       https://github.com/GenLouvain/GenLouvain

    narginchk(3,3);

    validateattributes(SC, {'numeric'}, ...
        {'real','2d','square','nonempty'}, mfilename, 'SC', 2);
    validateattributes(tau, {'numeric'}, ...
        {'real','scalar','finite','nonnegative'}, mfilename, 'tau', 3);

    SC = double(SC);
    N = size(SC,1);

    if any(~isfinite(SC(:)))
        error('apply_sc_prior:NonfiniteSC', ...
            'SC contains NaN or Inf values.');
    end

    is_single_layer = ismatrix(FC);
    if is_single_layer
        if ~isequal(size(FC), [N N])
            error('apply_sc_prior:FCSize', ...
                'A single FC layer must have size N-by-N, matching SC.');
        end
        n_layers = 1;
    else
        if ndims(FC) ~= 3 || size(FC,2) ~= N || size(FC,3) ~= N
            error('apply_sc_prior:FCSize', ...
                'FC must be N-by-N or L-by-N-by-N, with node dimensions matching SC.');
        end
        n_layers = size(FC,1);
    end

    if ~isnumeric(FC) || ~isreal(FC) || any(~isfinite(FC(:)))
        error('apply_sc_prior:InvalidFC', ...
            'FC must be a real, finite numeric matrix or tensor.');
    end

    SC = (SC + SC.')/2;
    SC(SC < 0) = 0;
    SC(1:N+1:end) = 0;

    degree = sum(SC,2);
    if any(degree <= 0)
        error('apply_sc_prior:ZeroDegreeSC', ...
            'Every SC node must have positive structural degree.');
    end

    Dinv2 = diag(1 ./ sqrt(degree));
    Lsc = eye(N) - Dinv2*SC*Dinv2;
    Lsc = (Lsc + Lsc.')/2;

    [U,Lambda] = eig(Lsc);
    lambda = max(real(diag(Lambda)),0);

    response = 1 ./ (1 + tau*lambda);
    G = U*diag(response)*U.';
    G = real((G + G.')/2);

    if is_single_layer
        SiFC = filter_one_layer(FC,G,N);
    else
        SiFC = zeros(size(FC),'double');
        for s = 1:n_layers
            A = squeeze(FC(s,:,:));
            SiFC(s,:,:) = filter_one_layer(A,G,N);
        end
    end
end

function Aout = filter_one_layer(A,G,N)
    A = double(A);
    A = (A + A.')/2;
    A(1:N+1:end) = 0;
    A(A < 0) = 0;

    Aout = G*A*G.';
    Aout = real((Aout + Aout.')/2);
    Aout(1:N+1:end) = 0;
    Aout(Aout < 0) = 0;
end
