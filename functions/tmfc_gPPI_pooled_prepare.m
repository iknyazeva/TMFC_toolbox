function [Cy,h0,keep,info] = tmfc_gPPI_pooled_prepare(job,SPM1,Y,UFp,warm)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% One-subject preparation for the gPPI pooled-Cy path (gppi_method=
% 'pooled', analysis 2): ONE whole-brain BOLD pass over the concatenated
% scan list of all analyzed runs (tmfc_gPPI_pooled_job +
% tmfc_LSS_load_data, called unmodified), the shared run-pooled Cy
% (tmfc_gPPI_pooled_cy) and, when warm start is requested, the h0 cache
% from a cold first-seed pooled ReML (OPTIONS_DESIGN Revision R4: h0 =
% hyperparameters of the FIRST estimated seed; the cache lives within one
% tmfc_gPPI/tmfc_gPPI_FIR call and is never reused across subjects).
% Computing h0 here, BEFORE the seed loop, keeps the warm semantics
% identical for the serial and parfor seed loops (h0 is broadcast, every
% seed including the first is warm-started; the first seed is estimated
% twice - once cold for h0, once warm in the loop).
%
% FORMAT [Cy,h0,keep,info] = tmfc_gPPI_pooled_prepare(job,SPM1,Y,UFp,warm)
%
%   job  - whole-brain load job (tmfc_gPPI_pooled_job; caller must have
%          rejected SPM_concat == 1 beforehand)
%   SPM1 - STRIPPED first-seed SPM (xX.W / xVi.V removed by the caller)
%   Y    - preloaded ROI data stack of the caller ([nScan_total x nROI])
%   UFp  - voxel-pooling threshold p (SPM stats.fmri.ufp)
%   warm - logical: estimate the h0 cache from a cold first-seed ReML
%
%   Cy   - shared run-pooled covariance [nScan_total x nScan_total]
%   h0   - warm-start hyperparameters ([] when warm is false)
%   keep - logical mask of retained Vi components ([] when warm is false)
%   info - evidence struct for the per-subject log/parsave summary
%
% License: GPL-3.0-or-later. Derived from the TMFC toolbox (R. Masharipov).

D = tmfc_LSS_load_data(job);
[Cy,nPool] = tmfc_gPPI_pooled_cy(SPM1,D,UFp);

h0 = []; keep = [];
if warm
    [~,~,h0] = tmfc_gPPI_pooled_estimate(SPM1,Y,Cy);
    % keep = [] ALWAYS (owner decision s16g, variant b): the h0 source is
    % a SIBLING SEED, and the LSS keep=(h0~=0) filter (parent->trial
    % semantics) would hard-drop components that seed 1 estimated as ~0
    % from every other seed's ReML, changing the model subspace
    % (catastrophic on FAST, s16f). Warm here means INITIALIZATION ONLY
    % (R4: same ML estimator, different start); h0(p) inside reml_block
    % applies the start values to the seed's own components.
    keep = [];
end

info = struct('nPool',nPool,'nScan',size(Cy,1),'warm',warm, ...
    'cy_trace',trace(Cy),'cy_diag_mean',mean(diag(Cy)));
end
