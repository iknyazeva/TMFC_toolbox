function fx = gppi_exact_case_meta(name)

% Case-static metadata for gppi_exact_compare: mirrors the fx fields of
% gppi_exact_fixture WITHOUT building any data (pure function of the case
% name; no MATLAB GLM work). Used by verify_gppi_exact_saved to validate
% already saved outputs of run_gppi_exact_tests.
%
% FORMAT fx = gppi_exact_case_meta(name)
%   name - fixture case name (see gppi_exact_fixture)
%
% License: GPL-3.0-or-later

cs = struct('cvi','AR(1)','method','gPPI','analysis',2,'window',32,'bins',8, ...
    'maxmem',2^31);
switch name
    case 'ar1'
    case 'fir_ar1',   cs.method = 'fir';
    case 'none',      cs.cvi = 'None';
    case 'fast',      cs.cvi = 'FAST';
    case 'rwls',      cs.cvi = 'wls';
    case 'fir_none',  cs.method = 'fir'; cs.cvi = 'None';
    case 'fir_fast',  cs.method = 'fir'; cs.cvi = 'FAST';
    case 'fir_rwls',  cs.method = 'fir'; cs.cvi = 'wls';
    case 'ar1_s2v',   cs.analysis = 3;
    case 'ar1_a1',    cs.analysis = 1;
    case 'fir_ar1_s2v', cs.method = 'fir'; cs.analysis = 3;
    case 'fir_ar1_a1',  cs.method = 'fir'; cs.analysis = 1;
    case 'none_s2v',  cs.cvi = 'None'; cs.analysis = 3;
    case 'none_a1',   cs.cvi = 'None'; cs.analysis = 1;
    case 'fir_none_s2v', cs.method = 'fir'; cs.cvi = 'None'; cs.analysis = 3;
    case 'fir_none_a1',  cs.method = 'fir'; cs.cvi = 'None'; cs.analysis = 1;
    case 'fast_s2v',   cs.cvi = 'FAST'; cs.analysis = 3;
    case 'fast_a1',    cs.cvi = 'FAST'; cs.analysis = 1;
    case 'rwls_s2v',   cs.cvi = 'wls';  cs.analysis = 3;
    case 'rwls_a1',    cs.cvi = 'wls';  cs.analysis = 1;
    case 'fir_fast_s2v', cs.method = 'fir'; cs.cvi = 'FAST'; cs.analysis = 3;
    case 'fir_fast_a1',  cs.method = 'fir'; cs.cvi = 'FAST'; cs.analysis = 1;
    case 'fir_rwls_s2v', cs.method = 'fir'; cs.cvi = 'wls';  cs.analysis = 3;
    case 'fir_rwls_a1',  cs.method = 'fir'; cs.cvi = 'wls';  cs.analysis = 1;
    case 'ar1_oldvoi'                                   % same compare contract
    otherwise, error('Unknown case %s',name);
end
if strcmp(cs.method,'fir'), fxfunc = 'tmfc_gPPI_FIR'; meth_dir = 'gPPI_FIR';
else, fxfunc = 'tmfc_gPPI'; meth_dir = 'gPPI'; end
roi_names = {'ROI_1','ROI_2','ROI 3'};
fx = struct('set_name','gppi_fx','subject','sub-01','roi_names',{roi_names}, ...
    'session_mode','strict','method',cs.method,'func',fxfunc, ...
    'meth_dir',meth_dir,'analysis',cs.analysis,'cvi',cs.cvi, ...
    'window',cs.window,'bins',cs.bins,'maxmem',cs.maxmem,'nRun',2, ...
    'expected_batches',numel(roi_names),'expected_matrices',4, ...
    'expected_betas_per_roi',4);   % strict mode: A/B x 2 runs = 4
end
