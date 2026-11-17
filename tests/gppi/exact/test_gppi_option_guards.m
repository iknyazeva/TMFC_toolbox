function test_gppi_option_guards(functions_dir)
% Preflight test for the exact-only option guards of tmfc_gPPI and
% tmfc_gPPI_FIR (owner review 2, 2026-11-17). SOURCE-READY ONLY: NEVER RUN
% in its authoring unit (zero MATLAB); it executes in the approved
% validation launch. GPL-3.0-or-later.
%
% The guards sit at the very top of both functions, BEFORE any tmfc field
% access, IO, directory creation or SPM initialisation. Therefore the two
% subtests below prove the explicit rejections fire before ANY heavy work:
%   1) gppi_method='pooled'  -> error <fun>:InvalidMethod
%   2) gppi_reml_warm_start=true -> error <fun>:WarmStartRemoved
% A defaults-free config is NOT probed here: the guards fall through by
% design and the normal legacy validation (later in the file) would run;
% the true default behaviour is validated end-to-end by the fixture cases
% (ar1/fir_ar1/... pass no option fields).
%
% The function THROWS on any unexpected outcome (no silent false return).
% No mkdir, no file loads, no SPM calls are made by this test itself.
%
% FORMAT test_gppi_option_guards(functions_dir)
%   functions_dir - the candidate exact-only functions/ folder

if nargin < 1, error('test_gppi_option_guards:inputs', 'provide functions_dir'); end
addpath(functions_dir);
for fun = {'tmfc_gPPI','tmfc_gPPI_FIR'}
    name = fun{1};
    w = which(name);
    assert(~isempty(w) && strncmp(w, [functions_dir filesep], numel(functions_dir)+1), ...
        'test_gppi_option_guards:which', ...
        '%s must resolve under %s (got %s)', name, functions_dir, w);

    % minimal struct: only the two fields the top guard may touch; any
    % behaviour beyond the guard is irrelevant because the guard MUST fire
    mk = @(opt) struct('defaults', opt);

    % 1) pooled rejected explicitly, before any tmfc access beyond defaults
    t = mk(struct('gppi_method','pooled'));
    e = get_err(@() feval(name, t, 1));
    assert(~isempty(e) && strcmp(e.identifier, [name ':InvalidMethod']), ...
        'test_gppi_option_guards:pooled', ...
        '%s: expected %s:InvalidMethod, got %s', name, name, disp_id(e));

    % 2) warm-start field rejected explicitly (even together with exact)
    t = mk(struct('gppi_method','exact','gppi_reml_warm_start',true));
    e = get_err(@() feval(name, t, 1));
    assert(~isempty(e) && strcmp(e.identifier, [name ':WarmStartRemoved']), ...
        'test_gppi_option_guards:warm', ...
        '%s: expected %s:WarmStartRemoved, got %s', name, name, disp_id(e));
end
fprintf('test_gppi_option_guards: PASS (both functions, both guards)\n');
end

function e = get_err(fn)
e = [];
try
    fn();
catch ME
    e = ME;
end
end

function s = disp_id(e)
if isempty(e), s = '(no error)'; else, s = e.identifier; end
end
