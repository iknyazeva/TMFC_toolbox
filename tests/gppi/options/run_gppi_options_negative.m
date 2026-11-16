function results = run_gppi_options_negative(out_dir,spm_dir)

% Negative checks for the gPPI acceleration pooled options slice : every check
% asserts a specific error/warning IDENTIFIER of the options validation
% (evaluation/gppi/OPTIONS_IMPL.md section 4). Cheap by design: the parse
% errors fire before any GLM work (minimal tmfc config); only the rWLS
% guard (needs a real wls parent fixture) and the exact+warm warning
% (needs one full exact run) do real work. Candidate subtree only.
%
% Checks:
%   gPPI: InvalidMethod, InvalidWarmStart, PooledUnsupported (analysis 1
%         and 3), PooledRWLSUnsupported (rwls fixture), PooledFastUnsupported
%         (fast fixture), and the warning PooledOptionsIgnored on an exact
%         run with warm set
%   FIR:  InvalidMethod, PooledUnsupported (analysis 1),
%         PooledFastUnsupported (fir_fast fixture)
%
% FORMAT results = run_gppi_options_negative(out_dir,spm_dir)
%
% License: GPL-3.0-or-later.

if nargin < 2, error('Provide out_dir and spm_dir.'); end
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
candidate = fullfile(repo,'functions');
exact_dir = fullfile(fileparts(here),'exact');

oldpath = path; restore = onCleanup(@() path(oldpath));
addpath(spm_dir); addpath(exact_dir); addpath(here);
set(0,'DefaultFigureVisible','off');
spm('defaults','fmri'); spm_jobman('initcfg');
addpath(genpath(candidate),'-begin');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

results = struct('check',{},'expected',{},'kind',{},'observed',{}, ...
    'pass',{},'error_msg',{});

% ---- fixtures: ar1 (parse checks + warning run) and rwls (wls parent)
[~,fx1] = gppi_exact_fixture('ar1',fullfile(out_dir,'fixture','ar1'));
proj1 = fullfile(out_dir,'neg_ar1');
if exist(proj1,'dir'), rmdir(proj1,'s'); end
mkdir(proj1); copyfile(fullfile(fx1.inputs_dir,'*'),proj1);

% ---- parse-level checks (fail before any subject work)
results = [results; neg_chk(@() tmfc_gPPI(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_method','foo')),1), ...
    'gpPI_invalid_method','tmfc_gPPI:InvalidMethod')];
results = [results; neg_chk(@() tmfc_gPPI(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_reml_warm_start','yes')),1), ...
    'gpPI_invalid_warm','tmfc_gPPI:InvalidWarmStart')];
results = [results; neg_chk(@() tmfc_gPPI(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_method','pooled','analysis',1)),1), ...
    'gpPI_pooled_a1','tmfc_gPPI:PooledUnsupported')];
results = [results; neg_chk(@() tmfc_gPPI(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_method','pooled','analysis',3)),1), ...
    'gpPI_pooled_a3','tmfc_gPPI:PooledUnsupported')];
results = [results; neg_chk(@() tmfc_gPPI_FIR(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_method','bogus')),1), ...
    'fir_invalid_method','tmfc_gPPI_FIR:InvalidMethod')];
results = [results; neg_chk(@() tmfc_gPPI_FIR(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_method','pooled','analysis',1)),1), ...
    'fir_pooled_a1','tmfc_gPPI_FIR:PooledUnsupported')];

% ---- rWLS guard: real wls parent fixture, pooled must be rejected
[~,fxw] = gppi_exact_fixture('rwls',fullfile(out_dir,'fixture','rwls'));
projw = fullfile(out_dir,'neg_rwls');
if exist(projw,'dir'), rmdir(projw,'s'); end
mkdir(projw); copyfile(fullfile(fxw.inputs_dir,'*'),projw);
results = [results; neg_chk(@() tmfc_gPPI(gppi_options_tmfc(fxw,projw, ...
    struct('gppi_method','pooled')),1), ...
    'gpPI_pooled_rwls','tmfc_gPPI:PooledRWLSUnsupported')];

% ---- FAST guard (guarded): FAST parents + pooled are
% rejected loudly; fixture 'fast' (gPPI) and 'fir_fast' (FIR) build real
% FAST parents via the unmodified gppi_exact_fixture
[~,fxf] = gppi_exact_fixture('fast',fullfile(out_dir,'fixture','fast'));
projf = fullfile(out_dir,'neg_fast');
if exist(projf,'dir'), rmdir(projf,'s'); end
mkdir(projf); copyfile(fullfile(fxf.inputs_dir,'*'),projf);
results = [results; neg_chk(@() tmfc_gPPI(gppi_options_tmfc(fxf,projf, ...
    struct('gppi_method','pooled')),1), ...
    'gpPI_pooled_fast','tmfc_gPPI:PooledFastUnsupported')];
[~,fxff] = gppi_exact_fixture('fir_fast',fullfile(out_dir,'fixture','fir_fast'));
projff = fullfile(out_dir,'neg_fir_fast');
if exist(projff,'dir'), rmdir(projff,'s'); end
mkdir(projff); copyfile(fullfile(fxff.inputs_dir,'*'),projff);
results = [results; neg_chk(@() tmfc_gPPI_FIR(gppi_options_tmfc(fxff,projff, ...
    struct('gppi_method','pooled')),1), ...
    'fir_pooled_fast','tmfc_gPPI_FIR:PooledFastUnsupported')];

% ---- exact + warm set: full exact run must WARN (and succeed)
lastwarn('','');   % reset (2-arg form; the 3-arg call was a runtime error)
sub = tmfc_gPPI(gppi_options_tmfc(fx1,proj1, ...
    struct('gppi_reml_warm_start',true)),1);
[~,wid] = lastwarn;
results(end+1) = struct('check','gpPI_exact_warm_warning', ...
    'expected','tmfc_gPPI:PooledOptionsIgnored','kind','warning', ...
    'observed',wid,'pass',strcmp(wid,'tmfc_gPPI:PooledOptionsIgnored') ...
    && sub == 1,'error_msg',''); %#ok<AGROW>
assert(sub == 1,'exact+warm run failed its sub_check');
assert(strcmp(wid,'tmfc_gPPI:PooledOptionsIgnored'), ...
    'expected warning tmfc_gPPI:PooledOptionsIgnored, got: %s',wid);

save(fullfile(out_dir,'results.mat'),'results');
fid = fopen(fullfile(out_dir,'results.json'),'w');
fwrite(fid,jsonencode(results)); fclose(fid);
for i = 1:numel(results)
    fprintf('%s: expected %s, observed %s -> %d\n', ...
        results(i).check,results(i).expected,results(i).observed,results(i).pass);
end
assert(all([results.pass]),'negative check failed');
end

function r = neg_chk(fn,name,expected)
% Run one must-fail call and build its record; the caller asserts all pass.
observed = ''; ok = false; em = '';
try
    fn();                                 % must NOT succeed
    em = 'no error raised';
catch ME
    observed = ME.identifier;
    ok = strcmp(ME.identifier,expected);
    if ~ok, em = getReport(ME,'basic'); end
end
fprintf('%s: expected %s, observed %s -> %d\n',name,expected,observed,ok);
r = struct('check',name,'expected',expected,'kind','error', ...
    'observed',observed,'pass',ok,'error_msg',em);
end
