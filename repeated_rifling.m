function bc = repeated_rifling(varargin)
% REPEATED_RIFLING  Build repeated-shot boundary conditions (BC) for the 2-D solver.
% Usage (case-insensitive):
%   bc = repeated_rifling('Nshots',7,'cool_gap',0.5,'Tamb_C',20,'use_30col',false,'plot',true);
%
% Returns struct bc with fields: t [s], h [Nt x Nz], Tg [K], Ttot [s], meta

% ---------- defaults ----------
cfg.Nshots    = 7;        % number of shots
cfg.gap_s     = 0.5;      % cool-down between shots [s]
cfg.Tamb_C    = 20;       % ambient during gaps [°C]
cfg.use_30col = false;    % keep 30 IB HTCs (true) else downselect to 6
cfg.scale     = [];       % optional per-zone HTC scale (1 x Nz)
cfg.plot      = true;     % quick BC figure
cfg.smoke_test= true;     % run a quiet 2-point test through the solver

% ---------- parse name–value ----------
if mod(numel(varargin),2)~=0, error('Use name–value pairs.'); end
for k=1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'nshots','shots'},           cfg.Nshots = double(val);
        case {'gap_s','cool_gap','gap'},   cfg.gap_s = double(val);
        case {'tamb_c','tamb','ambient_c'},cfg.Tamb_C = double(val);
        case {'use_30col','use30','keep30'}, cfg.use_30col = logical(val);
        case {'scale'},                    cfg.scale = val;
        case {'plot'},                     cfg.plot = logical(val);
        case {'smoke_test','smoketest'},   cfg.smoke_test = logical(val);
        otherwise, error('Unknown option "%s".', key);
    end
end
if cfg.Nshots<1 || fix(cfg.Nshots)~=cfg.Nshots
    error('Nshots must be a positive integer.');
end

% ---------- 1) single-shot interior ballistics ----------
if exist('interiorBallistics.m','file')~=2
    error('interiorBallistics.m not found on path.');
end
run('interiorBallistics.m');  % expects: timeVec (s), hHist (Nt x NzIB), Tgas (K)
if ~exist('timeVec','var') || ~exist('hHist','var') || ~exist('Tgas','var')
    error('interiorBallistics.m must define timeVec, hHist, Tgas.');
end
t1   = timeVec(:);
h1   = hHist;
Tg1K = Tgas(:);
if size(h1,1) ~= numel(t1) || numel(Tg1K) ~= numel(t1)
    error('IB arrays inconsistent: size(hHist,1) and numel(Tgas) must equal numel(timeVec).');
end
NzIB = size(h1,2);

% downselect or keep
if NzIB==30 && ~cfg.use_30col
    pick = [1,5,10,16,25,30];
    h1   = h1(:,pick);
    NzIB = 6;
elseif NzIB~=6 && NzIB~=30
    error('Unexpected number of IB zones: %d (expected 6 or 30).', NzIB);
end

% optional per-zone scale
if isempty(cfg.scale)
    cfg.scale = ones(1, NzIB);
else
    cfg.scale = cfg.scale(:).';
    if numel(cfg.scale) ~= NzIB
        error('scale must have length %d (zones).', NzIB);
    end
end

single = struct('t', t1, 'h', h1, 'Tg', Tg1K);

% ---------- 2) stitch repeated schedule ----------
bc = rifling_schedule_bc(single, cfg.Nshots, cfg.gap_s, ...
                         'Tamb_C', cfg.Tamb_C, 'zcols', NzIB, 'scale', cfg.scale);

% ---------- 3) quiet smoke test (optional) ----------
if cfg.smoke_test
    try
        test_bc = struct('t', bc.t(1:2), 'h', bc.h(1:2,:), 'Tg', bc.Tg(1:2));
        % Silence solver prints during the tiny probe run:
        evalc(['heat_transfer_2d_solver(''', 'Nr', ''',16, ''', 'Nz', ''',8, ', ...
               '''', 'dt_fd', ''',1e-3, ''', 'tEnd_fd', ''',1e-3, ''', 'theta', ''',1.0, ', ...
               '''', 'plot', ''',false, ''', 'debug', ''',false, ''', 'flux_check', ''',false, ', ...
               '''', 'bc', ''', test_bc);']);
    catch ME
        if contains(ME.message,'Unknown option "bc"') || contains(ME.message,'too many input arguments')
            error(['Your heat_transfer_2d_solver does not accept the ''bc'' option yet.', newline, ...
                   'Patch it (use the version you pasted with ''bc'') and rerun.']);
        else
            rethrow(ME);
        end
    end
end

% ---------- 4) quick BC plots ----------
if cfg.plot
    figure('Name','BC schedule','Color','w');
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile; plot(bc.t, bc.Tg - 273.15, 'k-','LineWidth',1.0); grid on;
    ylabel('T_g [°C]'); title(sprintf('Gas Temperature (N=%d, gap=%.3fs)', cfg.Nshots, cfg.gap_s));
    nexttile; 
    plot(bc.t, bc.h(:,1)/1e3, 'r-','LineWidth',1.0); hold on;
    idx = min(3,size(bc.h,2)); plot(bc.t, bc.h(:,idx)/1e3, 'b-','LineWidth',1.0);
    grid on; xlabel('t [s]'); ylabel('h [kW/m^2K]'); title('Sample HTCs (zones 1 & 3)');
    linkaxes(tl.Children,'x');
end
end

% ===== local helper =====
function bc = rifling_schedule_bc(single, Nshots, gap_s, varargin)
    p = inputParser;
    p.addParameter('Tamb_C', 20);
    p.addParameter('zcols', size(single.h,2));
    p.addParameter('scale',  ones(1, size(single.h,2)));
    p.parse(varargin{:});
    TambK = p.Results.Tamb_C + 273.15;

    t1   = single.t(:);
    h1   = single.h .* p.Results.scale(:).';
    Tg1K = single.Tg(:);
    Nz   = size(h1,2);

    % one block shot+gap
    t_shot   = t1;        h_shot   = h1;   Tg_shotK = Tg1K;
    if gap_s > 0
        t_gap   = [0; gap_s]; 
        h_gap   = zeros(2,Nz);
        Tg_gapK = ones(2,1) * TambK;
    else
        t_gap   = 0;   h_gap = zeros(1,Nz);   Tg_gapK = TambK;
    end
    block_t   = [t_shot;   t_shot(end) + t_gap];
    block_h   = [h_shot;   h_gap];
    block_TgK = [Tg_shotK; Tg_gapK];

    % tile N copies
    t = []; h = []; TgK = []; t0 = 0;
    for k = 1:Nshots
        t   = [t;   t0 + block_t]; %#ok<AGROW>
        h   = [h;   block_h];
        TgK = [TgK; block_TgK];
        t0  = t0 + (t1(end) + gap_s);
    end

    % dedup touching endpoints
    [t, idx] = unique(t,'stable');
    h   = h(idx,:);
    TgK = TgK(idx);

    bc = struct('t', t, 'h', h, 'Tg', TgK, 'Ttot', t(end), ...
                'meta', struct('Nshots', Nshots, 'gap_s', gap_s, 'Tamb_C', p.Results.Tamb_C));
end
