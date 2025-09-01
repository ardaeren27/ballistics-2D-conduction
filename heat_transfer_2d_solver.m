function out = heat_transfer_2d_solver(varargin)
% HEAT_TRANSFER_2D_SOLVER  Axisymmetric (r,z) transient heat conduction (FV, fully-implicit).
% Self-contained; mirrors the corrected 1-D solver's cylindrical FV & Robin (combined resistance).
%
% Uses interiorBallistics.m (timeVec [s], hHist [Nt x 6 or 30], Tgas [K]).
%
% Name–value options (matches your call):
%   'steel'         : 'DUPLEX' (default). (Other values fallback to DUPLEX tables)
%   'tCr_m'         : chromium thickness [m] (default 0)
%   'tCr_um'        : chromium thickness [µm] (overrides tCr_m)
%   'Nr','Nz'       : grid sizes in r,z
%   'dt_fd'         : fixed thermal time step (s)
%   'tEnd_fd'       : end time (s)
%   'theta'         : theta-scheme (default 1.0 fully implicit)
%   'plot'          : true/false
%   'debug'         : true/false
%   'debug_stride'  : print frequency
%   'keep_stride'   : snapshot stride (default 50)
%
% Example:
% out = heat_transfer_2d_solver('steel','DUPLEX','tCr_um',0, ...
%     'Nr',240*4,'Nz',120,'dt_fd',1e-4,'tEnd_fd',20e-3, ...
%     'debug',true,'debug_stride',10);

% -------------------- options --------------------
opt = struct('steel','DUPLEX','tCr_m',0,'tCr_um',[], ...
             'Nr',240,'Nz',120,'theta',1.0, ...
             'dt_fd',1e-5,'tEnd_fd',8e-3, ...
             'plot',true,'debug',false,'debug_stride',50, ...
             'keep_stride',50,'flux_check',true);
if mod(numel(varargin),2)~=0, error('Use name-value pairs.'); end
for k=1:2:numel(varargin)
    key = lower(varargin{k});
    val = varargin{k+1};
    switch key
        case 'steel',         opt.steel = upper(string(val));
        case 'tcr_m',         opt.tCr_m = double(val);
        case 'tcr_um',        opt.tCr_m = double(val)*1e-6;
        case 'nr',            opt.Nr    = double(val);
        case 'nz',            opt.Nz    = double(val);
        case 'theta',         opt.theta = double(val);
        case 'dt_fd',         opt.dt_fd = double(val);
        case 'tend_fd',       opt.tEnd_fd = double(val);
        case 'plot',          opt.plot = logical(val);
        case 'debug',         opt.debug = logical(val);
        case 'debug_stride',  opt.debug_stride = double(val);
        case 'keep_stride',   opt.keep_stride = double(val);
        case 'flux_check',    opt.flux_check  = logical(val);
        otherwise, error('Unknown option "%s".', key);
    end
end
thetaFD = opt.theta;
if opt.tCr_m < 0, error('tCr must be >= 0'); end

% -------------------- geometry (matches 1-D corrected sections) --------------------
Rin = 0.0175;                 % inner bore radius [m]
T0  = 20.0;                   % ambient/initial [°C]
hout = 9.2;                   % outer HTC [W/m^2-K]

sections.z_mid = [0.2649, 0.4605, 0.7075, 1.4805, 2.5305, 3.0655]; % [m]
sections.R_out = [0.055, 0.0555, 0.0575, 0.052, 0.044, 0.031];     % [m]
sections.names = {'S1','S3','S4','S16','S25','S30'};

% Barrel length: ensure it spans last section; add margin
L = max(sections.z_mid) + 0.20;
% Piecewise-linear outer radius along z
z_breaks = sections.z_mid(:);
Rout_breaks = sections.R_out(:);
Rout_fun = @(zz) interp1(z_breaks, Rout_breaks, zz, 'linear','extrap');

% Chromium (optional)
Rcr = Rin + opt.tCr_m;  % chromium outer radius

% -------------------- grids --------------------
z  = linspace(0, L, opt.Nz);
dz = z(2)-z(1);
Rout_z = Rout_fun(z);
tmax = max(Rout_z) - Rin;
r_max = Rin + tmax;
r  = linspace(Rin, r_max, opt.Nr);       % cell centers
% faces
r_face = zeros(opt.Nr+1,1); r_face(1)=Rin; r_face(end)=r_max;
for i=2:opt.Nr, r_face(i)=0.5*(r(i-1)+r(i)); end

% validity mask truncated by Rout(z)
valid = false(opt.Nr,opt.Nz);
i_end = zeros(1,opt.Nz);
for j=1:opt.Nz
    Ni = find(r_face(2:end) <= Rout_z(j)+1e-15, 1, 'last');
    if isempty(Ni), Ni=0; end
    valid(1:Ni,j) = true;  i_end(j) = Ni;
end
if any(i_end<1)
    warning('Some z-columns have <1 radial cell. Increase Nr or check geometry.');
end

% FV metrics (per valid cell)
V    = zeros(opt.Nr,opt.Nz);   % volume
A_rW = zeros(opt.Nr,opt.Nz);   % inner radial face area
A_rE = zeros(opt.Nr,opt.Nz);   % outer radial face area
A_zS = zeros(opt.Nr,opt.Nz);   % south axial face area
A_zN = zeros(opt.Nr,opt.Nz);   % north axial face area
for j=1:opt.Nz
    for i=1:i_end(j)
        rW = r_face(i);  rE = r_face(i+1);
        V(i,j)   = pi*(rE^2 - rW^2)*dz;
        A_rW(i,j)= 2*pi*rW*dz;
        A_rE(i,j)= 2*pi*rE*dz;
        A_zS(i,j)= pi*(rE^2 - rW^2);
        A_zN(i,j)= A_zS(i,j);
    end
end

% -------------------- ballistics (same file) --------------------
fprintf('Loading ballistics data...\n');
if exist('interiorBallistics.m','file')~=2
    error('interiorBallistics.m not found.');
end
run('interiorBallistics.m');  % creates: timeVec, hHist, Tgas (K)
if ~exist('timeVec','var') || ~exist('hHist','var') || ~exist('Tgas','var')
    error('Ballistics script must define timeVec, hHist, Tgas.');
end
t_ib = timeVec(:);
h_data = hHist;
Tgas = Tgas(:);
if median(Tgas)>200, Tgas = Tgas - 273.15; end  % °C

% Force to 6 columns if 30 provided
if size(h_data,2) ~= 6
    if size(h_data,2)==30
        warning('hHist has 30 cols; using [1,5,10,16,25,30] for the 6 sections.');
        h_data = h_data(:,[1,5,10,16,25,30]);
    else
        error('Expected 6 HTC columns, got %d.', size(h_data,2));
    end
end

min_len = min([numel(t_ib), size(h_data,1), numel(Tgas)]);
t_ib = t_ib(1:min_len);
h_data = h_data(1:min_len,:);
Tgas   = Tgas(1:min_len);

% Thermal time grid
t = (0:opt.dt_fd:opt.tEnd_fd).';
Nt = numel(t);

% Map each z-column to nearest section (like 1-D corrected)
zone_id = zeros(1,opt.Nz);
for j=1:opt.Nz
    [~, idx] = min(abs(z(j) - sections.z_mid));
    zone_id(j) = idx;  % 1..6
end

% Pick a representative column per section for reporting
jstar = zeros(1,6);
for k2=1:6
    [~, jstar(k2)] = min(abs(z - sections.z_mid(k2)));
end
fprintf('Δz to section mids [mm]: %s\n', sprintf('%.1f ', 1e3*abs(z(jstar)-sections.z_mid)));

% -------------------- materials (DUPLEX + Chromium) --------------------
T_props   = [20 100 200 300 400 500 600 700 800 900 1000]; % °C
k_steel_D = [13.3 15.3 17.0 17.8 18.1 18.7 20.0 20.8 21.4 23.7 25.9];
cp_steel_D= [417  462  515  548  567  576  582  588  602  629  668 ];
rho_steel_D=[7740 7720 7690 7650 7610 7580 7540 7490 7450 7400 7340];

% Chromium (coating)
k_chrom   = [94 90 85 80 75 70 66 62 59 56 53];
cp_chrom  = [450 480 520 560 590 610 630 650 670 690 710];
rho_chrom = [7200 7195 7185 7170 7155 7140 7125 7110 7095 7080 7065];

% choose steel set (fallback to DUPLEX tables)
switch upper(string(opt.steel))
    case "DUPLEX"
        kS_tab = k_steel_D; cpS_tab = cp_steel_D; rhoS_tab = rho_steel_D;
    otherwise
        warning('steel="%s" not recognized; using DUPLEX properties.', opt.steel);
        kS_tab = k_steel_D; cpS_tab = cp_steel_D; rhoS_tab = rho_steel_D;
end
kS_fun   = @(T) interp1(T_props, kS_tab,   T, 'linear','extrap');
cpS_fun  = @(T) interp1(T_props, cpS_tab,  T, 'linear','extrap');
rhoS_fun = @(T) interp1(T_props, rhoS_tab, T, 'linear','extrap');

kC_fun   = @(T) interp1(T_props, k_chrom,   T, 'linear','extrap');
cpC_fun  = @(T) interp1(T_props, cp_chrom,  T, 'linear','extrap');
rhoC_fun = @(T) interp1(T_props, rho_chrom, T, 'linear','extrap');

% -------------------- storage --------------------
T = T0*ones(opt.Nr,opt.Nz);
keep_idx = unique([1:opt.keep_stride:Nt, Nt]);
Ns = numel(keep_idx);
T_snap = nan(opt.Nr,opt.Nz,Ns);
snap_k = 0;
T_inner6 = T0*ones(Nt,6);

% first-on flux consistency flags
first_on_done    = false(1,6);
first_on_pending = false(1,6);

fprintf('FV 2-D: Nr=%d Nz=%d Nt=%d theta=%.2f dt=%.3g s  L=%.2f m  tCr=%.0f µm  steel=%s\n', ...
    opt.Nr, opt.Nz, Nt, thetaFD, opt.dt_fd, L, 1e6*opt.tCr_m, upper(string(opt.steel)));

% -------------------- time-march --------------------
for n = 1:Nt-1
    tN   = t(n);
    tNp1 = t(n+1);
    dt   = tNp1 - tN;

    % interpolate ballistics at t^{n+1}
    Tg   = interp1(t_ib, Tgas,  tNp1, 'linear','extrap');  % °C (scalar)
    h6   = interp1(t_ib, h_data, tNp1, 'linear','extrap'); % 1x6
    if numel(h6)==1, h6 = repmat(h6,1,6); end
    h_in_col = h6(zone_id); % 1xNz

    % mark first-on events
    for k2=1:6
        if ~first_on_done(k2) && (h6(k2) > 0)
            first_on_pending(k2) = true;
        end
    end

    % lagged properties at T^n
    k_n   = zeros(opt.Nr,opt.Nz);
    rhoCp = zeros(opt.Nr,opt.Nz);
    for j=1:opt.Nz
        for i=1:i_end(j)
            Tc = T(i,j);
            if r(i) <= Rcr + 1e-12 && opt.tCr_m>0
                k_n(i,j)   = kC_fun(Tc);
                rhoCp(i,j) = rhoC_fun(Tc)*cpC_fun(Tc);
            else
                k_n(i,j)   = kS_fun(Tc);
                rhoCp(i,j) = rhoS_fun(Tc)*cpS_fun(Tc);
            end
        end
    end

    % unknown numbering
    map = zeros(opt.Nr,opt.Nz);  map(valid) = 1:nnz(valid);
    Nunk = nnz(valid);

    rows = zeros(7*Nunk,1); cols = rows; vals = rows;
    rhs  = zeros(Nunk,1);   nzctr = 0;

    % assemble
    for j=1:opt.Nz
        Ni = i_end(j); if Ni==0, continue; end
        for i=1:Ni
            P = map(i,j); if P==0, continue; end

            % transient
            aP_time = rhoCp(i,j) * V(i,j) / dt;
            aP = aP_time;  bP = aP_time * T(i,j);

            % axial conduction (zero-flux at ends by omission)
            if j>1 && valid(i,j-1)
                S = map(i,j-1);
                kS = face_k(k_n(i,j-1), k_n(i,j));
                condS = thetaFD * (kS * A_zS(i,j) / dz);
                [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,S,-condS);
                aP = aP + condS;
            end
            if j<opt.Nz && valid(i,j+1)
                Nn = map(i,j+1);
                kNf = face_k(k_n(i,j+1), k_n(i,j));
                condN = thetaFD * (kNf * A_zN(i,j) / dz);
                [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,Nn,-condN);
                aP = aP + condN;
            end

            % radial faces
            if i==1
                % inner Robin at Rin using combined resistance
                drw = max(r(i) - Rin, 1e-12);
                hin = max(0, h_in_col(j));
                if hin > 0
                    Uin = 1 / (drw/max(k_n(i,j),eps) + 1/hin); % W/m^2K
                    aP = aP + thetaFD*(Uin * A_rW(i,j));
                    bP = bP + thetaFD*(Uin * A_rW(i,j))*Tg;
                end
                % east conduction to i=2
                if Ni>=2
                    E = map(2,j);
                    kE = face_k(k_n(2,j), k_n(1,j));
                    dE = max(r(2)-r(1), 1e-12);
                    condE = thetaFD * (kE * A_rE(1,j) / dE);
                    [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,E,-condE);
                    aP = aP + condE;
                end

            elseif i==Ni
                % outer Robin at Rout(z_j)
                dro = max(Rout_z(j) - r(i), 1e-12);
                Uout = 1 / (dro/max(k_n(i,j),eps) + 1/hout);
                Aout = 2*pi*Rout_z(j)*dz;
                aP = aP + thetaFD*(Uout * Aout);
                bP = bP + thetaFD*(Uout * Aout)*T0;

                % west conduction
                W = map(i-1,j);
                kW = face_k(k_n(i-1,j), k_n(i,j));
                dW = max(r(i)-r(i-1), 1e-12);
                condW = thetaFD * (kW * A_rW(i,j) / dW);
                [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,W,-condW);
                aP = aP + condW;

            else
                % interior radial conduction
                W = map(i-1,j);
                kW = face_k(k_n(i-1,j), k_n(i,j));
                dW = max(r(i)-r(i-1), 1e-12);
                condW = thetaFD * (kW * A_rW(i,j) / dW);
                [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,W,-condW);
                aP = aP + condW;

                E = map(i+1,j);
                kE = face_k(k_n(i+1,j), k_n(i,j));
                dE = max(r(i+1)-r(i), 1e-12);
                condE = thetaFD * (kE * A_rE(i,j) / dE);
                [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,E,-condE);
                aP = aP + condE;
            end

            % finalize
            [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,P,P,aP);
            rhs(P) = bP;
        end
    end

    % solve
    A = sparse(rows(1:nzctr), cols(1:nzctr), vals(1:nzctr), Nunk, Nunk);
    sol = A \ rhs;
    T(valid) = sol;

    % record inner wall at 6 sections
    for k2=1:6
        jj = jstar(k2);
        if i_end(jj)>=1, T_inner6(n+1,k2) = T(1,jj); else, T_inner6(n+1,k2) = NaN; end
    end

    % first-on flux consistency (diagnostic)
    if opt.flux_check
        for k2=1:6
            if first_on_pending(k2)
                jj = jstar(k2); if i_end(jj)<1, continue; end
                Tw = T(1,jj);  knw = k_n(1,jj);
                drw = max(r(1)-Rin, 1e-12);
                hin_now = max(0, h6(k2));
                qR = hin_now*(Tw - Tg);
                if hin_now>0
                    Uin = 1/(drw/max(knw,eps) + 1/hin_now);
                    qFV = Uin*(Tw - Tg);
                else
                    qFV = 0;
                end
                ratio = qFV / max(1e-16, qR);
                fprintf('S%d first-on: h=%0.2e  Tw=%.1f  Tg=%.1f  qR=%+.2e  qFV=%+.2e  ratio=%.3g\n', ...
                    k2, hin_now, Tw, Tg, qR, qFV, ratio);
                first_on_done(k2)=true; first_on_pending(k2)=false;
            end
        end
    end

    % debug print
    if opt.debug && (mod(n,opt.debug_stride)==0 || n==1 || n==Nt-1)
        Tin = T_inner6(n+1,:);
        Tmask = T;  Tmask(~valid)= NaN;
        Tmax = max(Tmask,[],'all','omitnan');
        Tmin = min(Tmask,[],'all','omitnan');
        fprintf('%5d  t=%.3f ms  dt=%.1f µs | Tin: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f | Tmin/Tmax: %6.2f / %6.2f\n', ...
            n, 1e3*tNp1, 1e6*dt, Tin, Tmin, Tmax);
    end

    % snapshots
    if any(n==keep_idx)
        snap_k = snap_k + 1;
        Tf = T; Tf(~valid)=NaN;
        T_snap(:,:,snap_k) = Tf;
    end
end

% final line for inner temperature
T_inner6(Nt,:) = T_inner6(Nt-1,:);

% -------------------- plots --------------------
if opt.plot
    figure('Name','Inner-surface T (2-D FV)','Color','w');
    hold on; grid on;
    plot(1e3*t, T_inner6(:,1),'LineWidth',1.4);
    plot(1e3*t, T_inner6(:,2),'LineWidth',1.4);
    plot(1e3*t, T_inner6(:,3),'LineWidth',1.4);
    plot(1e3*t, T_inner6(:,4),'LineWidth',1.4);
    plot(1e3*t, T_inner6(:,5),'LineWidth',1.4);
    plot(1e3*t, T_inner6(:,6),'LineWidth',1.4);
    xlabel('time [ms]'); ylabel('Inner-surface T [°C]');
    title('Inner-surface temperature at sections (2-D FV, DUPLEX)');
    legend(sections.names, 'Location','best');

    if snap_k>0
        Tf = T_snap(:,:,snap_k);
        figure('Name','Temperature field (final snapshot)','Color','w');
        imagesc(z, r*1e3, Tf); set(gca,'YDir','normal'); colorbar
        xlabel('z [m]'); ylabel('r [mm]');
        title('T field [°C] (final snapshot)');
        hold on; plot(z, Rout_z*1e3, 'k--','LineWidth',1.0); hold off
    end

    % BC check
    figure('Name','Boundary conditions','Color','w');
    subplot(2,1,1);
    plot(1e3*t_ib, Tgas, 'r-','LineWidth',1.5); grid on;
    xlabel('time [ms]'); ylabel('T_g [°C]'); title('Gas temperature');
    subplot(2,1,2);
    plot(1e3*t_ib, h_data(:,1)/1000, 'k-', 'LineWidth',1.2); hold on;
    plot(1e3*t_ib, h_data(:,3)/1000, 'b-', 'LineWidth',1.2);
    plot(1e3*t_ib, h_data(:,6)/1000, 'c-', 'LineWidth',1.2);
    xlabel('time [ms]'); ylabel('h [kW/m^2K]'); grid on;
    legend('Sec1','Sec3','Sec6','Location','best');
    title('Sample HTCs');
end

% -------------------- outputs --------------------
out = struct('t',t,'z',z,'r',r,'valid',valid, ...
             'T_snap',T_snap(:,:,1:max(1,snap_k)), ...
             'T_inner6',T_inner6, ...
             'Rout_z',Rout_z, ...
             'sections',sections, ...
             'opts',opt);
end

% ===== helpers =====
function kf = face_k(kL,kR)
    kf = 2*kL.*kR ./ max(kL + kR, eps);  % harmonic mean
end
function [rows,cols,vals,nzctr] = acc(rows,cols,vals,nzctr,i,j,v)
    nzctr = nzctr + 1; rows(nzctr)=i; cols(nzctr)=j; vals(nzctr)=v;
end
