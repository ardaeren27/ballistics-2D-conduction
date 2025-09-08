% Parameter definitions, should accept any number of arguments
% Declares geometry, body scripts may also need tweaking

function p = parameters(varargin)

% Usage examples:
%   p = parameters();                                   
%   p = parameters('steel','38HMJ','tCr_um',150);       % Ex: 150 µm chromium
% All units SI unless stated. Temperatures in °C for material property fns.

    % ---------- parse name-value inputs ----------
    opts = struct('steel','DUPLEX','tCr_m',0,'check',true);
    if mod(numel(varargin),2)~=0
        error('parameters: name-value pairs expected');
    end
    for i=1:2:numel(varargin)
        key = lower(varargin{i});
        val = varargin{i+1};
        switch key
            case 'steel'
                opts.steel = upper(string(val));
            case 'tcr_m'
                opts.tCr_m = double(val);
            case 'tcr_um'
                opts.tCr_m = double(val)*1e-6;
            case 'check'
                opts.check = logical(val);
            otherwise
                error('Unknown option "%s".', key);
        end
    end
    if opts.tCr_m < 0
        error('Chromium thickness must be >= 0.');
    end

    % ---------- geometry (from article's piecewise outer radius) ----------
    p.geom.Rin_bore_m    = 17.5e-3;  % bore radius (to gas)
    p.geom.z_break_m     = [0, 0.385, 0.535, 0.880, 2.081, 2.980, 3.150];
    p.geom.Rout_break_m  = 1e-3 * [55.00, 55.00, 57.00, 59.50, 44.07, 31.00, 31.00];
    p.geom.L_m           = p.geom.z_break_m(end);

    % chromium lining on the bore
    p.geom.tCr_m         = opts.tCr_m;                       % thickness knob
    p.geom.Rin_cr_in_m   = p.geom.Rin_bore_m;                % gas/Cr interface
    p.geom.Rout_cr_out_m = p.geom.Rin_bore_m + p.geom.tCr_m; % Cr/steel interface
    p.geom.Rin_steel_m   = p.geom.Rout_cr_out_m;             % steel inner radius

    % piecewise-linear outer radius r_out(z) [m]
    p.geom.Rout          = @Rout_piecewise;

    % convenient helpers
    p.geom.z_mid_m         = 0.5*(p.geom.z_break_m(1:end-1) + p.geom.z_break_m(2:end));
    p.geom.thickness_total = @(z) p.geom.Rout(z) - p.geom.Rin_bore_m;
    p.geom.thickness_cr    = @(z) 0.*z + p.geom.tCr_m;   % constant in z
    p.geom.thickness_steel = @(z) p.geom.Rout(z) - p.geom.Rin_steel_m;

    % sanity checks (optional)
    if opts.check
        tsteel_edges = p.geom.thickness_steel(p.geom.z_break_m(1:end-1));
        if any(tsteel_edges <= 0)
            error('Chromium thickness (%.0f µm) leaves no steel in some zones.', 1e6*p.geom.tCr_m);
        end
    end

    % ---------- materials ----------
    [rhoS, cpS, kS] = load_steel(opts.steel);     % ρ [kg/m3], cp [J/kgK], k [W/mK]
    [rhoC, cpC, kC] = load_chromium();

    p.materials.steel.tag = char(opts.steel);
    p.materials.steel.rho = rhoS;
    p.materials.steel.cp  = cpS;
    p.materials.steel.k   = kS;

    p.materials.chromium.rho = rhoC;
    p.materials.chromium.cp  = cpC;
    p.materials.chromium.k   = kC;

    % ---------- handy thermal-resistance helpers (per unit length) ----------
    % R'cyl = ln(r2/r1) / (2*pi*k)  (supports vectorized r2 via z)
    p.thermal.Rcyl_per_m = @(k, r1, r2) log(r2./r1) ./ (2*pi*k);
    p.thermal.Rlayers_per_m = @(Tmean, z) ...
        p.thermal.Rcyl_per_m(p.materials.chromium.k(Tmean), p.geom.Rin_bore_m,  p.geom.Rout_cr_out_m) + ...
        p.thermal.Rcyl_per_m(p.materials.steel.k(Tmean),    p.geom.Rin_steel_m, p.geom.Rout(z));

    % ---------- nested/local functions ----------
    function Rout_val = Rout_piecewise(z)
        % Piecewise-linear outer radius r_out(z) [m]; preserves input shape
        Rout_val = interp1( ...
            p.geom.z_break_m, ...
            p.geom.Rout_break_m, ...
            z, 'linear', 'extrap');
    end

end

% ===== material property libraries (local to file) ========================
function [rho_fun, cp_fun, k_fun] = load_steel(tag)
switch upper(tag)
    case '30HN2MFA'
        T_k=[54.2,149.1,250.0,352.0,453.3,553.6,651.1,704.4,723.0,743.3,763.0,783.1,802.8,822.9,842.8,904.9,1004.3];
        k_vals=[35.9,37.3,36.0,33.8,30.9,27.2,19.7,17.1,16.0,15.8,17.1,18.7,19.3,19.5,19.7,20.3,20.6];
        T_cp=[38,70,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,991];
        cp_vals=[0.440,0.462,0.475,0.492,0.505,0.517,0.528,0.539,0.550,0.560,0.569,0.579,0.589,0.598,0.607,0.616,0.625,0.634,0.642,0.658];
        T_rho=[50,100,200,400,600,700,720,725,730,735,740,750,765,770,775,785,800,900,1060];
        rho_vals=[7.77,7.75,7.72,7.65,7.59,7.55,7.55,7.55,7.55,7.55,7.55,7.57,7.58,7.59,7.59,7.59,7.58,7.53,7.46];
    case 'DUPLEX'
        T_k=[52.0,149.1,249.8,351.7,457.0,553.6,654.7,704.2,744.2,762.9,782.6,802.5,811.7,821.7,842.3,904.7,1004.1];
        k_vals=[13.3,15.3,17.0,17.8,18.1,18.7,20.0,20.8,21.4,21.6,21.8,22.2,22.3,22.4,22.7,23.7,25.9];
        T_cp=[38,70,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,991];
        cp_vals=[0.417,0.442,0.462,0.492,0.515,0.534,0.548,0.559,0.567,0.572,0.576,0.579,0.582,0.584,0.588,0.594,0.602,0.614,0.629,0.668];
        T_rho=[50,100,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1060];
        rho_vals=[7.74,7.72,7.69,7.67,7.65,7.63,7.61,7.60,7.58,7.56,7.54,7.52,7.49,7.47,7.45,7.43,7.40,7.34,7.32];
    case '38HMJ'
        T_k=[50.9,149.0,250.0,351.3,453.4,553.6,654.7,704.5,741.0,762.6,782.7,802.8,811.9,821.8,842.3,904.7,1004.2];
        k_vals=[30.0,33.6,34.4,33.0,30.7,27.4,22.5,19.4,16.4,19.3,20.9,23.2,24.4,25.4,26.3,27.8,28.9];
        T_cp=[38,70,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,991];
        cp_vals=[0.458,0.485,0.502,0.525,0.543,0.559,0.574,0.587,0.598,0.609,0.618,0.626,0.634,0.640,0.645,0.650,0.653,0.656,0.657,0.658];
        T_rho=[50,100,200,400,600,780,795,800,820,830,840,850,860,870,880,890,900,1000,1060];
        rho_vals=[7.66,7.65,7.62,7.55,7.48,7.42,7.41,7.42,7.43,7.43,7.43,7.43,7.43,7.42,7.42,7.42,7.41,7.37,7.34];
    otherwise
        error('Unknown steel type "%s"', tag);
end
k_fun  = @(Tc) interp1(T_k,  k_vals,  Tc, 'linear','extrap');
cp_fun = @(Tc) 1e3*interp1(T_cp, cp_vals, Tc, 'linear','extrap');  % kJ/kgK -> J/kgK
rho_fun= @(Tc) 1e3*interp1(T_rho, rho_vals, Tc, 'linear','extrap');% g/cm3 -> kg/m3
end

function [rho_fun, cp_fun, k_fun] = load_chromium()
T_ref  = [20 100 200 300 400 500 600 700 800 900 1000];
k_ref  = [94  90  85  80  75  70  66  62  59  56  53];
cp_ref = [450 480 520 560 590 610 630 650 670 690 710];
rho_ref= [7200 7195 7185 7170 7155 7140 7125 7110 7095 7080 7065];
k_fun  = @(Tc) interp1(T_ref,  k_ref,  Tc, 'linear','extrap');
cp_fun = @(Tc) interp1(T_ref,  cp_ref, Tc, 'linear','extrap');
rho_fun= @(Tc) interp1(T_ref,  rho_ref,Tc, 'linear','extrap');
end
