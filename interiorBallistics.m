%% Interior-ballistics two-phase solver (interior + transitional)
% Arda Eren

tic

% ────────────────────────────────────────────────────────────────────────
%% CONSTANTS (main script scope)
% ────────────────────────────────────────────────────────────────────────
mP    = 0.380;        % projectile mass only [kg]
mp    = 0.376;        % propellant mass [kg]
s     = 9.98e-4;      % bore cross-sectional area [m^2]
V0    = 373e-6;       % initial free volume [m^3]
lm    = 2.934;        % barrel length to muzzle from breech [m]
K     = 1.37;         % gas inertia factor (phi = K, dimensionless, constant)
f     = 1.071e6;      % propellant energy parameter [J/kg]
eta   = 1.064e-3;     % covolume parameter [-]
gamma = 1.2;          % specific heat ratio [-]
R     = 340;          % specific gas constant [J/(kg·K)]
rho_p = 1600;         % propellant density [kg/m^3]

% Shot-start (release) threshold
p_start   = 30e6;     % 30 MPa, p_0 in the article
released  = false;    % latched once p >= p_start
t_release = NaN;      % for debugging/plotting if you want

% Geometry per the article (absolute positions from breech start)
l0 = 0.216;                                   % chamber length to PS1 [m]
li = [l0 0.385 0.535 0.880 2.081 2.980];      % z-locations P1..P6 [m]
Di = 0.035*ones(1,6);                         % bore diameter per section [m]

% ────────────────────────────────────────────────────────────────────────
% TIME SETTINGS
% ────────────────────────────────────────────────────────────────────────
dt      = 1e-6;        % keep small; the release is stiff
tEnd    = 0.1;         % analysis end time, if no break condition is triggered
maxStep = ceil(tEnd/dt)+1;

% ────────────────────────────────────────────────────────────────────────
% P R E - A L L O C A T I O N
% ────────────────────────────────────────────────────────────────────────
timeVec     = zeros(maxStep,1);
location    = zeros(maxStep,1);
velocityP   = zeros(maxStep,1);   velocityP(1) = 0;   % sanity IC
zp          = zeros(maxStep,1);   zp(1) = 0.001;      % article IC
zetaHist    = zeros(maxStep,1);
theta       = zeros(maxStep,1);   theta(1) = R*3150;  % article IC
Tgas        = zeros(maxStep,1);   Tgas(1)  = 3150;    % article IC
rhoHist     = zeros(maxStep,1);
pressureEqn = zeros(maxStep,1);

% Convection coefficients per section, h(t) [W/m^2.K]
hHist = zeros(maxStep,6);

% helper state vector for RK4 during phase 1: [x  v  zp  θ]
state = [0 0 zp(1) theta(1)];
t_on  = NaN(1,6);   % first-on times for P1..P6 (debug purposes)

% For marking when z_p hits 1.0 in Phase-2
t_zp1 = NaN;
%% Phase 1 - Interior
% ────────────────────────────────────────────────────────────────────────
% ─────────────────────────   P H A S E  1   ─────────────────────────────
% ───────────────────────── (0 ≤ t ≤ t_exit) ─────────────────────────────
% ────────────────────────────────────────────────────────────────────────
n = 1;  t = 0;
t_exit = NaN;    % determined when base (tail) clears muzzle (x >= lm)
% run until base crosses the muzzle or loop-break trigger

%% 4th order Runge-Kutta Loop 
while t < tEnd
    % store previous state/time for crossing detection
    prev_state = state;    % [x v zp θ] at time t
    prev_t     = t;

    % ===== k1 =====
    [~,P1,dx1,dv1,dzp1,dth1] = thermoModel(prev_state(3), prev_state(4)/R, prev_state(2), prev_state(1));
    if ~released && P1 < p_start
        % hold projectile; no inertial work in energy eqn
        dx1 = 0; dv1 = 0;
        T1  = prev_state(4)/R;  z1 = max(prev_state(3), eps);
        dth1 = ((f - R*T1)*mp*dzp1) / (mp*z1);
    end
    k1 = dt*[dx1 dv1 dzp1 dth1];

    % ===== k2 =====
    mid = prev_state + 0.5*k1;
    mid(3) = min(max(mid(3),0),1);       % clamp zp at stage point
    [~,P2,dx2,dv2,dzp2,dth2] = thermoModel(mid(3), mid(4)/R, mid(2), mid(1));
    if ~released && P2 < p_start
        dx2 = 0; dv2 = 0;
        T2  = mid(4)/R;  z2 = max(mid(3), eps);
        dth2 = ((f - R*T2)*mp*dzp2) / (mp*z2);
    end
    k2 = dt*[dx2 dv2 dzp2 dth2];

    % ===== k3 =====
    mid = prev_state + 0.5*k2;
    mid(3) = min(max(mid(3),0),1);
    [~,P3,dx3,dv3,dzp3,dth3] = thermoModel(mid(3), mid(4)/R, mid(2), mid(1));
    if ~released && P3 < p_start
        dx3 = 0; dv3 = 0;
        T3  = mid(4)/R;  z3 = max(mid(3), eps);
        dth3 = ((f - R*T3)*mp*dzp3) / (mp*z3);
    end
    k3 = dt*[dx3 dv3 dzp3 dth3];

    % ===== k4 =====
    endS = prev_state + k3;
    endS(3) = min(max(endS(3),0),1);
    [~,P4,dx4,dv4,dzp4,dth4] = thermoModel(endS(3), endS(4)/R, endS(2), endS(1));
    if ~released && P4 < p_start
        dx4 = 0; dv4 = 0;
        T4  = endS(4)/R;  z4 = max(endS(3), eps);
        dth4 = ((f - R*T4)*mp*dzp4) / (mp*z4);
    end
    k4 = dt*[dx4 dv4 dzp4 dth4];

    % RK4 update
    new_state   = prev_state + (k1+2*k2+2*k3+k4)/6;
    new_state(3)= min(max(new_state(3),0),1);

    % Latch release if threshold reached this step 
    if ~released
        % pressure at step start and end (EOS)
        P_prev = (mp*prev_state(3)*R*(prev_state(4)/R)) / (V0 + s*prev_state(1) - mp*(1-prev_state(3))/rho_p - eta*mp*prev_state(3));
        P_new  = (mp*new_state(3)*R*(new_state(4)/R))   / (V0 + s*new_state(1)  - mp*(1-new_state(3))/rho_p  - eta*mp*new_state(3));
        if (P_prev < p_start) && (P_new >= p_start)
            % Optional: linearized release time within the step
            fracR     = (p_start - P_prev) / max(P_new - P_prev, eps);
            t_release = prev_t + max(0,min(1,fracR))*dt;
            released  = true;
        elseif P_new >= p_start
            released  = true;
            t_release = prev_t + dt;
        end
    end

    % detect if base crosses the muzzle within this step, section gating
    x_prev  = prev_state(1);
    x_new   = new_state(1);
    crossed = (x_prev < lm) && (x_new >= lm);

    if crossed
        % linearly interpolate within the step to the exact crossing
        frac       = (lm - x_prev) / max(x_new - x_prev, eps);   % ∈ (0,1]
        exit_state = prev_state + frac*(new_state - prev_state);
        exit_state(3) = min(max(exit_state(3),0),1);             % clamp zp at exit
        t     = prev_t + frac*dt;
        state = exit_state;

        % --- bookkeep at the exact exit moment ---
        location(n)    = state(1);
        velocityP(n)   = state(2);
        zp(n)          = state(3);
        theta(n)       = state(4);
        Tgas(n)        = theta(n)/R;

        pressureEqn(n) = (mp*zp(n)*R*Tgas(n)) / (V0 + s*location(n) - mp*(1-zp(n))/rho_p - eta*mp*zp(n));
        rhoHist(n)     = (mp*zp(n)) / (V0 + s*location(n) - mp*(1 - zp(n))/rho_p - eta*mp*zp(n));

        % HTC gating
        v_now    = velocityP(n);
        l_now    = location(n);
        tail_abs = l0 + max(l_now,0);
        wi       = (li./(l0 + l_now)) * v_now;
        h_now    = (6.1 ./ (Di.^0.2)) .* (rhoHist(n) * abs(wi)).^0.8;

        present  = (li <= tail_abs);
        h_now(~present) = 0;
        just_on  = present & isnan(t_on);  t_on(just_on) = t;
        hHist(n,:) = h_now;

        % finalize bookkeeping for exit moment
        timeVec(n) = t;
        t_exit     = t;
        n = n + 1;    % reserve next slot for Phase 2
        break
    else
        % no crossing; accept full step
        state = new_state;

        % --- bookkeep with the synced, updated state ---
        location(n)    = state(1);
        velocityP(n)   = state(2);
        zp(n)          = state(3);
        theta(n)       = state(4);
        Tgas(n)        = theta(n)/R;

        pressureEqn(n) = (mp*zp(n)*R*Tgas(n)) / (V0 + s*location(n) - mp*(1-zp(n))/rho_p - eta*mp*zp(n));
        rhoHist(n)     = (mp*zp(n)) / (V0 + s*location(n) - mp*(1 - zp(n))/rho_p - eta*mp*zp(n));

        % HTC gating
        v_now    = velocityP(n);
        l_now    = location(n);
        tail_abs = l0 + max(l_now,0);
        wi       = (li./(l0 + l_now)) * v_now;
        h_now    = (6.1 ./ (Di.^0.2)) .* (rhoHist(n) * abs(wi)).^0.8;

        present  = (li <= tail_abs);
        h_now(~present) = 0;
        just_on  = present & isnan(t_on);  t_on(just_on) = t;
        hHist(n,:) = h_now;

        % Debug prints
        fprintf('\nTime: %.8f s\n', t);
        fprintf('Location: %.6f m\n', location(n));
        fprintf('Velocity: %.6f m/s\n', velocityP(n));
        fprintf('Burn Ratio (zp): %.6f\n', zp(n));
        fprintf('Gas Temperature: %.2f K\n', Tgas(n));
        fprintf('Pressure: %.2f kPa\n', pressureEqn(n)/1000);

        % advance time & index
        t = prev_t + dt;
        timeVec(n) = t;
        n = n + 1;
    end
end

if isnan(t_exit) % Warning issue if a change breaks the physics
    warning('Projectile did not reach the muzzle before tEnd; forcing exit at tEnd.');
    t_exit = t;
end

% keep last piston position & velocity for plotting
location(n:end)  = lm;
velocityP(n:end) = 0;

%% Hand-off
% ────────────────────────────────────────────────────────────────────────
% H A N D O F F   C O N T I N U I T Y (seed Phase-2 state)
% ────────────────────────────────────────────────────────────────────────
zetaCur   = 0;                % article IC 
zp_cur    = zp(n-1);          % carry forward live z_p into Phase-2
T_phase2  = Tgas(n-1);        % Bulk-gas temp phase-2 seed

den2                = V0 + s*lm - mp*(1 - zp_cur)/rho_p - eta*mp*(zp_cur - zetaCur);
pressureEqn(n-1)    = (mp*(zp_cur - zetaCur)*R*T_phase2)/den2;
rhoHist(n-1)        = (mp*(zp_cur - zetaCur))/den2;
zetaHist(n-1)       = zetaCur;

% Seed h for Phase 2 using critical outflow speed (same mapping)
wcr            = sqrt(gamma*R*T_phase2) * sqrt(2/(gamma+1));
wi2            = (li./(l0 + lm)) * wcr;
hHist(n-1,:)   = (6.1 ./ (Di.^0.2)) .* (rhoHist(n-1) * abs(wi2)).^0.8;

%% Phase 2 - Transitional
% ────────────────────────────────────────────────────────────────────────
% ─────────────────────────  P H A S E   2  ──────────────────────────────
% ─────────────────────────  (t ≥ t_exit)  ───────────────────────────────
% ────────────────────────────────────────────────────────────────────────
% State vector during Phase-2: [zeta, theta, zp]
zeta        = zetaCur;
theta_here  = theta(n-1);
zp_here     = zp_cur;
%% 4th order Runge-Kutta Loop
while ((zp_here - zeta) > 1e-6) && t < tEnd
    % k1
    [~, ~, dzp1, dzeta1, dth1] = thermoModelv2(zp_here, theta_here/R, zeta);
    if zp_here >= 1-1e-12, dzp1 = 0; end
    dzp1 = max(dzp1,0);
    k1 = dt * [dzeta1, dth1, dzp1];

    % k2
    z_mid = min(max(zp_here + 0.5*k1(3),0),1);
    [~, ~, dzp2, dzeta2, dth2] = thermoModelv2(z_mid, (theta_here + 0.5*k1(2))/R, zeta + 0.5*k1(1));
    if z_mid >= 1-1e-12, dzp2 = 0; end
    dzp2 = max(dzp2,0);
    k2 = dt * [dzeta2, dth2, dzp2];

    % k3
    z_mid = min(max(zp_here + 0.5*k2(3),0),1);
    [~, ~, dzp3, dzeta3, dth3] = thermoModelv2(z_mid, (theta_here + 0.5*k2(2))/R, zeta + 0.5*k2(1));
    if z_mid >= 1-1e-12, dzp3 = 0; end
    dzp3 = max(dzp3,0);
    k3 = dt * [dzeta3, dth3, dzp3];

    % k4
    z_end = min(max(zp_here + k3(3),0),1);
    [~, ~, dzp4, dzeta4, dth4] = thermoModelv2(z_end, (theta_here + k3(2))/R, zeta + k3(1));
    if z_end >= 1-1e-12, dzp4 = 0; end
    dzp4 = max(dzp4,0);
    k4 = dt * [dzeta4, dth4, dzp4];

    % RK4 advance
    zeta       = zeta       + (k1(1)+2*k2(1)+2*k3(1)+k4(1))/6;
    theta_here = theta_here + (k1(2)+2*k2(2)+2*k3(2)+k4(2))/6;
    old_zp     = zp_here;
    zp_here    = zp_here + (k1(3)+2*k2(3)+2*k3(3)+k4(3))/6;
    zp_here    = min(max(zp_here,0),1);    % clamp

    % mark first time zp hits 1
    if isnan(t_zp1) && (old_zp < 1) && (zp_here >= 1)
        t_zp1 = t + dt;  % approx to end of step
    end

    % Bookkeeping at index n
    theta(n)    = theta_here;
    Tgas(n)     = theta_here/R;
    zp(n)       = zp_here;
    zetaHist(n) = zeta;

    den2             = V0 + s*lm - mp*(1 - zp_here)/rho_p - eta*mp*(zp_here - zeta);
    P2               = (mp*(zp_here - zeta)*R*Tgas(n))/den2;
    pressureEqn(n)   = max(P2, 101325);     % floor at 1 atm
    rhoHist(n)       = (mp*(zp_here - zeta))/den2;

    % Convection coefficient h_i(t) in Phase 2 (no gating; flow everywhere)
    wcr           = sqrt(gamma*R*Tgas(n)) * sqrt(2/(gamma+1));
    wi2           = (li./(l0 + lm)) * wcr;
    hHist(n,:)    = (6.1 ./ (Di.^0.2)) .* (rhoHist(n) * abs(wi2)).^0.8;

    % advance time/index
    t = t + dt;
    timeVec(n) = t;
    n = n + 1;

    % Debug prints (optional)
    fprintf('\nTime: %.8f s\n', t);
    fprintf('Zeta: %.6f\n', zetaHist(n-1));
    fprintf('Gas Temperature: %.2f K\n', Tgas(n-1));
    fprintf('Pressure: %.2f kPa\n', pressureEqn(n-1)/1000);

    % Early loop break as defined in article
    % May also be changed to atmospheric if needed
    if pressureEqn(n-1) <= 0.18e6
        break
    end
end
%% Memory handling
% truncate unused tail, convenient for later use
last        = min(n-1, numel(location));
timeVec     = timeVec(1:last);
location    = location(1:last);
velocityP   = velocityP(1:last);
zp          = zp(1:last);
zetaHist    = zetaHist(1:last);
Tgas        = Tgas(1:last);
pressureEqn = pressureEqn(1:last);
rhoHist     = rhoHist(1:last);
hHist       = hHist(1:last,:);

%% Rabbit + Plots
% Declaration of Success + Run-Time Indicator
% fprintf('Simulation finished in %.2f s (wall-clock).\n', toc);
fprintf(' {\\__/}\n');
fprintf('( • . •)\n');
fprintf('/ > 🥕>  Here are your plots\n');

% ────────────────────────────────────────────────────────────────────────
% P L O T S
% ────────────────────────────────────────────────────────────────────────
figure('Name','Interior Ballistics','Color','w');
tiledlayout(3,3,'TileSpacing','compact');

% Gas Temperature
nexttile
plot(timeVec, Tgas, 'LineWidth',1.2, "Color","black");
xlabel('t  [s]'); ylabel('T_{gas}  [K]'); title('Gas Temperature'); grid on

% Chamber Pressure
nexttile
yyaxis left
plot(timeVec, pressureEqn/1e6, 'LineWidth',1.2); ylabel('P [MPa] (linear)')
yyaxis right
semilogy(timeVec, pressureEqn/1e6, '--', 'LineWidth',1.0); ylabel('P [MPa] (log)')
xlabel('t [s]'); title('Chamber Pressure (linear + log)'); grid on
xline(t_exit,'k--','t_{exit}');
if ~isnan(t_zp1), xline(t_zp1,'r-.','z_p{=}1'); end
if ~isnan(t_release), xline(t_release,'g-.','release'); end

% Piston Velocity (FD of x for pre-exit visualization)
nexttile
idx = timeVec <= t_exit;
plot(timeVec(idx), [0; diff(location(idx))]/dt ,'LineWidth',1.2,"Color","red");
hold on; plot(timeVec(idx), velocityP(idx), 'LineWidth',1.2,"Color","red")
xlabel('t  [s]'); ylabel('v_P  [m/s]'); title('Piston Velocity'); grid on

% Propellant Burned
nexttile
plot(timeVec, zp, 'LineWidth',1.2,"Color","blue");
xlabel('t  [s]'); ylabel('z_p  [-]'); title('Propellant Burned'); grid on
if ~isnan(t_zp1), xline(t_zp1,'r-.','z_p{=}1'); end

% Gas Expelled
nexttile
plot(timeVec, zetaHist, 'LineWidth',1.2,"Color","magenta");
xlabel('t  [s]'); ylabel('\zeta  [-]'); title('Gas Expelled'); grid on

% Gas Density
nexttile
plot(timeVec, rhoHist, 'LineWidth',1.2,"Color","cyan");
xlabel('t  [s]'); ylabel('\rho  [kg/m^3]'); title('Gas Density'); grid on

% Convection coefficients + Tgas overlay
nexttile([1 3])
yyaxis left
plot(timeVec, hHist(:,1), 'k-',  'LineWidth', 1.1); hold on
plot(timeVec, hHist(:,2), 'r-',  'LineWidth', 1.1);
plot(timeVec, hHist(:,3), 'b-',  'LineWidth', 1.1);
plot(timeVec, hHist(:,4), 'g-',  'LineWidth', 1.1);
plot(timeVec, hHist(:,5), 'm-',  'LineWidth', 1.1);
plot(timeVec, hHist(:,6), 'y-',  'LineWidth', 1.1);
xlim([0, 0.007]);  ylim([0, 2.75e5]);
xlabel('t  [s]'); ylabel('h_i  [W/(m^2·K)]'); grid on
yyaxis right
plot(timeVec, Tgas, 'r--', 'LineWidth', 1.2);
ylim([0, 3.5e3]); ylabel('T_{gas}  [K]');
legend({'P1','P2','P3','P4','P5','P6','Tgas'}, 'Location','northwest');
title('Convection Coefficients & Gas Temperature'); grid on

% === Separate plot for convection coefficients + Tgas ===
figure;
tiledlayout(1,3);
nexttile([1 3]);
yyaxis left
plot(timeVec, hHist(:,1), 'k-',  'LineWidth', 1.1); hold on
plot(timeVec, hHist(:,2), 'r-',  'LineWidth', 1.1);
plot(timeVec, hHist(:,3), 'b-',  'LineWidth',  1.1);
plot(timeVec, hHist(:,4), 'g-',  'LineWidth',  1.1);
plot(timeVec, hHist(:,5), 'm-',  'LineWidth',  1.1);
plot(timeVec, hHist(:,6), 'y-',  'LineWidth',  1.1);
xlim([0, 0.007]);  ylim([0, 2.75e5]);
yticks([0 2.5e4 5e4 7.5e4 1e5 1.25e5 1.5e5 1.75e5 2e5 2.25e5 2.5e5 2.75e5]);
yticklabels({'0','2.5×10^4','5.0×10^4','7.5×10^4','1.0×10^5','1.25×10^5','1.50×10^5', ...
             '1.75×10^5','2.00×10^5','2.25×10^5','2.50×10^5','2.75×10^5'});
ylabel('Heat Transfer Coefficient HTC, W·m^{-2}·K^{-1}');
grid on
yyaxis right
plot(timeVec, Tgas, 'r--', 'LineWidth', 1.2);
ylim([0, 3.5e3]);
yticks([0 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500]);
yticklabels({'0','2.5×10^2','5.0×10^2','7.5×10^2','1.0×10^3','1.25×10^3','1.50×10^3', ...
             '1.75×10^3','2.00×10^3','2.25×10^3','2.50×10^3','2.75×10^3', ...
             '3.00×10^3','3.25×10^3','3.50×10^3'});
ylabel('Gas Temperature, K');
xlabel('Time, s');
title('Convection Coefficients & Gas Temperature');
legend({'P1','P2','P3','P4','P5','P6','Tgas'}, 'Location','northwest');
ax = gca; ax.YAxis(1).MinorTick = 'on'; ax.YAxis(2).MinorTick = 'on'; ax.XAxis.MinorTick = 'on';


%% ======================= Function definitions ============================
% Phase-1 Function Definition
function [rho, pressureEqn, diffLocation, diffVelocity, diffBurnRate, diffTheta] = thermoModel(zp, Tgas, velocityP, location)
    % Coefs and Params (self-contained)
    mP = 0.380; mp = 0.376; s = 9.98e-4; V0 = 373e-6; K = 1.37;  % phi = K (constant)
    f = 1.071e6; eta = 1.064e-3; gamma = 1.2; R = 340; rho_p = 1600;
    r1 = 0.597e-9; S1 = 134.4e-6; LAMBDA1 = 75.2e-9; Kappa1 = 0.755; lambda1 = 0.159;

    % Live chamber state
    pressureEqn = (mp*zp*R*Tgas) / (V0 + s*location - (mp/rho_p)*(1-zp) - eta*mp*zp);
    rho         = (mp*zp)        / (V0 + s*location - mp*(1-zp)/rho_p - eta*mp*zp);

    % Effective accelerating mass
    m_eff = mP;                           % moving lumped mass [kg]
    phi   = 1.37;                         % gas-inertia factor (dimensionless)

    % Kinematics and kinetics
    diffLocation = velocityP;
    diffVelocity = pressureEqn*s / (phi * m_eff);

    % Burning law (non-negative) and energy with consistent inertial work term
    diffBurnRate = (S1/LAMBDA1)*r1*pressureEqn*sqrt(1 + 4*zp*lambda1/Kappa1);
    diffBurnRate = max(diffBurnRate,0);
    if zp >= 1, diffBurnRate = 0; end    % freeze burn at full conversion

    % Full energy eqn (caller may override when locked)
    diffTheta = ((f - R*Tgas)*mp*diffBurnRate - (gamma - 1)*phi*m_eff*velocityP.*diffVelocity) / max(mp*zp, eps);
end

% Phase-2 Function Definition (post-exit)
function [rho2, pressureEqn2, diffZp, diffZeta, diffTheta2] = thermoModelv2(zp, Tgas, zeta)
    % Coefs and Params
    mp = 0.376; s = 9.98e-4; V0 = 373e-6; lm = 2.934; f = 1.071e6;
    eta = 1.064e-3; gamma = 1.2; R = 340; rho_p = 1600;
    r1 = 0.597e-9; S1 = 134.4e-6; LAMBDA1 = 75.2e-9; Kappa1 = 0.755; lambda1 = 0.159;

    % Equation of state in Phase-2 (post-exit), with live z_p:
    denom        = V0 + s*lm - (mp/rho_p)*(1 - zp) - eta*mp*(zp - zeta);
    pressureEqn2 = (mp*(zp - zeta)*R*Tgas) / denom;
    rho2         = (mp*(zp - zeta)) / denom;

    % Choked outflow factor (isentropic)
    CF_CHOKE = sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
    diffZeta = (s*pressureEqn2)/(mp*sqrt(R*Tgas)) * CF_CHOKE;

    % Post-exit burning (non-negative) until z_p -> 1
    if zp < 1 - 1e-12
        diffZp = (S1/LAMBDA1)*r1*pressureEqn2*sqrt(1 + 4*zp*lambda1/Kappa1);
        diffZp = max(diffZp,0);
    else
        diffZp = 0;
    end

    % Energy equation in theta-form
    gas_mass_frac = max(zp - zeta, eps);
    diffTheta2 = ( (f - R*Tgas)*diffZp - (gamma - 1)*R*Tgas*diffZeta ) / gas_mass_frac;
end

%% Final sanity check prints, optional
function S = summarize_ballistics(timeVec, pressureEqn, Tgas, location, velocityP, zp, zetaHist, t_release, t_exit)
    S.t_release = t_release;
    S.t_exit    = t_exit;
    [S.P_peak, iP] = max(pressureEqn);
    S.t_P_peak  = timeVec(iP);
    S.Tg_peak   = max(Tgas);
    S.u_exit    = interp1(timeVec, velocityP, t_exit, 'linear','extrap');
    S.zp_exit   = interp1(timeVec, zp, t_exit, 'linear','extrap');
    S.t_zp1     = timeVec(find(zp>=1-1e-6,1,'first'));
    S.impulse   = trapz(timeVec, pressureEqn);        % Pa·s
    S.t_on      = nan(1, size(hHist,2));
    for k = 1:size(hHist,2)
        ii = find(hHist(:,k)>0,1,'first');
        if ~isempty(ii), S.t_on(k) = timeVec(ii); end
    end
end
