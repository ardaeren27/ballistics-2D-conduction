%% main_repeated_rifling.m
% Build repeated BC and launch the patched 2-D solver (bc-enabled).

clear; clc;

% ------------ knobs ------------
steel      = 'DUPLEX';
tCr_um     = 0;
Nr         = 240*4;
Nz         = 120;
dt_fd      = 5e-4;
thetaFD    = 1.0;
debug      = true;
debug_strd = 10;
keep_stride= 20;

% repeated schedule
Nshots   = 7;
cool_gap = 0.1;   % seconds
Tamb_C   = 20;
use_30col= false; % false unless new array specified

% ------------ build BC ------------
bc = repeated_rifling('Nshots',Nshots,'cool_gap',cool_gap,'Tamb_C',Tamb_C, ...
                      'use_30col',use_30col,'plot',true,'smoke_test',true);

% ------------ launch solver ------------
tEnd_fd = bc.t(end);
fprintf('Launching 2-D FV solver for %.3f s total...\n', tEnd_fd);

out = heat_transfer_2d_solver('steel',steel,'tCr_um',tCr_um, ...
      'Nr',Nr,'Nz',Nz,'dt_fd',dt_fd,'tEnd_fd',tEnd_fd, ...
      'theta',thetaFD,'debug',debug,'debug_stride',debug_strd, ...
      'keep_stride',keep_stride,'plot',true, ...
      'bc',bc);

% ------------ quick post ------------
figure('Name','Inner-surface T','Color','w'); grid on; hold on;
plot(1e3*out.t, out.T_inner6,'LineWidth',1.2);
xlabel('time [ms]'); ylabel('T_{inner} [°C]');
title('Inner-surface temperature at sections');
legend(out.sections.names,'Location','best');
