% % DS beamformer for CA and central sensor
% % Effect of radius 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all ; clc ; close all ;
% 
% design_name = 'CCA_design' ; 
% design = load(design_name) ;
% r_p = design.r_p ;
% phi_p_m = design.phi_p_m ;
% 
% r_p = r_p([1:2:5]') ; % 0, 10, 20 cm
% phi_p_m = phi_p_m([1:2:5]') ; % 0, 10, 20 cm
% 
% 
% c = 340 ; Ts = 1/16000 ; FS = 1/Ts ; 
% f_max = FS / 2 ;
% lambda_min = c / f_max ;
% 
% % % choose frequency
% f = [0 : FS/256 : FS/2]' ; % Hz
% f = f/FS ;
% 
% theta_d = 45 ; % Evevation - DOA of the SOI between [0,90]
% phi_d = 45 ; % Azimuth - DOA of the SOI between [0,90]
% 
% 
% % DS
% %-------------------------------------------------------------------------------------------
% 
% b_theta = zeros( length(r_p), length(f) ) ;
% b_phi = zeros( length(r_p), length(f) ) ;
% D = zeros( length(r_p), length(f) ) ;
% W = zeros( length(r_p), length(f) ) ;
% 
% for idx_r_p = 1 : length(r_p) 
% 
%     [ d ] = d_CCA( r_p(idx_r_p), phi_p_m(idx_r_p), theta_d, phi_d, f, c, Ts ) ;
%     [ Gamma_distance ] = GammaDistance_CCA( r_p(idx_r_p), phi_p_m(idx_r_p) ) ;
%     
%     [ h ] = DS_CCA( d ) ;
%     
%     % Metrics
%     F_low = 1000 ; % Hz
%     F_high = 3000 ; % Hz
%     [ D(idx_r_p,:) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
% 
%     power_level_diff = 6 ;
%     [ b_theta(idx_r_p,:), b_phi(idx_r_p,:) ] = BW_CCA(h, r_p(idx_r_p), phi_p_m(idx_r_p), theta_d, phi_d, f, c, Ts, power_level_diff) ;
% 
%     [ W(idx_r_p,:) ] = WNG_CCA(h, d) ;
% 
% end
% 
% save('varying_radius_metrics', 'r_p', 'phi_p_m', 'b_theta', 'b_phi', 'D', 'W') ;
% 
% exit ;
% 
% return ;

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc ; close all ;

load(['varying_radius_metrics']) ;

% % choose frequency
Ts = 1/16000 ; FS = 1/Ts ; % Hz
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;

choice = 'normal' ;
if strcmp( choice , 'logscale' )
    frequencies = log2(f*FS) / log2(2) ; %log2
    fig_ticks = frequencies([1,4+1, 16+1, 32+1, 64+1, end]) ;
    fig_labels = 2.^( fig_ticks ) / 1000 ;
elseif strcmp( choice , 'normal' )
    frequencies = f*FS ;
    fig_ticks = linspace(0, FS/2, 5)' ;
    fig_labels = fig_ticks/1000 ;
end

figure();
subplot(1,4,1) ; 
values = b_theta' ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$b_{\theta}(f)$']) ; 
xlabel('$f$ (kHz)') ; ylabel('degrees') ; %axis('tight') ;
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ min(ylim)-0, 360]) ; 
yticks( sort(  unique( [ 0, round( linspace(min(ylim), max(ylim), 4) ) ] )  ) ) ;
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,2) ; 
values = b_phi' ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$b_{\phi}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('degrees') ; %axis('tight') ; 
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ min(ylim)-0, 360]) ; 
yticks( sort(  unique( [ 0, round( linspace(min(ylim), max(ylim), 4) ) ] )  ) ) ;
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,3) ; 
values = 10*log10( D' ) ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$\mathcal{D}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('dB') ; %axis('tight') ;
xlim([ min(f*FS)+10, max(f*FS)]) ;  
ylim([ 0, 20]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,4) ; 
values = 10*log10( W' ) ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$\mathcal{W}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('dB') ; %axis('tight') ;
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ 0, 22]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
hleg = legend([num2str(r_p*100)],'Interpreter','Latex' ); 
title(hleg, ['CA (Nyquist)~$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = 45^o$~radius (cm) = '],'Interpreter','Latex');
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

