% % Modified Kaiser Window beamformer for Concentric Circular Array (UCCA)
% % Using Gradient Descent Algorithm 
% % Convergence - number of iterations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all ; clc ; close all ;
% 
% design_name = 'CCA_design' ; 
% design = load(design_name) ;
% r_p = design.r_p ;
% phi_p_m = design.phi_p_m ;
% 
% clear design ;
% 
% active_rings = zeros(1,length(phi_p_m)) ;
% for p = 1 : length(phi_p_m) 
%     active_rings(p) = not( isempty(phi_p_m{p}) ) + 0 ; 
% end
% r_p = r_p( find(active_rings) ) ;
% phi_p_m = phi_p_m( find(active_rings) ) ;
% P = length(r_p) ;
% 
% c = 340 ; Ts = 1/16000 ; FS = 1/Ts ; 
% 
% % % choose frequency
% f = [0 : FS/256 : FS/2]' ; % Hz
% f = f/FS ;
% 
% theta_d = 45 ; % Elevation - DOA of the SOI between [0,90]
% phi_d = 45 ; % Azimuth - DOA of the SOI between [0,90]
% 
% % Modified Kaiser Window
% %-------------------------------------------------------------------------------------------
% M_active = 100 ;
% init_beta_prime = cell(1,P) ;
% for p = 1 : P
%     if r_p(p) == 0
%         init_beta_prime{p} = nan*ones(1,length(f)) ;
%     else
%         init_beta_prime{p} = 0.5*ones(1,length(f)) ;
%     end
% end
% 
% init_w = cell(1,P) ;
% for p = 1 : P
%     init_w{p} = ( 1/P ) * ones(1, length(f)) ;
% end
% 
% num_iterations = 41 ;
% gradient_stepsize = 0.2*ones(1,length(f)) ;
% mu_beta_prime = 0.1*ones(1,length(f)) ; mu_w = 0.1*ones(1,length(f)) ;
% theta_BW = 40 ; 
% phi_BW = 40 ;
% power_level_diff = 6 ;
% 
% try
%     
% [ h, beta_prime, w, theta_range, phi_range, Delta_b_theta, Delta_b_phi, D, W] = GradDes_Kaiser_CCA( M_active, init_beta_prime, init_w, num_iterations, gradient_stepsize, mu_beta_prime, mu_w, theta_BW, phi_BW, power_level_diff, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
% 
% save([design_name, '_optimize_bw_D'], 'h', 'beta_prime', 'w', 'theta_range', 'phi_range', 'Delta_b_theta', 'Delta_b_phi', 'D', 'W' ) ;
% 
% exit ;
% 
% catch
%     exit ;
% end
% 
% return ;

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc ; close all ;

design_name = 'CCA_design' ; 
load([design_name, '_optimize_bw_D']) ;

% % choose frequency
Ts = 1/16000 ; FS = 1/Ts ; % Hz
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;
P = length(h) ;
num_iterations = 41 ;

theta_BW = 40 ; 
phi_BW = 40 ;

idx_freq = find( or(f*FS == 1000, f*FS == 6000 ) ) ; 
freq = f(idx_freq)*FS/1000 ;

figure();
subplot(1,4,1) ; plot([0:(num_iterations-1)]', Delta_b_theta(:,idx_freq)+theta_BW ) ; title(['$b_{\theta}(f)$']) ;
xlabel('iterations ($n$)') ; ylabel('degrees') ; 
xlim([ 0, (num_iterations-1)]) ; 
ylim([ 0, 360]) ; 
yticks( sort(  unique( [ 0, theta_BW, round( linspace(min(ylim), max(ylim), 5) ) ] )  ) ) ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,2) ; plot([0:(num_iterations-1)]', Delta_b_phi(:,idx_freq)+phi_BW ) ; title(['$b_{\phi}(f)$']) ;
xlabel('iterations ($n$)') ; ylabel('degrees') ; 
xlim([ 0, (num_iterations-1)]) ; 
ylim([ 0, 360]) ; 
yticks( sort(  unique( [ 0, phi_BW, round( linspace(min(ylim), max(ylim), 5) ) ] )  ) ) ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,3) ; plot([0:(num_iterations-1)]', 10*log10( D(:,idx_freq) ) ) ; title(['$\mathcal{D}(f)$']) ;
xlabel('iterations ($n$)') ; ylabel('dB') ;
xlim([ 0, (num_iterations-1)]) ; 
ylim([ min(ylim)-0, max(ylim)+0]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,4) ; plot([0:(num_iterations-1)]', 10*log10( W(:,idx_freq) ) ) ; title(['$\mathcal{W}(f)$']) ;
xlabel('iterations ($n$)') ; ylabel('dB') ;
xlim([ 0, (num_iterations-1)]) ; 
ylim([ min(ylim)-0, max(ylim)+0]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
hleg = legend(num2str(freq)); 
title(hleg, ['CCA (Nyquist)~$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = 45^o$~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$~$f$ (kHz) = '],'Interpreter','Latex');
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

return;
%-------------------------------------------------------------------------------------------------------------------------

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
values = Delta_b_theta(end,:)' + theta_BW ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$b_{\theta}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('degrees') ; axis('tight') ;
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ 0, 360]) ; 
yticks( sort(  unique( [ 0, theta_BW, round( linspace(min(ylim), max(ylim), 5) ) ] )  ) ) ;
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,2) ; 
values = Delta_b_phi(end,:)' + phi_BW ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$b_{\phi}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('degrees') ; axis('tight') ; 
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ 0, 360]) ; 
yticks( sort(  unique( [ 0, phi_BW, round( linspace(min(ylim), max(ylim), 5) ) ] )  ) ) ;
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,4,3) ; 
values = 10*log10( D(end,:)' ) ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$\mathcal{D}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('dB') ; axis('tight') ;
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
values = 10*log10( W(end,:)' ) ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
plot(frequencies, values ) ; title(['$\mathcal{W}(f)$']) ;
xlabel('$f$ (kHz)') ; ylabel('dB') ; axis('tight') ;
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ 0, 22]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
hleg = legend(''); 
title(hleg, ['CCA (Nyquist)~$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = 45^o$~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$'],'Interpreter','Latex');
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');


return;
%-------------------------------------------------------------------------------------------------------------------------

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

linestyles = {':+r', '-b', '--c', '-.m', ':g'} ;

figure() ;


for p = 1 : length(w)
values = w{p} ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
subplot(2,1,1) ; plot(frequencies, values, linestyles{p} ) ; hold on ; axis('tight') ;
end
xlabel('$f$ (kHz)') ; ylabel('$w_p (f)$') ;
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ min(ylim)-0, max(ylim)+0]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 2 ) ) ;
hleg = legend(num2str([0:P-1]')); 
title(hleg, ['CCA (Nyquist)~$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = 45^o$~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$~$p$ (ring) = '],'Interpreter','Latex');
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

for p = 2 : length(beta_prime)
values = beta_prime{p} ;
values = movmean(movmean(movmedian(values,11), 11), 11) ;
subplot(2,1,2) ; plot(frequencies, values, linestyles{p} ) ; hold on ; 
end
xlabel('$f$ (kHz)') ; ylabel('$\beta_p^{\prime} (f)$') ; axis('tight') ;
xlim([ min(f*FS)+10, max(f*FS)]) ; 
ylim([ min(ylim)-0, max(ylim)+0]) ; 
yticks( round( linspace(min(ylim), max(ylim), 5), 2 ) ) ;
hleg = legend(num2str([1:P-1]')); 
title(hleg, ['CCA (Nyquist)~$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = 45^o$~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$~$p$ (ring) = '],'Interpreter','Latex');
b=gca;
set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
