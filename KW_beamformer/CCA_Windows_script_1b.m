% % Modified Kaiser Window beamformer for Circular Array (CA)
% % Effect of beta kaiser 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all ; clc ; close all ;
% 
% design_name = 'CCA_design' ; 
% design = load(design_name) ;
% r_p = design.r_p ;
% phi_p_m = design.phi_p_m ;
% 
% r_p = r_p(3) ; % 10 cm
% phi_p_m = phi_p_m(3) ; % 10 cm
% 
% M_tot = size(phi_p_m{1},1) ;
% P = 1 ;
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
% [ d ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
% [ Gamma_distance ] = GammaDistance_CCA( r_p, phi_p_m ) ;
% 
% % Modified Kaiser Window
% %-------------------------------------------------------------------------------------------
% M_active = 100 ;
% beta_range = 100.^([-inf,-1:1/20:2]') ; 
% ring_weight_p = { ones(1, length(f)) } ; % 1
% 
% b_theta = zeros( length(beta_range), length(f) ) ;
% b_phi = zeros( length(beta_range), length(f) ) ;
% D = zeros( length(beta_range), length(f) ) ;
% W = zeros( length(beta_range), length(f) ) ;
% sensors_weights = zeros( length(beta_range), M_tot ) ;
% 
% for idx_beta_range = 1 : length(beta_range) 
%     
%     beta_kaiser_p = cell(1,P) ;
%     for p = 1 : P
%         beta_kaiser_p{p} = beta_range(idx_beta_range)*ones(1,length(f)) ;
%     end
%     
%     [ h ] = Modified_Bidir_Kaiser_CCA( M_active, beta_kaiser_p, ring_weight_p, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
% 
%     % Sensor-weights
%     tmp_weights = abs( h{1}(:,1) ) ;
%     sensors_weights(idx_beta_range,:) = tmp_weights ;
% 
%     % Metrics
%     F_low = 1000 ; % Hz
%     F_high = 3000 ; % Hz
%     [ D(idx_beta_range,:) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
% 
%     power_level_diff = 6 ;
%     [ b_theta(idx_beta_range,:), b_phi(idx_beta_range,:) ] = BW_CCA(h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
% 
%     [ W(idx_beta_range,:) ] = WNG_CCA(h, d) ;
% 
% end
% 
% save('varying_beta_metrics', 'r_p', 'phi_p_m', 'beta_range', 'sensors_weights', 'b_theta', 'b_phi', 'D', 'W') ;
% 
% exit ;
% 
% return ;

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc ; close all ;
load('varying_beta_metrics') ;

beta_range = 100.^([-inf,-1:1/20:2]') ; 
idx_beta = 1 + [21,31,41] ;
beta = beta_range( idx_beta ) ;
beta_prime = 0.5 * log10( beta ) ; % log_100 (beta)

sensors_angles = phi_p_m{1} * 180/pi ;

figure();
plot( sensors_angles, sensors_weights(idx_beta,:)', '-o' ) ;
xticks([-180:45:180]) ; 
xlim([-180,180]) ; %axis('tight') ;
ylim([ min(ylim)-0, max(ylim)+0]) ; yticks( unique(round( linspace(min(ylim), max(ylim), 5),2 )) ) ;

title(['Kaiser window on a CA']) ; 
xlabel('degrees') ; ylabel('sensor weight') ; 
hleg = legend(num2str(beta_prime)); 
title(hleg, '$\beta_1^{\prime} (f) = \log_{100} \beta_1 (f)$','Interpreter','Latex');

b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');




% % choose frequency
Ts = 1/16000 ; FS = 1/Ts ; % Hz
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;
beta_range = 100.^([-inf,-1:1/20:2]') ; 

idx_freq = 1 + 2.^[3,4] ;
freq = f(idx_freq)*FS/1000 ;

figure();
for idx_row = 1 : length(freq)
    
    idx_tmp = idx_freq(idx_row) ;
    
    idx_plot = 4*(idx_row-1) + 1 ;
    subplot(length(freq),4,idx_plot) ; plot(log10(beta_range)/2, b_theta(:, idx_tmp ) ) ; 
    xticks([-1:0.5:2]) ; xlim([-0.2,1.2]) ; 
    ylim([ min(ylim)-0, max(ylim)+0]) ; yticks( unique(round( linspace(min(ylim), max(ylim), 3) )) ) ;
    if idx_row == 1
        title(['$b_{\theta}(f)$']) ; %axis('tight') ;
    elseif idx_row == length(freq)
        xlabel('$\beta_1^{\prime} (f) = \log_{100} \beta_1 (f)$') ; 
    end
    ylabel('degrees') ;
            
    b=gca;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
    
    


    idx_plot = 4*(idx_row-1) + 2 ;
    subplot(length(freq),4,idx_plot) ; plot(log10(beta_range)/2, b_phi(:, idx_tmp ) ); 
    xticks([-1:0.5:2]) ; xlim([-0.2,1.2]) ; 
    ylim([ min(ylim)-0, max(ylim)+0]) ; yticks( unique(round( linspace(min(ylim), max(ylim), 3) )) ) ;
    if idx_row == 1
        title(['$b_{\phi}(f)$']) ; %axis('tight') ;
    elseif idx_row == length(freq)
        xlabel('$\beta_1^{\prime} (f) = \log_{100} \beta_1 (f)$') ; 
    end
    ylabel('degrees') ;
        
    b=gca;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
    
    


    idx_plot = 4*(idx_row-1) + 3 ;
    subplot(length(freq),4,idx_plot) ; plot(log10(beta_range)/2, 10*log10( D(:, idx_tmp ) ) ) ; 
    xticks([-1:0.5:2]) ; xlim([-0.2,1.2]) ; 
    ylim([ min(ylim)-0, max(ylim)+0]) ; yticks( round( linspace(min(ylim), max(ylim), 3),3 ) ) ;
    if idx_row == 1
        title(['$\mathcal{D} (f)$']) ; %axis('tight') ;
    elseif idx_row == length(freq)
        xlabel('$\beta_1^{\prime} (f) = \log_{100} \beta_1 (f)$') ; 
    end
    ylabel('dB') ;
        
    b=gca;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

    
    

    idx_plot = 4*(idx_row-1) + 4 ;
    subplot(length(freq),4,idx_plot) ; plot(log10(beta_range)/2, 10*log10( W(:, idx_tmp ) ) ) ; 
    xticks([-1:0.5:2]) ; xlim([-0.2,1.2]) ; 
    ylim([ min(ylim)-0, max(ylim)+0]) ; yticks( round( linspace(min(ylim), max(ylim), 3),3 ) ) ;
    if idx_row == 1
        title(['$\mathcal{W} (f)$']) ; %axis('tight') ;
    elseif idx_row == length(freq)
        xlabel('$\beta_1^{\prime} (f) = \log_{100} \beta_1 (f)$') ; 
    end
    ylabel('dB') ;
        
    hleg = legend(num2str(freq(idx_row))); title(hleg, '$f $ (kHz)','Interpreter','Latex');
        
    b=gca;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');


end
