%SINDy average relative errors
close all; 
F1 = figure

RR
WRR
IT
WIT
SP
WSP
SPGL1
WSPGL1

semilogy(...
    IT_stats(:,1),IT_stats(:,2),'b*:', IT_stats(:,1),WIT_stats(:,2),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
hold on
semilogy(...
    IT_stats(:,1),SP_stats(:,2),'o:', IT_stats(:,1),WSP_stats(:,2),'o-',...
    'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
semilogy(...
    IT_stats(:,1),RR_stats(:,2),'rs:',IT_stats(:,1),WRR_stats(:,2),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
semilogy(...
    IT_stats(:,1),SPGL1_stats(:,2),'kd:',IT_stats(:,1),WSPGL1_stats(:,2),'kd-',...
    'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')

   
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);

grid on
% 
 L = legend('STLS','STLS$_{RPMO}$','SP','SP$_{RPMO}$','STRidge','STRidge$_{RPMO}$','BPDN','BPDN$_{RPMO}$','location','best')
 set(L,'interpreter','latex','fontsize',28); 
xlim([IT.N(1),IT.N(end)])
% ylim([min([LASSO_stats(:,2); WLASSO_stats(:,2); IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2); SP_stats(:,2); WSP_stats(:,2) ]),...
%     max([LASSO_stats(:,2); WLASSO_stats(:,2); IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2);SP_stats(:,2); WSP_stats(:,2) ])])
ylim([min([IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2); SP_stats(:,2); WSP_stats(:,2) ]),...
    max([ IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2);  SP_stats(:,2); WSP_stats(:,2) ])])

xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
ylabel('$\bar{e}_{ref}$','interpreter','latex','fontsize',24) 

% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
set(gca,'FontSize',28)



%%

%save_name = [func2str(func), '_1percent_err']
save_name = [func2str(func), '_1percent_err']
saveas(F1,save_name,'epsc')



%% The standard deviation of relative errors
close all; 
F1 = figure

semilogy(...
    IT_stats(:,1),IT_stats(:,end-2),'b*:', IT_stats(:,1),WIT_stats(:,end-2),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
hold on
semilogy(...
    IT_stats(:,1),SP_stats(:,end-2),'o:', IT_stats(:,1),WSP_stats(:,end-2),'o-',...
    'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
semilogy(...
    IT_stats(:,1),RR_stats(:,end-2),'rs:',IT_stats(:,1),WRR_stats(:,end-2),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
semilogy(...
    IT_stats(:,1),SPGL1_stats(:,end-2),'kd:',IT_stats(:,1),WSPGL1_stats(:,end-2),'kd-',...
    'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')

   
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);

grid on

L = legend('STLS','STLS$_{RPMO}$','SP','SP$_{RPMO}$','STRidge','STRidge$_{RPMO}$','BPDN','BPDN$_{RPMO}$','location','best')
set(L,'interpreter','latex','fontsize',28); 
xlim([IT.N(1),IT.N(end)])
ylim([min([IT_stats(:,end-2);WIT_stats(:,end-2); RR_stats(:,end-2);WRR_stats(:,end-2); SPGL1_stats(:,end-2);WSPGL1_stats(:,end-2); SP_stats(:,end-2); WSP_stats(:,end-2) ]),...
    max([ IT_stats(:,end-2);WIT_stats(:,end-2); RR_stats(:,end-2);WRR_stats(:,end-2); SPGL1_stats(:,end-2);WSPGL1_stats(:,end-2);  SP_stats(:,end-2); WSP_stats(:,end-2) ])])

xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
ylabel('$\sigma_{e_{ref}}$','interpreter','latex','fontsize',24) 

% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
set(gca,'FontSize',28)
%%

%save_name = [func2str(func), '_1percent_err_std']
save_name = [func2str(func), '_1percent_err_std']
saveas(F1,save_name,'epsc')



%% This is for SINDy (relative fit  error)
%close all; 
F1 = figure
semilogy(...
    IT_stats(:,1),IT_stats(:,4),'b*:', IT_stats(:,1),WIT_stats(:,4),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
hold on
semilogy(...
    IT_stats(:,1),SP_stats(:,4),'o:', IT_stats(:,1),WSP_stats(:,4),'o-',...
    'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
semilogy(...
    IT_stats(:,1),RR_stats(:,4),'rs:',IT_stats(:,1),WRR_stats(:,4),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
semilogy(...
    IT_stats(:,1),SPGL1_stats(:,4),'kd:',IT_stats(:,1),WSPGL1_stats(:,4),'kd-',...
    'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')

   
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);    
 %L = legend('STLS','STLS$_{RPMO}$','SP','SP$_{RPMO}$','STRidge','STRidge$_{RPMO}$','BPDN','BPDN$_{RPMO}$','location','best')

grid on


%set(L,'interpreter','latex','fontsize',28); 
xlim([IT.N(1),IT.N(end)])
ylim([min([IT_stats(:,4);WIT_stats(:,4); RR_stats(:,4);WRR_stats(:,4); SPGL1_stats(:,4);WSPGL1_stats(:,4); SP_stats(:,4); WSP_stats(:,4) ]),...
    max([ IT_stats(:,4);WIT_stats(:,4); RR_stats(:,4);WRR_stats(:,4); SPGL1_stats(:,4);WSPGL1_stats(:,4);  SP_stats(:,4); WSP_stats(:,4) ])])

xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
ylabel('$\bar{e}_{fit}$','interpreter','latex','fontsize',24) 

% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
set(gca,'FontSize',28)
%%
%save_name = [func2str(func), '_1percent_fit']
save_name = [func2str(func), '_1percent_fit']
saveas(F1,save_name,'epsc')




%%
close all
%This is for SINDy Validation error relative fit standard dev
F1 = figure
semilogy(...
    IT_stats(:,1),IT_stats(:,end),'b*:', IT_stats(:,1),WIT_stats(:,end),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
hold on
semilogy(...
    IT_stats(:,1),SP_stats(:,end),'o:', IT_stats(:,1),WSP_stats(:,end),'o-',...
    'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
semilogy(...
    IT_stats(:,1),RR_stats(:,end),'rs:',IT_stats(:,1),WRR_stats(:,end),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
semilogy(...
    IT_stats(:,1),SPGL1_stats(:,end),'kd:',IT_stats(:,1),WSPGL1_stats(:,end),'kd-',...
    'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')

   
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);

grid on

L = legend('STLS','STLS$_{RPMO}$','SP','SP$_{RPMO}$','STRidge','STRidge$_{RPMO}$','BPDN','BPDN$_{RPMO}$','location','best')
set(L,'interpreter','latex','fontsize',28); 
xlim([IT.N(1),IT.N(end)])
ylim([min([IT_stats(:,end);WIT_stats(:,end); RR_stats(:,end);WRR_stats(:,end); SPGL1_stats(:,end);WSPGL1_stats(:,end); SP_stats(:,end); WSP_stats(:,end) ]),...
    max([ IT_stats(:,end);WIT_stats(:,end); RR_stats(:,end);WRR_stats(:,end); SPGL1_stats(:,end);WSPGL1_stats(:,end);  SP_stats(:,end); WSP_stats(:,end) ])])

xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
ylabel('$\sigma_{e_{fit}}$','interpreter','latex','fontsize',24) 

% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
set(gca,'FontSize',28)

%%
%save_name = [func2str(func), '_1percent_fit_std']
save_name = [func2str(func), '_1percent_fit_std']
saveas(F1,save_name,'epsc')


%%
% 
% 
% close all
% 
% ylim_min = min([WSPGL1_stats(:,8); SPGL1_stats(:,8)])
% ylim_max = max([WSPGL1_stats(:,8); SPGL1_stats(:,8)])
% figure
% 
% semilogy(WSPGL1_stats(:,1),WSPGL1_stats(:,8),'^-',SPGL1_stats(:,1),SPGL1_stats(:,8),'s-','markersize',12,'linewidth',3);
% grid on
% xlim([IT.N(1),IT.N(end)]);ylim([ylim_min, ylim_max])
% xlabel('Time step  $i = 1, \ldots, N$','interpreter','latex','fontsize',24)
% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
% set(gca,'FontSize',28)
% 
% 
% L = legend('$\texttt{cond}(\mathbf{W\Theta})$','$\texttt{cond}(\mathbf{\Theta})$','location','northeast')
% set(L,'interpreter','latex','fontsize',28)
% 
% %%
% save_name = [func2str(func)]
% 
% print_name = [save_name,'_cond']
% print(print_name,'-dpng','-r300')
% 
% 
%%
close all



%ylim_min = min([WSPGL1_stats(:,5) ; SPGL1_stats(:,5); WSPGL1_stats(:,6); SPGL1_stats(:,6)])
%ylim_max = max([WSPGL1_stats(:,5); SPGL1_stats(:,5); WSPGL1_stats(:,6); SPGL1_stats(:,6)])


F1 = figure

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
semilogy(IT_stats(:,1),IT_stats(:,5),'r^:',IT_stats(:,1),IT_stats(:,6),'rs:', 'markersize',36,'linewidth',10, 'markerfacecolor',[0,200,0]/255,'color',[0,200,0]/255);
hold on


semilogy(WIT_stats(:,1),WIT_stats(:,5),'k^-', WIT_stats(:,1),WIT_stats(:,6),'bs-',  'markersize',36,'linewidth',10, 'color',[0,20,137]*1.2/255);
semilogy(WSP_stats(:,1),WSP_stats(:,5),'k^-',WSP_stats(:,1),WSP_stats(:,6),'ks-', 'Color',[250,70,22]/255 ,'markersize',36,'linewidth',10);
semilogy(WRR_stats(:,1),WRR_stats(:,5),'k^-',WRR_stats(:,1),WRR_stats(:,6),'ks-','markersize',36,'linewidth',10,'color',[.8118    0.7216    0.4863]);
semilogy(WSPGL1_stats(:,1),WSPGL1_stats(:,5),'b^-',WSPGL1_stats(:,1),WSPGL1_stats(:,6),'bs-', 'markersize',36,'linewidth',10,'Color','k')
hold on

grid on
xlim([IT.N(1),IT.N(end)]);%ylim([ylim_min, ylim_max])
   
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
%L = legend('STLS','STLS$_{opt}$','SP','SP$_{opt}$','STRidge','STRidge$_{opt}$','BPDN','BPDN$_{opt}$','location','best')
xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
grid on


L = legend('$\mu(\mathbf{\Theta)}$','$\mu_{avg}(\mathbf{\Theta)}$','$\mu^{STLS}(\mathbf{D)}$','$\mu^{STLS}_{avg}(\mathbf{D)}$','$\mu^{SP}(\mathbf{D)}$','$\mu^{SP}_{avg}(\mathbf{D)}$','$\mu^{STRidge}(\mathbf{D)}$','$\mu^{STRidge}_{avg}(\mathbf{D)}$','$\mu^{BPDN}(\mathbf{D)}$','$\mu^{BPDN}_{avg}(\mathbf{D)}$','location','eastoutside')
set(L,'interpreter','latex','fontsize',60)
set(gca,'FontSize',60)

%%
save_name = [func2str(func)]
print_name = [save_name,'_mu']
saveas(F1,print_name,'epsc')


%%

% %%
% 
% save_name = 'lorenz63_1percent_val_std'
% saveas(F1,save_name,'epsc')
% 
% 
% 
%%
%SINDy average relative errors
close all; 
F1 =figure

RR
WRR
IT
WIT
SP
WSP
SPGL1
WSPGL1


F2 = boxplot([IT_err,RR_err,SP_err,SPGL1_err])

% semilogy(...
%     IT_stats(:,1),IT_stats(:,2),'b*:', IT_stats(:,1),WIT_stats(:,2),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
% hold on
% semilogy(...
%     IT_stats(:,1),SP_stats(:,2),'o:', IT_stats(:,1),WSP_stats(:,2),'o-',...
%     'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
% semilogy(...
%     IT_stats(:,1),RR_stats(:,2),'rs:',IT_stats(:,1),WRR_stats(:,2),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
% semilogy(...
%     IT_stats(:,1),SPGL1_stats(:,2),'kd:',IT_stats(:,1),WSPGL1_stats(:,2),'kd-',...
%     'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')
% 
%    
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
    
%semilogy(sample_sizes,SP_stats(:,1),'*-',sample_sizes,GDSP_stats(:,1),'o-', sample_sizes, DSP_stats(:,1),'d-',sample_sizes,ORACLE_stats(:,1),'h-','markersize',12,'linewidth',3);
%L = legend('LASSO','LASSO$_{opt}$','IT','IT$_{opt}$','STRidge','STRidge$_{opt}$','SPG','SPGL$1_{opt}$','SP','WSP','location','best')
%L = legend('STLS','STLS$_{opt}$','SP','SP$_{opt}$','STRidge','STRidge$_{opt}$','BPDN','BPDN$_{opt}$','location','best')

grid on


%set(L,'interpreter','latex','fontsize',28); 
%xlim([IT.N(1),IT.N(end)])
% ylim([min([LASSO_stats(:,2); WLASSO_stats(:,2); IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2); SP_stats(:,2); WSP_stats(:,2) ]),...
%     max([LASSO_stats(:,2); WLASSO_stats(:,2); IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2);SP_stats(:,2); WSP_stats(:,2) ])])
%ylim([min([IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2); SP_stats(:,2); WSP_stats(:,2) ]),...
%    max([ IT_stats(:,2);WIT_stats(:,2); RR_stats(:,2);WRR_stats(:,2); SPGL1_stats(:,2);WSPGL1_stats(:,2);  SP_stats(:,2); WSP_stats(:,2) ])])

%xlabel('interpreter','latex','fontsize',24)
ylabel('$e_{ref}$','interpreter','latex','fontsize',24) 


set(gca,'FontSize',28,'XTickLabel',{'STLS','STRidge','SP','BPDN'});  
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
set(gca,'Yscale','log')
set(F2,'LineWidth',4)
set(F2,'MarkerSize',12)

%%
save_name = [func2str(func)]
print_name = [save_name,'_large_rel_err']
saveas(F1,print_name,'epsc')


%%

%set(gca,'defaulttextinterpreter','latex');
%set(gca,'defaultLegendInterpreter','latex');
% 
% 
% %%
% %This is for SINDy Validation error
% %close all; 
% F1 = figure
% semilogy(...
%     IT_stats(:,1),IT_stats(:,3),'b*:', IT_stats(:,1),WIT_stats(:,3),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
% hold on
% semilogy(...
%     IT_stats(:,1),SP_stats(:,3),'o:', IT_stats(:,1),WSP_stats(:,3),'o-',...
%     'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
% semilogy(...
%     IT_stats(:,1),RR_stats(:,3),'rs:',IT_stats(:,1),WRR_stats(:,3),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
% semilogy(...
%     IT_stats(:,1),SPGL1_stats(:,3),'kd:',IT_stats(:,1),WSPGL1_stats(:,3),'kd-',...
%     'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')
% 
%    
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
% L = legend('STLS','STLS$_{opt}$','SP','SP$_{opt}$','STRidge','STRidge$_{opt}$','BPDN','BPDN$_{opt}$','location','best')
% 
% grid on
% 
% 
% set(L,'interpreter','latex','fontsize',28); 
% xlim([IT.N(1),IT.N(end)])
% ylim([min([IT_stats(:,3);WIT_stats(:,3); RR_stats(:,3);WRR_stats(:,3); SPGL1_stats(:,3);WSPGL1_stats(:,3); SP_stats(:,3); WSP_stats(:,3) ]),...
%     max([ IT_stats(:,3);WIT_stats(:,3); RR_stats(:,3);WRR_stats(:,3); SPGL1_stats(:,3);WSPGL1_stats(:,3);  SP_stats(:,3); WSP_stats(:,3) ])])
% 
% xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
% ylabel('Validation Error $e_{val}$','interpreter','latex','fontsize',24) 
% 
% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
% set(gca,'FontSize',28)
% 
% 
% %%
% 
% save_name = 'lorenz63_1percent_val'
% saveas(F1,save_name,'epsc')
% 
% %%
% %This is for SINDy Validation standard deviation
% %close all; 
% F1 = figure
% semilogy(...
%     IT_stats(:,1),IT_stats(:,end-1),'b*:', IT_stats(:,1),WIT_stats(:,end-1),'b*-', 'markersize',16,'linewidth',4, 'markerfacecolor',[0,20,137]*1.2/255,'color',[0,20,137]*1.2/255);
% hold on
% semilogy(...
%     IT_stats(:,1),SP_stats(:,end-1),'o:', IT_stats(:,1),WSP_stats(:,end-1),'o-',...
%     'Color',[250,70,22]/255 ,'markersize',16,'linewidth',4,'MarkerFacecolor',[250,70,22]/255);
% semilogy(...
%     IT_stats(:,1),RR_stats(:,end-1),'rs:',IT_stats(:,1),WRR_stats(:,end-1),'rs-','markerfacecolor',[.8118    0.7216    0.4863],'markersize',16,'linewidth',4,'color',[.8118    0.7216    0.4863]);
% semilogy(...
%     IT_stats(:,1),SPGL1_stats(:,end-1),'kd:',IT_stats(:,1),WSPGL1_stats(:,end-1),'kd-',...
%     'markerfacecolor','k', 'markersize',16,'linewidth',4,'Color','k')
% 
%    
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
% L = legend('STLS','STLS$_{opt}$','SP','SP$_{opt}$','STRidge','STRidge$_{opt}$','BPDN','BPDN$_{opt}$','location','best')
% 
% grid on
% 
% 
% set(L,'interpreter','latex','fontsize',28); 
% xlim([IT.N(1),IT.N(end)])
% ylim([min([IT_stats(:,end-1);WIT_stats(:,end-1); RR_stats(:,end-1);WRR_stats(:,end-1); SPGL1_stats(:,end-1);WSPGL1_stats(:,end-1); SP_stats(:,end-1); WSP_stats(:,end-1) ]),...
%     max([ IT_stats(:,end-1);WIT_stats(:,end-1); RR_stats(:,end-1);WRR_stats(:,end-1); SPGL1_stats(:,end-1);WSPGL1_stats(:,end-1);  SP_stats(:,end-1); WSP_stats(:,end-1) ])])
% 
% xlabel('Number of samples $N$','interpreter','latex','fontsize',24)
% ylabel('Standard Deviation $\sigma_{e_{val}}$','interpreter','latex','fontsize',24) 
% 
% if noisy_data 
%     title_string = ['Noisy data ','$\eta = \;$', num2str(eta)]
% elseif noisy_measurements
%     title_string = ['Noisy derivatives ','$\eta = \;$', num2str(eta)]
% else
%     title_string = ['Noiseless ','$\eta = \;$', num2str(eta)]
% end
%     
% title(title_string,'interpreter','latex','fontsize',24);
% set(gca,'FontSize',28)
% 
% 
