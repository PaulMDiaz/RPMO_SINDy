clear all
load('120_MS_movie_data.mat')


%%
%X0 = X0_ref_scaled;
close all

max_dt =  max([max(dTA_WSPGL1),max(dTA_WSP),max(dTA_WRR),max(dTA_WIT),max(dTA_SPGL1),max(dTA_SP),max(dTA_RR),max(dTA_IT)]);

N_frames = 180;
Time_percents = linspace(0.001,1,N_frames);
fontsize = 28
val_max = size(X_data_val,1);
%fig = figure('Renderer', 'painters', 'Position', [10 10 1400 800])
fig = figure('Position', [10 10 1400 800])
%hold off
% axis equal
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% C = fig.CData;
% fig.CDataMapping = 'scaled';
F(N_frames) = struct('cdata',[],'colormap',[]);
%%

for frame_count = 1:N_frames
    Tper = Time_percents(frame_count)

    if strcmp(precond_type,'robust')
    %T = [T0:dt:2*TF]; 
    
    Tmax = length(dTA_WSPGL1);
    WSPGL1_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_WSPGL1_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,1)
    color_line3(X_WSPGL1(1:floor(Tper*Tmax),1),X_WSPGL1(1:floor(Tper*Tmax),2),X_WSPGL1(1:floor(Tper*Tmax),3),dTA_WSPGL1(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on;
    title_string = [{'BPDN$_{opt}$'}, {['$e_{ref}=$',num2str(WSPGL1_stats(1:n,2))]}, {['$e_{val}=$',num2str(WSPGL1_validation_error)]}];
%     title_string = [{'BPDN$_{opt}$'}, {['$e_{ref}=$',num2str(WSPGL1_stats(1:n,2)), ' $e_{val}=$',num2str(WSPGL1_validation_error)]}, ...
%                {'$\lambda_{avg} = $', num2str(mean(sum(abs(WSPGL1_Lyap(1:floor(Tper*val_max),:)))))} ...
%            ];
    %title_string = ['BPDN$_{opt}$,  $e_{val}=$',num2str(WSPGL1_validation_error)];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_WSPGL1);
    view([-45,45,45])

    
    Tmax = length(dTA_WRR);
    WRR_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_WRR_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,2)
    color_line3(X_WRR(1:floor(Tper*Tmax),1),X_WRR(1:floor(Tper*Tmax),2),X_WRR(1:floor(Tper*Tmax),3),dTA_WRR(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
     title_string = [{'STRidge$_{RPMO}$'}, {['$e_{ref}=$',num2str(WRR_stats(1:n,2))]}, {[' $e_{val}=$',num2str(WRR_validation_error)]}];
%  title_string = [{'STRidge$_{opt}$'}, {['$e_{ref}=$',num2str(WRR_stats(1:n,2)), ' $e_{val}=$',num2str(WRR_validation_error)]}, ...
%                 {'$\lambda_{avg} = $', num2str(mean(sum(abs(WRR_Lyap(1:floor(Tper*val_max),:)))))} ...
%             ];    
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)
    set(gca,'FontSize',fontsize)
    length(TT_WRR);
    view([-45,45,45])

    %T = [0,TF]; 
    Tmax = length(dTA_WIT);
    WIT_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_WIT_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,3)
    color_line3(X_WIT(1:floor(Tper*Tmax),1),X_WIT(1:floor(Tper*Tmax),2),X_WIT(1:floor(Tper*Tmax),3),dTA_WIT(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
     title_string = [{'STLS$_{RPMO}$'}, {['$e_{ref}=$',num2str(WIT_stats(1:n,2))]},{[' $e_{val}=$',num2str(WIT_validation_error)]}];
%  title_string = [{'STLS$_{opt}$'}, {['$e_{ref}=$',num2str(WIT_stats(1:n,2)), ' $e_{val}=$',num2str(WIT_validation_error)]}, ...
%                 {'$\lambda_{avg} = $', num2str(mean(sum(abs(WIT_Lyap(1:floor(Tper*val_max),:)))))} ...
%             ];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_WIT);
    view([-45,45,45])

    %T = [0,TF];
    Tmax = length(dTA_WSP);
    WSP_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_WSP_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,4)
    color_line3(X_WSP(1:floor(Tper*Tmax),1),X_WSP(1:floor(Tper*Tmax),2),X_WSP(1:floor(Tper*Tmax),3),dTA_WSP(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'SP$_{RPMO}$'}, {['$e_{ref}=$',num2str(WSP_stats(1:n,2))]},{['$e_{val}=$',num2str(WSP_validation_error)]}];
    %title_string = ['SP$_{opt}$,  $e_{val}=$',num2str(WSP_validation_error)];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_WSP);
    view([-45,45,45])
    %T = [0,TF]; 

    Tmax = length(dTA_SPGL1);
    subplot(2,4,5)
    if size(X_data_val,1) == size(X_SPGL1_val,1)
        SPGL1_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_SPGL1_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
        color_line3(X_SPGL1(1:floor(Tper*Tmax),1),X_SPGL1(1:floor(Tper*Tmax),2),X_SPGL1(1:floor(Tper*Tmax),3),dTA_SPGL1(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    else
        SPGL1_validation_error = Inf;
    end    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'BPDN'}, {['$e_{ref}=$',num2str(SPGL1_stats(1:n,2))]}, {['$e_{val}=$',num2str(SPGL1_validation_error)]}];
%      title_string = [{'BPDN'}, {['$e_{ref}=$',num2str(SPGL1_stats(1:n,2)), ' $e_{val}=$',num2str(SPGL1_validation_error)]}, ...
%                 {'$\lambda_{avg} = $', num2str(mean(sum(abs(SPGL1_Lyap(1:floor(Tper*val_max),:)))))} ...
%             ];
    %title_string = ['BPDN,  $e_{val}=$',num2str(SPGL1_validation_error)];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_SPGL1);
    view([-45,45,45])

    %T = [0,TF]; 

    Tmax = length(dTA_RR);
    RR_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_RR_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,6)
    color_line3(X_RR(1:floor(Tper*Tmax),1),X_RR(1:floor(Tper*Tmax),2),X_RR(1:floor(Tper*Tmax),3),dTA_RR(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'STRidge'}, {['$e_{ref}=$',num2str(RR_stats(1:n,2))]}, {['$e_{val}=$',num2str(RR_validation_error)]}];
    %title_string = ['STRidge,  $e_{val}=$',num2str(RR_validation_error)];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_RR);
    view([-45,45,45])

    %%%
    Tmax = length(dTA_IT);
    IT_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_IT_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,7)
    color_line3(X_IT(1:floor(Tper*Tmax),1),X_IT(1:floor(Tper*Tmax),2),X_IT(1:floor(Tper*Tmax),3),dTA_IT(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'STLS'}, {['$e_{ref}=$',num2str(IT_stats(1:n,2))]}, {['$e_{val}=$',num2str(IT_validation_error)]} ];
    %title_string = ['STLS,  $e_{val}=$',num2str(IT_validation_error)];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_IT);
    view([-45,45,45])


    %T = [0,TF]; 
    Tmax = length(dTA_SP);
    SP_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_SP_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(2,4,8)
    color_line3(X_SP(1:floor(Tper*Tmax),1),X_SP(1:floor(Tper*Tmax),2),X_SP(1:floor(Tper*Tmax),3),dTA_SP(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'SP'}, {['$e_{ref}=$',num2str(SP_stats(1:n,2))]}, {['$e_{val}=$',num2str(SP_validation_error)]}];
    %title_string = ['SP,  $e_{val}=$',num2str(SP_validation_error)];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_SP);
    view([-45,45,45])

    else 

    %%%

    
    Tmax = length(dTA_SPGL1);
    SPGL1_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_SPGL1_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(1,4,1)
    color_line3(X_SPGL1(1:floor(Tper*Tmax),1),X_SPGL1(1:floor(Tper*Tmax),2),X_SPGL1(1:floor(Tper*Tmax),3),dTA_SPGL1(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'BPDN'}, {['$e_{ref}=$',num2str(SPGL1_stats(1:n,2))]}, {['$e_{val}=$',num2str(SPGL1_validation_error)]}];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_SPGL1);
    view([-45,45,45])

    %T = [0,TF]; 
    
    Tmax = length(dTA_RR);
    RR_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_RR_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(1,4,2)
    color_line3(X_RR(1:floor(Tper*Tmax),1),X_RR(1:floor(Tper*Tmax),2),X_RR(1:floor(Tper*Tmax),3),dTA_RR(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'STRidge'}, {['$e_{ref}=$',num2str(RR_stats(1:n,2))]}, {['$e_{val}=$',num2str(RR_validation_error)]} ];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_RR);
    view([-45,45,45])


    %T = [0,TF]; 

    Tmax = length(dTA_IT);
    IT_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_IT_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(1,4,3)
    color_line3(X_IT(1:floor(Tper*Tmax),1),X_IT(1:floor(Tper*Tmax),2),X_IT(1:floor(Tper*Tmax),3),dTA_IT(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'STLS'}, {['$e_{ref}=$',num2str(IT_stats(1:n,2))]}, {['$e_{val}=$',num2str(IT_validation_error)]}];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_IT);
    view([-45,45,45])

    %%%
    Tmax = length(dTA_SP);
    SP_validation_error = norm(X_data_val(1:floor(Tper*val_max),:)-X_SP_val(1:floor(Tper*val_max),:))/norm(X_data_val(1:floor(Tper*val_max),:));
    subplot(1,4,4)
    color_line3(X_SP(1:floor(Tper*Tmax),1),X_SP(1:floor(Tper*Tmax),2),X_SP(1:floor(Tper*Tmax),3),dTA_SP(1:floor(Tper*Tmax))/max_dt,'LineWidth',4);
    
    view(27,16)
    axis equal;
    grid on ;
    title_string = [{'SP'}, {['$e_{ref}=$',num2str(SP_stats(1:n,2))]},  {['$e_{val}=$',num2str(SP_validation_error)]}];
    title(title_string,'interpreter','latex','fontsize',fontsize);
    xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
    set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',fontsize)

    set(gca,'FontSize',fontsize)
    length(TT_SP)
    view([-45,45,45])
    end
    drawnow
    F(frame_count) = getframe(fig);
    
    
end
%%
close all
pause(1)

save_name = ['120_',func2str(func),'_movie.mat']
save(save_name)

%%
%clear
%load('lorenz_movie.mat')
% hf = figure('Position', [10 10 1400 800]);
% pause(0.5)
% movie(hf,F,10,2);


%%
save_name = ['120_',func2str(func),'.gif']

movie2gif(F, save_name, 'LoopCount', Inf, 'DelayTime', 0)


