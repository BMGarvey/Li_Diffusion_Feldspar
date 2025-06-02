%Diffy Boi function to numerically model diffusion coefficients

function [Slow_DC, Fast_DC, Slow_Res, Fast_Res, Fast_Co] = Diff_2Front_NoPlot(excel_file, excel_sheet, Diffs, FC_Left, Inflection, Core)

    %Loads in data and assigns columns
    diff_data=xlsread(excel_file,excel_sheet);
    obs_Li=diff_data(:,1);
    obs_pos=diff_data(:,2);
    %convert distance to m
    obs_pos_m=obs_pos*1e-6;
    
    %Your input parameter guesses
    %Number of slow and fast guesses need to match for statistics to work
    Diff=Diffs;%Slow path guess
    F_Diff=Diffs; %Fast path guess
    F_C_left=FC_Left; %Fast path boundary condition (y-int)
    inflection=Inflection; %Distance (in um) that separates fast and slow regimes
    core=Core;
    
    %Setting Variables
    %duration, seconds
    time=diff_data(1,4);
    %Boundary conditions for first fit (Slow path)
    C_left=diff_data(1,8);
    C_right=diff_data(1,9);
    %Step size
    dx=obs_pos_m(1)-obs_pos_m(2);
    %Profile length
    profile_length=max(obs_pos_m);
    %Number of steps
    steps=length(obs_pos_m)-1;
    %Interface location
    [interface, min_index]=min(obs_pos_m);
    %dt in sec
    dt=.01;
    
    %Moving average to find inflection point
    Profile_Avg=movmean(obs_Li, 5, 'omitnan');
    % inflection=ischange(Profile_Avg, 'linear', MaxNumChanges=1);
    point=interp1(obs_pos, 1:length(obs_pos), inflection, 'nearest');
    
    % figure(1)
    % plot(obs_pos,obs_Li, 'ok')
    % hold on
    % plot(obs_pos,Profile_Avg, '-m')
    % xline(obs_pos(point)+2, 'k--')
    % hold off
    
    %Extra Nans to add onto the Slow and Fast Li for correct matrix math
    Slow_Nan=NaN(point-1,1);
    Fast_Nan=NaN(length(obs_Li)-point,1);
    
    %Split the data in slow and fast sections according to inflection point
    Slow_Li=obs_Li(point:end);
    Fast_Li=obs_Li(1:point);
    
    %Add in the Nans
    Slow_Li=vertcat(Slow_Nan,Slow_Li);
    Fast_Li=vertcat(Fast_Li,Fast_Nan);
    
    %Make array to store slow and fast portions of model (keep for now)
    % Slow_Model=zeros(length(Slow_Li), length(Diff));
    % Fast_Model=zeros(length(Fast_Li), length(F_Diff));
    
    %Make array to store residuals of different D lines
    NLR_Slow=zeros(length(obs_pos),length(Diff));
    NLR_Fast=zeros(length(obs_pos),length(F_Diff));

    %Make array for full transect from mirrored half transect
    full=min_index*2;
    test=ones(full,1)*core;
    test(1)=C_left;
    test(end)=C_left;

    %Array to store slow model runs
    D_Slow=zeros(full,length(Diff));
    
    % overall
    for q=1:length(Diff)
    %r must be less than 0.5
    r=(dt*Diff(q))/(dx^2);
    num_cols=time/dt;
    new_column=zeros(full,num_cols+2);
    new_column(:,1)=test;
    new_column(1,:)=C_left;
    new_column(end,:)=C_left;
        for j=1:num_cols+1
        % first equation
            for i = 2:full-1
            c_before=new_column(i-1, j);
            c_current=new_column(i, j);
            c_after=new_column(i+1, j);
            new_column(i,j+1) = c_current + r*((c_after-(2*c_current)+c_before));
            end
        end
    model_S=flip(new_column(:,end));
    D_Slow(:,q)=model_S(:);
    end
    
    % %Simple NL regression
    % for p=1:size(D_Slow,2)
    %         NLR=(((Slow_Li-D_Slow(:,p))./Slow_Li).^2);
    %         NLR_Slow(:,p)=NLR(:);
    % end
    % NLR_Slow_Sum=sum(NLR_Slow, 1, 'omitmissing');
    %Run again for the fast path diffusion
    
    %Find point where fast path is dominant
    %Dist_ind=interp1(obs_pos_m,1:length(obs_pos_m),S_Dist,'nearest');
    % Dist_ind=point;
    %Fast path x-axis
    % Fast_pos=obs_pos_m(1:Dist_ind);
    
    %Slow path x-axis
    % Slow_pos=obs_pos_m(Dist_ind:end);
    
    %Data point where the slow-fast transition is (used for plotting)
    % Change_point=(Fast_pos(end)*1e6)+2;
    
    %Display the cutoff point
    % figure(2)
    % plot(Fast_pos*1e6,obs_Li(1:Dist_ind), 'ro')
    % hold on
    % plot(Slow_pos*1e6,obs_Li(Dist_ind:end), 'bo')
    % xline(obs_pos(point)+2, 'k--')
    % title('Selected points for slow and fast path modelling')
    % legend('Fast path', 'Slow path', '')
    % ylabel('Li (ppm)')
    % xlabel('Distance (um)')
    % hold off
    
    %Make array for full transect from mirrored half transect
    test2=ones(full,1)*Core;
    test2(1)=F_C_left;
    test2(end)=F_C_left;

    %Array to store fast model runs
    D_Fast=zeros(full,length(F_Diff));

    %Fast path calculation
    for qq=1:length(F_Diff)
    r_f=(dt*F_Diff(qq))/(dx^2);
    num_cols=time/dt;
    new_column_F=zeros(full,num_cols+2);
    new_column_F(:,1)=test2;
    new_column_F(1,:)=F_C_left;
    new_column_F(end,:)=F_C_left;
    
        for j = 1:num_cols+1
        % second equation
        %for i = 2:Dist_ind-1
            for i = 2:full-1
            c_before=new_column_F(i-1, j);
            c_current=new_column_F(i, j);
            c_after=new_column_F(i+1, j);
            new_column_F(i,j+1)=c_current+r_f*((c_after-(2*c_current)+c_before));
            end
        end
    model_F=flip(new_column_F(:,end));
    D_Fast(:,qq)=model_F(:);
    end
   
    %Cut the slow model to compare to observed data with lsq in next step
    D_Slow(min_index+1:end,:) = [];
    D_Slow=flip(D_Slow);

    %Cut the fast model to compare to observed data with lsq in next step
    D_Fast(min_index+1:end,:) = [];
    D_Fast=flip(D_Fast);

    %Simple lsq for slow path
    for p=1:size(D_Slow,2)
            NLR=(((Slow_Li-D_Slow(:,p))./Slow_Li).^2);
            NLR_Slow(:,p)=NLR(:);
    end
    NLR_Slow_Sum=sum(NLR_Slow, 1, 'omitmissing');

    %Simple lsq for fast path
    for p=1:size(D_Fast,2)
            NLR=(((Fast_Li-D_Fast(:,p))./Fast_Li).^2);
            NLR_Fast(:,p)=NLR(:);
    end
    NLR_Fast_Sum=sum(NLR_Fast, 1, 'omitmissing');
    
    %Now Select the proper D values based on the lowest residual sums
    [Slow_Res, Slow_index]=min(NLR_Slow_Sum);
    [Fast_Res, Fast_index]=min(NLR_Fast_Sum);
    
    %Plotting
    % figure(3)
    % plot(obs_pos, D_Slow, 'b-')
    % hold on
    % plot(obs_pos, D_Fast, 'r-')
    % plot(obs_pos, obs_Li, 'ko')
    % xline(obs_pos(point)+2, 'k--')
    % xlabel('Distance (um)')
    % ylabel('Li (ppm)')
    % title('2 path finite difference forward model runs')
    % %legend('Slow model runs','','','','','','','','','','','','','','','','','','','','','','','','','',"Fast model runs",'','','','','','','','','','','','','','','','','','','','','','','','', 'Slow-Fast regime change')
    % hold off
    % 
    %Plotting
    % figure(4)
    % plot(obs_pos, D_Slow(:,Slow_index), 'b-')
    % hold on
    % plot(obs_pos, D_Fast(:,Fast_index), 'r-')
    % xline(obs_pos(point)+2, 'k--')
    % plot(obs_pos, obs_Li, 'ko')
    % % errorbar(obs_pos, obs_Li, Two_SE, 'vertical', 'k')
    % legend('Slow path', 'Fast path', '')
    % txt1=['D Slow: ' num2str(Diff(Slow_index)) '  |  ' 'D Fast: ' num2str(F_Diff(Fast_index))];
    % txt1=join(txt1);
    % annotation('textbox', [0.15 0.7 0.5 0.2], 'String', txt1, 'FitBoxToText','on', 'BackgroundColor', 'white', 'FontSize', 10)
    % ylabel('Li (ppm)')
    % xlabel('Distance (um)')
    % title('2 path diffusion model')
    % hold off
    
    % col_slow=new_column(:,end);
    % col_fast=new_column_F(:,end);
    
    % Initializing outputs for Slow and Fast diffusion coefficients
    Slow_DC=Diff(:,Slow_index);
    Fast_DC=F_Diff(:,Fast_index);

    %Pulling the Y-Int associated with the fast path run
    Fast_Path=D_Fast(:,Fast_index);
    Fast_Co=Fast_Path(end);

    % figure(5)
    % subplot(2,1,1)
    % plot(Fast_pos*1e6,obs_Li(1:Dist_ind), 'ro')
    % hold on
    % plot(Slow_pos*1e6,obs_Li(Dist_ind:end), 'bo')
    % title('Selected points for slow and fast path modelling')
    % legend('Fast path', 'Slow path')
    % ylabel('Li (ppm)')
    % xlabel('Distance (um)')
    % hold off
    % subplot(2,1,2)
    % plot(obs_pos, D_Slow(:,Slow_index), 'b' )
    % hold on
    % %plot(obs_pos,flip(right_length(:,end)), 'r')
    % plot(obs_pos, D_Fast(:,Fast_index), 'r')
    % plot(obs_pos, obs_Li, 'ko')
    % legend('Slow path', 'Fast path')
    % ylabel('Li (ppm)')
    % xlabel('Distance (um)')
    % title('2 path diffusion model')
    % hold off
    
end 