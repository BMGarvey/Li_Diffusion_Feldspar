%Numerical Li D single path with two fronts to account for core composition
%raised
clear all

D1=2.49e-14;
D2=1.00e-12;
delta=D2-D1;

%Binary diffusion things
DNa=0.000794*exp(-268000/(8.314*1373));%Start with a lit. number then go from there (i think this is an60)
DNa=1e-14;%10^-13.02;
D2=4e-9;%10^-10.49;

%Loads in data and assigns columns
diff_data=xlsread('LiAn-1000.xlsx');
obs_Li=diff_data(:,1);
obs_pos=diff_data(:,2);
%convert distance to m
obs_pos_m=obs_pos*1e-6;

%Converting profile to mol % for effective binary diffusion
XNa=0.04568; %An60 0.3867 %Ab 0.9851   %An95 0.04568
XCa=0.9540; %An60 0.6133 %Ab 0.004409 %An95 0.9540
Li_wt=obs_Li(:,1)./10000;
Li_mol=(Li_wt(:,1)./6.94)./((Li_wt(:,1)./6.94)+XNa+XCa);

%Setting Variables
%duration, seconds
time=diff_data(1,4);
%Boundary conditions for first fit (Slow path)
C_left=diff_data(1,8);
%Step size
dx=obs_pos_m(1)-obs_pos_m(2);
%Profile length
profile_length=max(obs_pos_m);
%Interface location
[interface, min_index]=min(obs_pos_m);
%dt in sec
dt=.01;
%Original concentration in feldspar
Core=8;%diff_data(1,10); 50

%convert core and C_left to mol
Core_wt=Core/10000;
Core_mol=(Core_wt./6.94)./((Core_wt(:,1)./6.94)+XNa+XCa);
Left_wt=C_left/10000;
Left_mol=(Left_wt./6.94)./((Left_wt(:,1)./6.94)+XNa+XCa);

%Functional form concentration dependence

%Make array for full transect from mirrored half transect
full=min_index*2;
test=ones(full,1)*Core;
test(1)=C_left;
test(end)=C_left;

%Array to store slow model runs
D_Slow=zeros(full,1);
number=0;

% overall
        number=number+1
        num_cols=time/dt;
        new_column=zeros(full,num_cols+1);
        new_column(:,1)=test;
        new_column(1,:)=C_left;
        new_column(end,:)=C_left;
        for j=1:num_cols
            % first equation
            for i = 2:full-1
            c_before=new_column(i-1, j);
            c_current=new_column(i, j);
            c_after=new_column(i+1, j);
            %Calculate D step
            %Linear
            % D_current=10^((-0.028*c_current)-11.386);
            % D_after=10^((-0.028*c_after)-11.386);
            %Exponential
            D_current=10^(-11.4*exp(0.002*c_current));
            D_after=10^(-11.4*exp(0.002*c_after));
            % D_current=10^(-8.2589*exp(0.0031*c_current)+2.9975*exp(-0.0910*c_current));
            % D_after=10^(-8.2589*exp(0.0031*c_after)+2.9975*exp(-0.0910*c_after));
            % D_current=D2-(delta*exp(-0.01*c_current))
            % D_after=D2-(delta*exp(-0.01*c_after))
            %Effective Binary Diffusion
            % D_current=(XNa*D2*DNa)/((c_current*D2)+(XNa*DNa))
            % D_after=(XNa*D2*DNa)/((c_after*D2)+(XNa*DNa))
            r=(dt*D_current)/(dx^2);
            term1=r*((c_after-(2*c_current)+c_before));
            term2=dt*((D_after-D_current)/dx)*((c_after-c_current)/dx);
            new_column(i,j+1)=c_current+term1+term2;
            end
        end
    model_S=(new_column(:,end));
    D_Slow(:,number)=model_S(:);

%Cut the model to compare to observed data with lsq in next step
D_Slow(min_index+1:end,:)=[];
D_Slow=flip(D_Slow);

% %Simple NL regression
% NLR_Slow=zeros(size(D_Slow,1), size(D_Slow,2));
% for p=1:size(D_Slow,2)
%         NLR=(((obs_Li-D_Slow(:,p))./obs_Li).^2);
%         NLR_Slow(:,p)=NLR(:);
% end
% NLR_Slow_Sum=sum(NLR_Slow, 1, 'omitmissing');

% %Now select the best fit based on the lowest residual sums
% [Slow_Res, Slow_index]=min(NLR_Slow_Sum);

% %Select the D values
% Best=Combos(Slow_index,:);
% D1=Best(1,1);
% D2=Best(1,2);

%Plotting
% figure(1)
% plot(obs_pos, D_Slow, 'm-')
% hold on
% plot(obs_pos, obs_Li, 'ko')
% xlabel('Distance from rim (um)')
% ylabel('Li (ppm)')
% title('Model runs')
% hold off

%Plotting
figure(2)
plot(obs_pos, D_Slow, 'm-')
hold on
yline(Core, 'k--')
plot(obs_pos, obs_Li, 'ko')
legend('Diffusion Model', '', 'FontSize', 11)
ylabel('Li (ppm)', 'FontSize', 14)
xlabel('Distance (um)', 'FontSize', 14)
title('Concentration dependent diffusion model', 'FontSize', 16)
hold off

%Make array for full transect from mirrored half transect
full=min_index*2;
test=ones(full,1)*Core_mol;
test(1)=Left_mol;
test(end)=Left_mol;

%Array to store slow model runs
D_Slow=zeros(full,1);
number=0;

% overall
bar=waitbar(0, 'Starting');
        number=number+1;
        num_cols=time/dt;
        new_column=zeros(full,num_cols+1);
        new_column(:,1)=test;
        new_column(1,:)=Left_mol;
        new_column(end,:)=Left_mol;
        for j=1:num_cols
            % first equation
            for i = 2:full-1
            c_before=new_column(i-1, j);
            c_current=new_column(i, j);
            c_after=new_column(i+1, j);
            %Calculate D step
            %Linear
            % D_current=10^((-0.028*c_current)-11.386)
            % D_after=10^((-0.028*c_after)-11.386)
            %Exponential
            % D_current=10^(-11.4*exp(0.002*c_current))
            % D_after=10^(-11.4*exp(0.002*c_after))
            % D_current=10^(-8.2589*exp(0.0031*c_current)+2.9975*exp(-0.0910*c_current))
            % D_after=10^(-8.2589*exp(0.0031*c_after)+2.9975*exp(-0.0910*c_after))
            % D_current=D2-(delta*exp(-0.01*c_current))
            % D_after=D2-(delta*exp(-0.01*c_after))
            %Effective Binary Diffusion
            D_current=(XNa*D2*DNa)/((c_current*D2)+(XNa*DNa));
            D_after=(XNa*D2*DNa)/((c_after*D2)+(XNa*DNa));
            r=(dt*D_current)/(dx^2);
            term1=r*((c_after-(2*c_current)+c_before));
            term2=dt*((D_after-D_current)/dx)*((c_after-c_current)/dx);
            new_column(i,j+1)=c_current+term1+term2;
            end
            waitbar(j/num_cols, bar, sprintf('We runnin: %d %%', floor(j/num_cols*100)));
        end
    model_S=(new_column(:,end));
    D_Slow(:,number)=model_S(:);

%Cut the model to compare to observed data with lsq in next step
D_Slow(min_index+1:end,:)=[];
D_Slow=flip(D_Slow);

% %Simple NL regression
% NLR_Slow=zeros(size(D_Slow,1), size(D_Slow,2));
% for p=1:size(D_Slow,2)
%         NLR=(((obs_Li-D_Slow(:,p))./obs_Li).^2);
%         NLR_Slow(:,p)=NLR(:);
% end
% NLR_Slow_Sum=sum(NLR_Slow, 1, 'omitmissing');

% %Now select the best fit based on the lowest residual sums
% [Slow_Res, Slow_index]=min(NLR_Slow_Sum);

% %Select the D values
% Best=Combos(Slow_index,:);
% D1=Best(1,1);
% D2=Best(1,2);

%Plotting
% figure(1)
% plot(obs_pos, D_Slow, 'm-')
% hold on
% plot(obs_pos, obs_Li, 'ko')
% xlabel('Distance from rim (um)')
% ylabel('Li (ppm)')
% title('Model runs')
% hold off

%Plotting
figure(3)
plot(obs_pos, D_Slow, 'm-')
hold on
yline(Core_mol, 'k--')
plot(obs_pos, Li_mol, 'ko')
legend('Diffusion Model', '', 'FontSize', 11)
ylabel('Li (mol %)', 'FontSize', 14)
xlabel('Distance (um)', 'FontSize', 14)
title('Effective binary diffusion model', 'FontSize', 16)
hold off

close(bar)