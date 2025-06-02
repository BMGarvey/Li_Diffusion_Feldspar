%Numerical model to calculate timescales for single path diffusion

%This code is for a diffusion profile that looks to have 1 path of
%diffusion (i.e. single erf look) but first uses the slow parameters then
%uses the fast parameters to get a bracketed timescale.

clear all

%Your slow path arhhenius parameters
S_Do=10^-5.3;
S_Ea=178130;
T=827+273.14;

%Your fast path arhhenius parameters
F_Do=10^-6.0;
F_Ea=144050;

%Calculates D for both paths
S_D=S_Do*exp(-S_Ea/(8.314*T));
F_D=F_Do*exp(-F_Ea/(8.314*T));

%Loads in data and assigns columns
diff_data=xlsread('LiAn-BOBB1.xlsx','Sheet1');
obs_Li=diff_data(:,1);
obs_pos=diff_data(:,2);
%convert distance to m
obs_pos_m=obs_pos*1e-6;

%Visualize the profile
figure(1)
plot(obs_pos,obs_Li, 'ok')
ylabel('Li (ppm)')
xlabel('Distance (um)')

%Setting Variables

%dt in sec (can change this to satisfy stablility of model. Smaller dt=longer
%run time though
dt=2;
%duration, seconds
time=0:dt:518400;

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

%Make array to store all time iterations (slow and fast)
t_Slow=zeros(length(obs_pos),length(time));

%Make array to store residuals of different time lines
NLR_Slow=zeros(length(obs_pos),length(time));

%Setting initial conditions
Initial=ones(min_index,1)*C_right;
Initial(1)=C_left;

%First numerical model for the slow path + plotting it

%r must be less than 0.5
r=(dt*S_D)/(dx^2);
num_cols=length(time);
new_column_S=zeros(min_index,num_cols);
new_column_S(:,1)=Initial;
new_column_S(1,:)=C_left;
new_column_S(end,:)=C_right;
    for j=1:num_cols
    % first equation
        for i = 2:min_index-1
        c_before=new_column_S(i-1, j);
        c_current=new_column_S(i, j);
        c_after=new_column_S(i+1, j);
        new_column_S(i,j+1)=c_current+r*((c_after-(2*c_current)+c_before));
        end
    end
new_column_S=flip(new_column_S);

%Simple NL regression
for p=1:size(new_column_S,2)
        NLR=(((obs_Li-new_column_S(:,p))./obs_Li).^2);
        NLR_Slow(:,p)=NLR(:);
end
NLR_Slow_Sum=sum(NLR_Slow, 1, 'omitmissing');

%Now Select the proper time value based on the lowest residual sum
[Slow_Res, Slow_index]=min(NLR_Slow_Sum);

%Plotting
figure(2)
plot(obs_pos, new_column_S(:,Slow_index), 'm-')
hold on
plot(obs_pos, obs_Li, 'ko')
legend('Slow path', '')
txt1=['Time: ' num2str(time(Slow_index))];
txt1=join(txt1);
annotation('textbox', [0.15 0.7 0.5 0.2], 'String', txt1, 'FitBoxToText','on', 'BackgroundColor', 'white', 'FontSize', 10)
ylabel('Li (ppm)')
xlabel('Distance (um)')
title('Slow path diffusion model')
hold off

%Second numerical model for fast path + plotting it

%r must be less than 0.5
r_f=(dt*F_D)/(dx^2);
num_cols=length(time);
new_column_F=zeros(min_index,num_cols);
new_column_F(:,1)=Initial;
new_column_F(1,:)=C_left;
new_column_F(end,:)=C_right;
    for j=1:num_cols
    % first equation
        for i = 2:min_index-1
        c_before=new_column_F(i-1, j);
        c_current=new_column_F(i, j);
        c_after=new_column_F(i+1, j);
        new_column_F(i,j+1)=c_current+r_f*((c_after-(2*c_current)+c_before));
        end
    end
new_column_F=flip(new_column_F);

%Simple NL regression
for p=1:size(new_column_F,2)
        NLR=(((obs_Li-new_column_F(:,p))./obs_Li).^2);
        NLR_Fast(:,p)=NLR(:);
end
NLR_Fast_Sum=sum(NLR_Fast, 1, 'omitmissing');

%Now Select the proper time value based on the lowest residual sum
[Fast_Res, Fast_index]=min(NLR_Fast_Sum);

%Plotting
figure(3)
plot(obs_pos, new_column_F(:,Fast_index), 'm-')
hold on
plot(obs_pos, obs_Li, 'ko')
legend('Fast path', '')
txt1=['Time: ' num2str(time(Fast_index))];
txt1=join(txt1);
annotation('textbox', [0.15 0.7 0.5 0.2], 'String', txt1, 'FitBoxToText','on', 'BackgroundColor', 'white', 'FontSize', 10)
ylabel('Li (ppm)')
xlabel('Distance (um)')
title('Fast path diffusion model')
hold off

%Initializing outputs for Slow and Fast times
Slow_Time=time(Slow_index);
Fast_Time=time(Fast_index);

Avg_Time=(Slow_Time+Fast_Time)/2;