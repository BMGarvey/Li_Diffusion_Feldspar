%Numerical Li D single path with two fronts to account for core composition
%raised
clear all

%Diffusivity max and min from separate fitting
Dslow=3e-14;
Dfast=3e-12;

%Loads in data and assigns columns
diff_data=xlsread('LiAn-1000.xlsx');
obs_Li=diff_data(:,1);
obs_pos=diff_data(:,2);
%convert distance to m
obs_pos_m=obs_pos*1e-6;

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
dt=.1;
%Original concentration in feldspar
Core=8;%diff_data(1,10);

%Make array for full transect from mirrored half transect
full=min_index*2;
test=ones(full,1)*Core;
test(1)=C_left;
test(end)=C_left;

%Array to store slow model runs
D_Slow=zeros(full,100);
number=0;

%Matrix to track the combos of D values
[n,m]=ndgrid(Dfast,Dslow);
Combos=[m(:),n(:)];

% overall
for k=1:length(Dslow)
    tracking=0;
    for h=1:length(Dfast)
        number=number+1
        num_cols=time/dt;
        new_column=zeros(full,num_cols+1);
        new_column(:,1)=test;
        new_column(1,:)=C_left;
        new_column(end,:)=C_left;
        delta=log10(Dfast(h))-log10(Dslow(k));
        for j=1:num_cols
            % first equation
            for i = 2:full-1
            c_before=new_column(i-1, j);
            c_current=new_column(i, j);
            c_after=new_column(i+1, j);
            %Calculate D step
            D_current=10^(log10(Dfast(h))-(delta*((c_current-Core)/(C_left-Core))));
            D_after=10^(log10(Dfast(h))-(delta*((c_after-Core)/(C_left-Core))));
            r=(dt*D_current)/(dx^2);
            term1=r*((c_after-(2*c_current)+c_before));
            term2=dt*((D_after-D_current)/dx)*((c_after-c_current)/dx);
            new_column(i,j+1)=c_current+term1+term2;
            end
        end
    model_S=(new_column(:,end));
    D_Slow(:,number)=model_S(:);
    end
end

%Cut the model to compare to observed data with lsq in next step
D_Slow(min_index+1:end,:)=[];
D_Slow=flip(D_Slow);

%Simple NL regression
NLR_Slow=zeros(size(D_Slow,1), size(D_Slow,2));
for p=1:size(D_Slow,2)
        NLR=(((obs_Li-D_Slow(:,p))./obs_Li).^2);
        NLR_Slow(:,p)=NLR(:);
end
NLR_Slow_Sum=sum(NLR_Slow, 1, 'omitmissing');

%Now select the best fit based on the lowest residual sums
[Slow_Res, Slow_index]=min(NLR_Slow_Sum);

%Select the D values
Best=Combos(Slow_index,:);
D1=Best(1,1);
D2=Best(1,2);

%Plotting
figure(1)
plot(obs_pos, D_Slow, 'm-')
hold on
plot(obs_pos, obs_Li, 'ko')
xlabel('Distance from rim (um)')
ylabel('Li (ppm)')
title('Model runs')
hold off

%Plotting
figure(2)
plot(obs_pos, D_Slow(:,Slow_index), 'm-')
hold on
yline(Core, 'k--')
plot(obs_pos, obs_Li, 'ko')
legend('Diffusion Model', '', 'FontSize', 11)
ylabel('Li (ppm)', 'FontSize', 14)
xlabel('Distance (um)', 'FontSize', 14)
txt1=['D Slow: ' num2str(D1) '  |  ' 'D Fast: ' num2str(D2)];
txt1=join(txt1);
annotation('textbox', [0.15 0.7 0.5 0.2], 'String', txt1, 'FitBoxToText','on', 'BackgroundColor', 'white', 'FontSize', 10)
title('Concentration dependent diffusion model', 'FontSize', 16)
hold off
