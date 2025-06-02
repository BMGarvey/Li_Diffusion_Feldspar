%Testing Monte Carlo diffusion code

%Best guess of the y-int of the fast path (Co). The model will then run
%through plus and minus 5 ppm around this point.
Yint=17;
FC_Left=Yint-5:1:Yint+5;

%Best guess of the inflection point that separates the slow and fast
%regimes. This is mostly used to get the best fit line.
Inflection=38; %In microns

%Core concentration of your crystal in ppm
Core=8;

%Number of iterations to run the Monte Carlo simulation
iter=100;

%Need to initialize space for all the fast path y-int run
Yint_runs=zeros(length(FC_Left),6);

%Need to initialize space for all the runs of the Monte Carlo
D_MC=zeros(iter,4);

count=0;
bar=waitbar(0, 'Starting');
for z=1:iter
    count=count+1;   
    %tic

%First need to initialize the random diffusivities
    D_exp=-1*randi([10 15],1,200);
    D_float=1+(10-1).*rand(1,200);

    Diffs=zeros(1,150);

    for b=1:length(Diffs)
        Diffs(1,b)=D_float(1,b).*10.^D_exp(1,b);
    end
    % FC_Left=Yint-5:1:Yint+5;
    % FC_Left=20:10:50;
count2=0;
    for p=FC_Left(1):1:FC_Left(end)
        count2=count2+1;
        [Slow_DC, Fast_DC, Slow_Res, Fast_Res, Fast_Co]=Diff_2Front_NoPlot('LiAn-1000C.xlsx', 'Sheet1', Diffs, p, Inflection, Core);

        %Store the y-int results
        Yint_runs(count2,1)=Slow_DC;
        Yint_runs(count2,2)=Slow_Res;
        Yint_runs(count2,3)=Fast_DC;
        Yint_runs(count2,4)=Fast_Res;
        Yint_runs(count2,5)=Slow_Res+Fast_Res;
        Yint_runs(count2,6)=Fast_Co;

        %Now select the best fit from the range of y-ints given
        % [Y_Res, Y_index]=min(Yint_runs(:,5));
        % 
        % %Take everything from the best fit and store it
        % D_MC(count,1)=count;
        % D_MC(count,2)=Yint_runs(Y_index,1); %Ds
        % D_MC(count,3)=Yint_runs(Y_index,3); %Df
        % D_MC(count,4)=Yint_runs(Y_index,6); %Fast path Co used
    end
        %Now select the best fit from the range of y-ints given
        [Y_Res, Y_index]=min(Yint_runs(:,5));

        %Take everything from the best fit and store it
        D_MC(count,1)=count;
        D_MC(count,2)=Yint_runs(Y_index,1); %Ds
        D_MC(count,3)=Yint_runs(Y_index,3); %Df
        D_MC(count,4)=Yint_runs(Y_index,6); %Fast path Co used
    %toc
     waitbar(count/iter, bar, sprintf('We runnin: %d %%', floor(count/iter*100)));
end

figure(1)
subplot(1,2,1)
histogram(D_MC(:,2), 'facecolor', 'b')
hold on
xlabel('Diffusivity')
ylabel('Frequency')
title('Slow Path Diffusivities')
hold off
subplot(1,2,2)
histogram(D_MC(:,3), 'facecolor', 'r')
hold on
xlabel('Diffusivity')
ylabel('Frequency')
title('Fast Path Diffusivities')
hold off

figure(2)
histogram(D_MC(:,4),'facecolor', 'g')
hold on
xlabel('Fast path C_O')
ylabel('Frequency')
title('Fast Path C_O Distribution')

%Confidence intervals
Ds=mean(D_MC(:,2));
SEMs=std(D_MC(:,2))/sqrt(length(D_MC(:,2)));
TSs=tinv([0.025 0.975],length(D_MC(:,2))-1);
CIs=Ds+TSs*SEMs;
Uncertain_Slow=Ds-CIs;

Df=mean(D_MC(:,3));
SEMf=std(D_MC(:,3))/sqrt(length(D_MC(:,3)));
TSf=tinv([0.025 0.975],length(D_MC(:,3))-1);
CIf=Df+TSf*SEMf;
Uncertain_Fast=Df-CIf;

%Just to visualize the distribution of fast Co, I dont think the statistics
%are that reliable
FY=mean(D_MC(:,4));
SEMy=std(D_MC(:,4))/sqrt(length(D_MC(:,4)));
TSy=tinv([0.025 0.975],length(D_MC(:,4))-1);
CIy=Fast_Co+TSy*SEMy;
Uncertain_FastCo=FY-CIy;

%Find the best match fast path Co to the D_fast mean
nearestfast=interp1(D_MC(:,3), 1:length(D_MC(:,3)), Df, 'nearest');
Fast_Co=D_MC(nearestfast,4);

%Now Run the model once with the best fit coefficients to get figures
[Slow_DC, Fast_DC]=Diff_2Front('LiAn-1000C.xlsx', 'Sheet1', Ds, Df, Fast_Co, Inflection, Core);

