%COURTEMANCHE_FIG3_LSQ plots least squares regression for simplified Kpoly
%calculation and approximate values from courtemanche figure 3
% 
%   Approximate values from [profilin] = 2.5uM data from figure 3
% 
%   Makes two plots, one using the actual PRM locations, the other using
%   the values reported in the courtemanche paper.

K_court=[9.5, 8, 2.5]; %Kpoly from courtemanche graph
%K_court=[9.5, 8, 2, 2.5]; %Kpoly from courtemanche graph

N_court=[18, 28, 65]; %courtemanche n
%N_court=[18, 28, 37, 65]; %courtemanche n

%N_act=[25, 36, 47, 75]; %how our model defines n
N_act=[25, 36, 75]; %how our model defines n

%P_occ= [0.9185, 0.9215, 0.9139, 0.9075] %P_occ for single
P_occ= [0.9185, 0.9215, 0.9075] %P_occ for single


YData=K_court./(1-P_occ) %move (1-P_occ) to other side so that lsq fxn will work

fun_court = @(c,N_court)c(1).*(N_court.^(-3/2)./(N_court.^(-3/2)+c(2))); %using N_court

c0=[10,0]
lb=[0,-0.002]
ub=[]
options = optimoptions('lsqcurvefit')
options.StepTolerance = 1e-14;
options.FunctionTolerance = 1e-20;
options.OptimalityTolerance = 1e-14;
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1e3;

c = lsqcurvefit(fun_court,c0,N_court,YData,lb,ub,options)

figure

times = linspace(N_court(1),N_court(end));
plot(N_court,YData,'ko',times,(fun_court(c,times)),'b-')
legend('Data','Fitted curve')
title('Polymerization rate vs distance from PRM to FH2')
xlabel('n (courtemanche)');
ylabel('Kpoly/(1-P_{occ})');

subtitle1= "$K_{poly}=c1\cdot\left(1-P_{occ}\right)\cdot\left(\frac{n^{-\frac{3}{2}}}{n^{-\frac{3}{2}}+c2}\right)$"+"  c1="+c(1)+"  c2="+c(2)
hl=subtitle(subtitle1)
set(hl ,'Interpreter','latex')

fun_act = @(c,N_act)c(1).*(N_act.^(-3/2)./(N_act.^(-3/2)+c(2))); %using N_act

c = lsqcurvefit(fun_act,c0,N_act,YData,lb,ub,options)

figure

times = linspace(N_act(1),N_act(end));
plot(N_act,YData,'ko',times,(fun_court(c,times)),'b-')
legend('Data','Fitted curve')
title('Polymerization rate vs distance from PRM to FH2')
xlabel('n');
ylabel('Kpoly/(1-P_{occ})');

subtitle1= "$K_{poly}=c1\cdot\left(1-P_{occ}\right)\cdot\left(\frac{n^{-\frac{3}{2}}}{n^{-\frac{3}{2}}+c2}\right)$"+"  c1="+c(1)+"  c2="+c(2)
hl=subtitle(subtitle1)
set(hl ,'Interpreter','latex')

