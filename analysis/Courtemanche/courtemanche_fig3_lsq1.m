%%Get values
pdf_name = "Courtemanche LSQ plots.pdf"

vals=readtable("Courtemanche plot digitizer data F3.csv");


KPoly=[vals.K_18,vals.K_28,vals.K_37,vals.K_65];
Errors=[vals.ErrPer_18,vals.ErrPer_28,vals.ErrPer_37,vals.ErrPer_65];
Nlabels=["n (courtemanche)","n"];


%K_court=[9.5, 8, 2.5]; %Kpoly from courtemanche graph
%K_court=[9.5, 8, 2, 2.5]; %Kpoly from courtemanche graph

%N_court=[18, 28, 65]; %courtemanche n
N_court=[18, 28, 37, 65]; %courtemanche n

N_act=[25, 36, 47, 75]; %how our model defines n
%N_act=[25, 36, 75]; %how our model defines n


P_occ= [0.9185, 0.9215, 0.9139, 0.9075]; %P_occ for single
%P_occ= [0.9185, 0.9215, 0.9075] %P_occ for single


%% LSQ options
c0=[10,0];
lb=[0,-0.002];
ub=[];
options = optimoptions('lsqcurvefit');
options.StepTolerance = 1e-14;
options.FunctionTolerance = 1e-20;
options.OptimalityTolerance = 1e-14;
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1e3;

%% initialize arrays
C1_ncourt=[];
C1_n=[];
C2_ncourt=[];
C2_n=[];

%% run LSQ for each
for LOOP=1:length(vals.x_profilin_UM) %for each [profilin]
CP=vals.x_profilin_UM(LOOP);
K=KPoly(LOOP,:);
Err=Errors(LOOP,:); % percent errors from table
YData=K./(1-P_occ); %move (1-P_occ) to other side so that lsq fxn will work
AbsErr=((Err.*(K+5))./(1-P_occ))/2 %the absolute error on one side (for plotting)

NMat=[N_court;N_act];

for i=1:length(K)
    if isnan(K(i))
        NMat(:,i)=[];
        Err(i)=[];
        YData(i)=[];
        AbsErr(i)=[];
    end
end
K=K(~isnan(K));

C1=[];
C2=[];
            
for NLOOP=1:2 %for our N and court N
    N=NMat(NLOOP,:);
    
    fun = @(c,N)c(1).*(N.^(-3/2)./(N.^(-3/2)+c(2))); 
    
    c = lsqcurvefit(fun,c0,N,YData,lb,ub,options)
    
    C1=[C1,c(1)];
    C2=[C2,c(2)];

    figure
    
    errorbar(N,YData,AbsErr,"o")

    hold on
    times = linspace(N(1)-2,N(end)+2);
    plot(N,YData,'ko',times,(fun(c,times)),'b-')
    legend('Data','Fitted curve')
    title('Polymerization rate vs distance from PRM to FH2')
    xlabel(Nlabels(NLOOP));
    ylabel('Kpoly/(1-P_{occ})');

    subtitle1= "[P]="+CP+" $K_{poly}=c1\cdot\left(1-P_{occ}\right)\cdot\left(\frac{n^{-\frac{3}{2}}}{n^{-\frac{3}{2}}+c2}\right)$"+"  c1="+c(1)+"  c2="+c(2)
    hl=subtitle(subtitle1)
    set(hl ,'Interpreter','latex')
    
    saveas(gcf, append('temp.pdf'))
    append_pdfs(pdf_name, append('temp.pdf'))
end

C1_ncourt=[C1_ncourt;C1(1)];
C1_n=[C1_n;C1(2)];
C2_ncourt=[C2_ncourt;C2(1)];
C2_n=[C2_n;C2(2)];

end

T = table(C1_ncourt,C2_ncourt,C1_n,C2_n);

T2=[vals(:,1) T vals(:,[2:width(vals)])];

delete 'temp.pdf'

