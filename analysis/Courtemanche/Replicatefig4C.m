%REPLICATEFIG4C reads in values from csv file and uses them to replicate
%courtemanche figure 4c scatter plot with error bars
% 
%   Reads in data from "Courtemanche plot digitizer data F4C.csv"
%
%   See also REPLICATEFIG3, REPLICATEFIG4A.

values=readmatrix("Courtemanche plot digitizer data F4C.csv");
PConc=values(:,1);
K_B=values(:,2);
Er_B=values(:,3);
K_C=values(:,6);
Er_C=values(:,7);
K_D=values(:,10);
Er_D=values(:,11);

hold off
errorbar(PConc,K_B,(Er_B.*(K_B+5))/2,"MarkerFaceColor","r","color","r","Marker","square","MarkerSize",8);
hold on
errorbar(PConc,K_C,((Er_C.*(K_C+5))/2),"MarkerFaceColor","b","color","b","Marker","diamond");
errorbar(PConc,K_D,(Er_D.*(K_D+5))/2,"MarkerFaceColor","k","color","k","Marker","o");


title('Figure 4C')
xlabel('[Profilin] μM');
ylabel('Kpoly');
axis([0 18 -0.6 10.1]);
legend("pPB(18)","pPC(18)","pPD(18)");