%REPLICATEFIG3 reads in values from csv file and uses them to replicate
%courtemanche figure 3 scatter plot with error bars
% 
%   Reads in data from "Courtemanche plot digitizer data F3.csv"
%
%   See also REPLICATEFIG4A, REPLICATEFIG4C.

values=readmatrix("Courtemanche plot digitizer data F3.csv");
PConc=values(:,1);
K_18=values(:,2);
Er_18=values(:,3);
K_65=values(:,6);
Er_65=values(:,7);
K_37=values(:,10);
Er_37=values(:,11);
K_28=values(:,14);
Er_28=values(:,15);

hold off
errorbar(PConc,K_65,(Er_65.*(K_65+5))/2,"MarkerFaceColor","r","color","r","Marker","o");
hold on
errorbar([0,0.5,1,5,10,17.5],rmmissing(K_37),rmmissing((Er_37.*(K_37+5))/2),"MarkerFaceColor","b","color","b","Marker","o");
errorbar(PConc,K_28,(Er_28.*(K_28+5))/2,"MarkerFaceColor","#065f18","color","#065f18","Marker","o");
errorbar(PConc,K_18,(Er_18.*(K_18+5))/2,"MarkerFaceColor","k","color","k","Marker","o");


title('Figure 3')
xlabel('[Profilin] μM');
ylabel('Kpoly');
axis([0 20 -0.5 10.1]);
legend("pPD(65)","pPD(37)","pPD(28)","pPD(18)");