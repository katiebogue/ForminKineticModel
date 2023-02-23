values=readmatrix("Courtemanche plot digitizer data F4A.csv");
PConc=values(:,1);
K_B=values(:,2);
Er_B=values(:,3);
K_C=values(:,6);
Er_C=values(:,7);
K_D=values(:,10);
Er_D=values(:,11);

idx = ~isnan(K_C);

close all
hold off
errorbar(PConc,K_B,(Er_B.*(K_B+5))/2,"MarkerFaceColor","r","color","r","Marker","square","MarkerSize",8);
hold on
errorbar(PConc(idx),K_C(idx),((Er_C(idx).*(K_C(idx)+5))/2),"MarkerFaceColor","b","color","b","Marker","diamond");
errorbar(PConc(idx),K_D(idx),(Er_D(idx).*(K_D(idx)+5))/2,"MarkerFaceColor","k","color","k","Marker","o");


title('Figure 4A')
xlabel('[Profilin] μM');
ylabel('Kpoly');
axis([0 18 -0.3 10.6]);
legend("pPB(37)","pPC(37)","pPD(37)");