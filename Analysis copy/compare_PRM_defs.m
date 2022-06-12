close all
figure(1);
PRM_bar = bar(PRMs(:,:));
 set(gca,'xtick',[1:34], 'xticklabel',["Diap1--Human";"Diap2--Human";"Diap3--Human";"Diap1--Mouse";"Diap2--Mouse";"Diap3--Mouse";"Diap1--Rat";"Diap3--Rat";"DAAM1--Human";"DAAM2--Human";"DAAM1--Mouse";"DAAM2--Mouse";"CAPU--FruitFly";"FMN1--Human";"FMN2--Human";"FMN1--Mouse";"FMN2--Mouse";"INF2--Mouse";"FHOD1--Human";"FHOD3--Human";"FHOD1--Mouse";"FHOD3--Mouse";"BNR1--Yeast";"CDC12P--Yeast";"BNI1P--Yeast";"FHODB--FruitFly";"Delphinin--Human";"Delphinin--Mouse";"FMNL1--Human";"FMNL2--Human";"FHDC1--Mouse";"DM7--FruitFly";"FMNL1--Mouse";"Delphinin--Zebrafish"]);
 xtickangle(90);
 set(PRM_bar(1), 'FaceColor','#0072BD');
 set(PRM_bar(2), 'FaceColor','#EDB120');
 set(PRM_bar(3), 'FaceColor','#7E2F8E');
 set(PRM_bar(4), 'FaceColor','#77AC30');
 set(PRM_bar(5), 'FaceColor','#A2142F');
 set(PRM_bar(6), 'FaceColor','#4DBEEE');
 legend( '3 w/ int','4 w/ int','6 w/ int','3 w/o int','4 w/o int','6 w/o int');
 xlabel('Formins');
 ylabel('# PRMs');
 
 title('Number of PRMs per Formin')

 
%% 
figure(2);
KP1_bar = bar(KP1(:,:));
 set(gca,'xtick',[1:34], 'xticklabel',["Diap1--Human";"Diap2--Human";"Diap3--Human";"Diap1--Mouse";"Diap2--Mouse";"Diap3--Mouse";"Diap1--Rat";"Diap3--Rat";"DAAM1--Human";"DAAM2--Human";"DAAM1--Mouse";"DAAM2--Mouse";"CAPU--FruitFly";"FMN1--Human";"FMN2--Human";"FMN1--Mouse";"FMN2--Mouse";"INF2--Mouse";"FHOD1--Human";"FHOD3--Human";"FHOD1--Mouse";"FHOD3--Mouse";"BNR1--Yeast";"CDC12P--Yeast";"BNI1P--Yeast";"FHODB--FruitFly";"Delphinin--Human";"Delphinin--Mouse";"FMNL1--Human";"FMNL2--Human";"FHDC1--Mouse";"DM7--FruitFly";"FMNL1--Mouse";"Delphinin--Zebrafish"]);
 xtickangle(90);
 set(PRM_bar(1), 'FaceColor','#0072BD');
 set(PRM_bar(2), 'FaceColor','#EDB120');
 set(PRM_bar(3), 'FaceColor','#7E2F8E');
 set(PRM_bar(4), 'FaceColor','#77AC30');
 set(PRM_bar(5), 'FaceColor','#A2142F');
 set(PRM_bar(6), 'FaceColor','#4DBEEE');
 legend( '3 w/ int','4 w/ int','6 w/ int','3 w/o int','4 w/o int','6 w/o int');
 xlabel('Formins');
 ylabel('log_2(kpoly)');
 ylim([0 12]);
 
 title('Polymerization rate per Formin')
 subtitle('single filament')
 
 %% 
figure(3);
KP2_bar = bar(KP2(:,:));
 set(gca,'xtick',[1:34], 'xticklabel',["Diap1--Human";"Diap2--Human";"Diap3--Human";"Diap1--Mouse";"Diap2--Mouse";"Diap3--Mouse";"Diap1--Rat";"Diap3--Rat";"DAAM1--Human";"DAAM2--Human";"DAAM1--Mouse";"DAAM2--Mouse";"CAPU--FruitFly";"FMN1--Human";"FMN2--Human";"FMN1--Mouse";"FMN2--Mouse";"INF2--Mouse";"FHOD1--Human";"FHOD3--Human";"FHOD1--Mouse";"FHOD3--Mouse";"BNR1--Yeast";"CDC12P--Yeast";"BNI1P--Yeast";"FHODB--FruitFly";"Delphinin--Human";"Delphinin--Mouse";"FMNL1--Human";"FMNL2--Human";"FHDC1--Mouse";"DM7--FruitFly";"FMNL1--Mouse";"Delphinin--Zebrafish"]);
 xtickangle(90);
 set(PRM_bar(1), 'FaceColor','#0072BD');
 set(PRM_bar(2), 'FaceColor','#EDB120');
 set(PRM_bar(3), 'FaceColor','#7E2F8E');
 set(PRM_bar(4), 'FaceColor','#77AC30');
 set(PRM_bar(5), 'FaceColor','#A2142F');
 set(PRM_bar(6), 'FaceColor','#4DBEEE');
 legend( '3 w/ int','4 w/ int','6 w/ int','3 w/o int','4 w/o int','6 w/o int');
 xlabel('Formins');
 ylabel('log_2(kpoly)');
 ylim([0 12]);
 
 title('Polymerization rate per Formin')
 subtitle('double filament')
 
%% 
figure(4);
KP3_bar = bar(KP3(:,:));
 set(gca,'xtick',[1:34], 'xticklabel',["Diap1--Human";"Diap2--Human";"Diap3--Human";"Diap1--Mouse";"Diap2--Mouse";"Diap3--Mouse";"Diap1--Rat";"Diap3--Rat";"DAAM1--Human";"DAAM2--Human";"DAAM1--Mouse";"DAAM2--Mouse";"CAPU--FruitFly";"FMN1--Human";"FMN2--Human";"FMN1--Mouse";"FMN2--Mouse";"INF2--Mouse";"FHOD1--Human";"FHOD3--Human";"FHOD1--Mouse";"FHOD3--Mouse";"BNR1--Yeast";"CDC12P--Yeast";"BNI1P--Yeast";"FHODB--FruitFly";"Delphinin--Human";"Delphinin--Mouse";"FMNL1--Human";"FMNL2--Human";"FHDC1--Mouse";"DM7--FruitFly";"FMNL1--Mouse";"Delphinin--Zebrafish"]);
 xtickangle(90);
 set(PRM_bar(1), 'FaceColor','#0072BD');
 set(PRM_bar(2), 'FaceColor','#EDB120');
 set(PRM_bar(3), 'FaceColor','#7E2F8E');
 set(PRM_bar(4), 'FaceColor','#77AC30');
 set(PRM_bar(5), 'FaceColor','#A2142F');
 set(PRM_bar(6), 'FaceColor','#4DBEEE');
 legend( '3 w/ int','4 w/ int','6 w/ int','3 w/o int','4 w/o int','6 w/o int');
 xlabel('Formins');
 ylabel('log_2(kpoly)');
 ylim([0 12]);
 
 title('Polymerization rate per Formin')
 subtitle('N terminal dimerized')

%% 
figure(5);
KPdimer_bar = bar(dimer(:,:));
 set(gca,'xtick',[1:34], 'xticklabel',["Diap1--Human";"Diap2--Human";"Diap3--Human";"Diap1--Mouse";"Diap2--Mouse";"Diap3--Mouse";"Diap1--Rat";"Diap3--Rat";"DAAM1--Human";"DAAM2--Human";"DAAM1--Mouse";"DAAM2--Mouse";"CAPU--FruitFly";"FMN1--Human";"FMN2--Human";"FMN1--Mouse";"FMN2--Mouse";"INF2--Mouse";"FHOD1--Human";"FHOD3--Human";"FHOD1--Mouse";"FHOD3--Mouse";"BNR1--Yeast";"CDC12P--Yeast";"BNI1P--Yeast";"FHODB--FruitFly";"Delphinin--Human";"Delphinin--Mouse";"FMNL1--Human";"FMNL2--Human";"FHDC1--Mouse";"DM7--FruitFly";"FMNL1--Mouse";"Delphinin--Zebrafish"]);
 xtickangle(90);
 set(PRM_bar(1), 'FaceColor','#0072BD');
 set(PRM_bar(2), 'FaceColor','#EDB120');
 set(PRM_bar(3), 'FaceColor','#7E2F8E');
 set(PRM_bar(4), 'FaceColor','#77AC30');
 set(PRM_bar(5), 'FaceColor','#A2142F');
 set(PRM_bar(6), 'FaceColor','#4DBEEE');
 legend( '3 w/ int','4 w/ int','6 w/ int','3 w/o int','4 w/o int','6 w/o int');
 xlabel('Formins');
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)');
 ylim([-1.9 0.1]);

 
 title('Change in Polymerization rate per Formin')
 