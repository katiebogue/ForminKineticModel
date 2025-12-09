saveTF=true;
savefigfolder="/Users/katiebogue/MATLAB/GitHub/Data/polymer-c_data/statheatmaps/calc_delloc/prmslab"; polymer_heatmaps_paper
%set(groot, 'DefaultTextFontSize', 50);
%%
ltvar="Prvec0";
calcTF=true;
lookuptab=0;
FH2_dist="35.5";
types=["dimer","double","ratio"];

xlocs=[
6.1839697
6.1839697
8.825
8.825
0
17.75
17.75
8.333
8.333
8.875
8.875];

ylocs=[
12.5771768
0
11.154
0
0
0
16.667
0
16.667
0
16.667];

for i=1:length(ylocs)
    yloc=ylocs(i);
    xloc=xlocs(i);
    for j=1:length(types)
        type1=types(j);
        if type1=="dimer"
            limits=[-10 1];
        elseif type1=="double"
            limits=[-10 1];
        elseif type1=="ratio"
            limits=[-10 10];
        end
        polymerstatheatmap(ltvar,lookuptab,FH2_dist,type1,calcTF,saveTF,savefigfolder,limits,400,400,0,0,xloc,yloc)
    end
end
%%
types=["dimer","double","ratio"];
ltvars=["Prvec0","POcclude"];
lookuptabs=[lt_16_67,lt_35_5,lt_0_16,lt_0,lt_0_35,lt_35_16];
FH2_dist_vals=["16.67","35.5","0, 16.67","0","0, 35.5","35.5, 16.67"];
calcTFs=[true,true,false,true,false,false];

for i=1:length(types)
    type1=types(i);
    if type1=="dimer"
        limits=[-10 1];
    elseif type1=="double"
        limits=[-4.5 0.5];
    elseif type1=="ratio"
        limits=[-6 6];
    end
    for j=1:length(ltvars)
        ltvar=ltvars(j);
        for k=1:length(lookuptabs)
            lookuptab=lookuptabs(k);
            FH2_dist=FH2_dist_vals(k);
            calcTF=false;
            polymerstatheatmap(ltvar,lookuptab,FH2_dist,type1,calcTF,saveTF,savefigfolder,limits,Experiment1)

            if calcTFs(k)
                if ltvar=="Prvec0"
                    calcTF=true;
                    polymerstatheatmap(ltvar,lookuptab,FH2_dist,type1,calcTF,saveTF,savefigfolder,limits,Experiment1)
                end
            end
        end

    end
end