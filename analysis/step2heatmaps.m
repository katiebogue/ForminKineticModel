% STEP2HEATMAPS generates heatmaps of polymerization rate ratios (dimer/nondimer)
%for sweeps across possible probability density and occlusion probability
%ratios for a set parameter regime.

set(groot,'defaultfigureposition',[400 250 700 550]) % helps prevent cut offs in figs


kcap=@(p_occ,c_PA,k_cap) (p_occ).*c_PA.*k_cap;
kdel=@(p_r,G,k_del,pocc_base) (1-pocc_base).*p_r.*k_del.*G;
rcap=@(r_cap,r_cap_exp,size) r_cap.*exp(-size.*r_cap_exp);

kpolyfun= @(p_occ,c_PA,k_cap,p_r,G,k_del,pocc_base,r_cap,r_cap_exp,size) 1./((1./kdel(p_r,G,k_del,pocc_base)) + ((kdel(p_r,G,k_del,pocc_base) + rcap(r_cap,r_cap_exp,size))./(kdel(p_r,G,k_del,pocc_base).*kcap(p_occ,c_PA,k_cap))));
kpoly_ratio= @(p_occ_double, p_occ_dimer, p_r_double, p_r_dimer, pocc_base_double, pocc_base_dimer,c_PA,k_cap,G,k_del,r_cap,r_cap_exp,size) kpolyfun(p_occ_dimer,c_PA,k_cap,p_r_dimer,G,k_del,pocc_base_dimer,r_cap,r_cap_exp,size) ./ kpolyfun(p_occ_double,c_PA,k_cap,p_r_double,G,k_del,pocc_base_double,r_cap,r_cap_exp,size);

c_PA=1;
k_cap=10000000;
G=1;
k_del=1;
r_cap=10000;
r_cap_exp=0.8;
size=10;
pocc_base_double=0.75;
pocc_base_dimer=0.75;

p_occ_double=0.5;
p_r_double=10^5;

kpoly_ratio_sweep= @(p_occ_ratio, p_r_ratio) log2(kpoly_ratio(p_occ_double, p_occ_double.*(2.^p_occ_ratio), p_r_double, p_r_double.*(2.^p_r_ratio), pocc_base_double, pocc_base_dimer,c_PA,k_cap,G,k_del,r_cap,r_cap_exp,size));

kpoly_lab="log_{2}(k_{poly} N terminal dimerized/k_{poly} double)";

p_r_ratio_sweep_vals=-5:0.05:5;

p_occ_ratio_sweep_vals=[-1:0.005:-0.0000000001 0:0.005:log2(1/p_occ_double)];



% Generate all combinations of i and j using ndgrid
[I, J] = ndgrid(p_occ_ratio_sweep_vals, p_r_ratio_sweep_vals);

% Flatten the grids
I_flat = I(:);
J_flat = J(:);

% Preallocate kpval array
Kpvals = arrayfun(@(x, y) kpoly_ratio_sweep(x, y), I_flat, J_flat);

% Combine into result matrix
vals = [I_flat, J_flat, Kpvals];

%vals(imag(vals)~=0) = nan;
x=array2table(vals);
figure;
h=heatmap(x,"vals1","vals2",'ColorVariable',"vals3","ColorMethod","none","GridVisible","off");

load('customcolorbar_red_blue.mat');
h.Colormap=CustomColormap;
% minn=min((x.vals3(all([x.vals3~=0 x.vals3~=-Inf],2))));
% h.ColorLimits=[minn abs(minn)];
h.ColorLimits=[-5 5];
h.NodeChildren(3).YDir='normal';


% Convert each number in the array into a string
CustomXLabels = string(p_occ_ratio_sweep_vals);
CustomYLabels = string(p_r_ratio_sweep_vals);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(p_occ_ratio_sweep_vals,0.5) ~= 0) = " ";
CustomYLabels(mod(p_r_ratio_sweep_vals,1) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;
s = struct(h); 
s.XAxis.TickLabelRotation = 0;   % horizontal

lab_rates=strcat("kcap: ",num2str(k_cap)," kdel: ",num2str(k_del)," rcap: ",num2str(r_cap)," rcap_{exp}: ",num2str(r_cap_exp));
lab_formin=strcat("c_{PA}: ",num2str(c_PA)," G: ",num2str(G)," PRM size: ",num2str(size));
lab_sim=strcat("p_{occ} double: ",num2str(p_occ_double)," p_{r} double: ",num2str(p_r_double)," p_{occ}^{0} double: ",num2str(pocc_base_double)," p_{occ}^{0} dimer: ",num2str(pocc_base_dimer));

h.Title = {kpoly_lab,lab_rates,lab_formin,lab_sim};
h.YLabel = "log_{2}(p_{r} N terminal dimerized/p_{r} double)";
h.XLabel = "log_{2}(1-p_{occ} N terminal dimerized/1-p_{occ} double)";

lab_rates=strcat("kcap_",num2str(k_cap)," kdel_",num2str(k_del)," rcap_",num2str(r_cap)," rcapexp_",num2str(r_cap_exp));
lab_formin=strcat("cPA_",num2str(c_PA)," G_",num2str(G)," PRM size_",num2str(size));
lab_sim=strcat("p_occ_dob_",num2str(p_occ_double)," p_r_dob_",num2str(p_r_double)," p_occ_0_dob_",num2str(pocc_base_double)," p_occ_0_dim_",num2str(pocc_base_dimer));
exportgraphics(gca,strcat("step2heatmap_",lab_rates,lab_formin,lab_sim,".png"),"Resolution",750)
exportgraphics(gca,strcat("step2heatmap_",lab_rates,lab_formin,lab_sim,".eps"),'BackgroundColor','none','ContentType','vector');


