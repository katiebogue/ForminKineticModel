function out_tab= PRM_size_sweep_info(kpoly_sweep_table,resultspath)
colors=distinguishable_colors(20,'w');
out_tab=[];
for id=1224:-1:1050
    disp(id)
    strid=num2str(id);
    tab1=kpoly_sweep_table(kpoly_sweep_table.id==id,:);
    %id_info=tab1{1,8:13};
    for FH1_length=10:20:600
        tab=tab1(tab1.FH1_length==FH1_length,:);
        %[G,ID]=findgroups(tab.PRM_loc);
        h=gscatter(tab.PRM_size,log2(tab.kpoly_ratio),tab.PRM_loc,colors,'.',15,'on','PRM size','Kpoly ratio');
        lgd=legend('Location','northeastoutside');
        lgd.Title.String = 'PRM to FH2 dist';
        set(h, 'linestyle', '-');
        title(strcat("FH1 length:",num2str(FH1_length)," id:",strid))
        title1=strcat("FH1length_",num2str(FH1_length),"_id_",strid,".fig");
        %figuresave(gcf,opts,title1)
        if id>1150
            exportgraphics(gcf,fullfile(resultspath,"RESULTS_6.pdf"),'Append',true);
        elseif id>1100
            exportgraphics(gcf,fullfile(resultspath,"RESULTS_7.pdf"),'Append',true);
        else
            exportgraphics(gcf,fullfile(resultspath,"RESULTS_8.pdf"),'Append',true);
        end
        savefig(gcf,fullfile(resultspath,"figures",title1));
        saveas(gcf,fullfile(resultspath,"figures",append(extractBefore(title1,strlength(title1)-3),'.png')));
        % for i=1:length(ID)
        %     ind=tab.PRM_loc==ID(i);
        %     ydata=log2(tab.kpoly_ratio(ind));
        %     xdata=tab.PRM_size(ind);
        %     ydata = ydata( ~any( isnan( ydata ) | isinf( ydata ), 2 ),: );
        %     xdata = xdata( ~any( isnan( ydata ) | isinf( ydata ), 2 ),: );
        %     if length(ydata)>3
        %         max_val=max(ydata);
        %         min_val=min(ydata);
        %         max_x=xdata(ydata==max_val);
        %         min_x=xdata(ydata==min_val);
        %         if length(max_x)>1 
        %             max_x=max(max_x);
        %         end
        %         if length(min_x)>1 
        %             min_x=min(min_x);
        %         end
        %         if max_x>min_x
        %             max_x_boundry=min(xdata((abs(ydata - max_val))<= 0.05*(std(ydata))));
        %             min_x_boundry=max(xdata((abs(ydata - min_val))<= 0.05*(std(ydata))));
        %         else
        %             max_x_boundry=max(xdata((abs(ydata - max_val))<= 0.05*(std(ydata))));
        %             min_x_boundry=min(xdata((abs(ydata - min_val))<= 0.05*(std(ydata))));
        %         end
        %         entry=[id,FH1_length,ID(i),id_info,max_val,min_val,max_x,max_x_boundry,min_x,min_x_boundry];
        %         out_tab=[out_tab; entry];
        %     end
        % end
        clear h
    end
    
end
% out_tab=array2table(out_tab);
% out_tab.Properties.VariableNames=["id","FH1_length","PRM_loc","k_cap","k_del","r_cap","r_del","k_rel","model","max_val","min_val","max_x","max_x_boundry","min_x","min_x_boundry"];
% save('PRM_sweep_table.mat',"out_tab",'-mat')
end
