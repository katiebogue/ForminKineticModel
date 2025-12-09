% CHECKCOUNTS checks for zero counts in prvec values of the object in the
% workspace lt

prvec_types=["Prvec0","Prvec0_halfup","Prvec0_halfup_op","Prvec0_op","Prvec0_up","Prvec0_up_op","Prvec_cen","Prvec_cen_halfup","Prvec_cen_up","Prvec_offcen","Prvec_offcen_halfup","Prvec_offcen_halfup_op","Prvec_offcen_op","Prvec_offcen_up","Prvec_offcen_up_op"];
loctypes=["","_0_100000","_0_250000","_0_500000","_0_750000","_1_000000","_2_000000","_4_000000","_10_000000","_20_000000","_36_000000"];


for i=1:length(prvec_types)
    for j=1:length(loctypes)
        valtab=lt.stattable(strcat(prvec_types(i),"_sum",loctypes(j)),"double");
        valtab2=lt.stattable(strcat(prvec_types(i),loctypes(j)),"double");

        zero_locs.a=(valtab2.a(:,3)==0);
        zero_locs.b=(valtab2.b(:,3)==0);

        if (valtab2.a(zero_locs.a,3)~=0)
            errlocs.a=(valtab2.a(zero_locs.a,3)~=0);
            errs.a=valtab2.a(errlocs.a,1:2);
        else
            errs.a=0;
        end
        if (valtab2.b(zero_locs.b,3)~=0)
            errlocs.b=(valtab2.b(zero_locs.b,3)~=0);
            errs.b=valtab2.b(errlocs.b,1:2);
        else
            errs.b=0;
        end

        if errs.a || errs.b
            errstruct.(strcat(prvec_types(i),loctypes(j)))=errs;
        end
    end
end