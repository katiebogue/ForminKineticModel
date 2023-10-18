classdef Formin <handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name
        pp_length_vec
        p_length_vec
        pp_index_vec
        fh1_length
        o
        c_PA
    end

    methods
        function obj = Formin(name,pp_length_vec,p_length_vec, pp_index_vec,fh1_length,lt,c_PA)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = name;
            obj.pp_length_vec=pp_length_vec;
            obj.p_length_vec=p_length_vec;
            obj.pp_index_vec=pp_index_vec;
            obj.fh1_length=fh1_length;
            obj.o=polymerstats(obj,lt);
            obj.c_PA=c_PA;
        end

        function rate = kdob(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            % scale k_paf_rev with e^-L
            rate=obj.kdb(k_paf,k_pab,k_paf_rev.*exp(-1.*obj.pp_length_vec),r_PF_rev,r_paf_rev);
        end

        function rate = kdim(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            % scale k_paf_rev with e^-L
            rate=obj.kdm(k_paf,k_pab,k_paf_rev.*exp(-1.*obj.pp_length_vec),r_PF_rev,r_paf_rev);
        end

        function rate = kdob_rPFrev(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            % scale k_paf_rev with e^-L; scale r_PF_rev with e^-L
            rate=obj.kdb(k_paf,k_pab,k_paf_rev.*exp(-1.*obj.pp_length_vec),r_PF_rev.*exp(-1.*obj.pp_length_vec),r_paf_rev);
        end

        function rate = kdim_rPFrev(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            % scale k_paf_rev with e^-L; scale r_PF_rev with e^-L
            rate=obj.kdm(k_paf,k_pab,k_paf_rev.*exp(-1.*obj.pp_length_vec),r_PF_rev.*exp(-1.*obj.pp_length_vec),r_paf_rev);
        end

        function rate = kdob_rPFrev_kpaf(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            % scale k_paf with L; scale r_PF_rev with e^-L
            rate=obj.kdb(k_paf.*obj.pp_length_vec,k_pab,k_paf_rev,r_PF_rev.*exp(-1.*obj.pp_length_vec),r_paf_rev);
        end

        function rate = kdim_rPFrev_kpaf(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            % scale k_paf with L; scale r_PF_rev with e^-L
            rate=obj.kdm(k_paf.*obj.pp_length_vec,k_pab,k_paf_rev,r_PF_rev.*exp(-1.*obj.pp_length_vec),r_paf_rev);
        end

        function rate = kdb(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            ratea = sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-obj.o.p_occ2a_0(:)).*(1.0e33.*obj.o.p_r2a(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-obj.o.p_occ2a_0(:)).*(1.0e33.*obj.o.p_r2a(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*obj.c_PA.*(1-obj.o.p_occ2a(:))).* (k_pab.*(1-obj.o.p_occ2a_0(:)).*(1.0e33.*obj.o.p_r2a(:)/(27*6.022e23))) .* r_PF_rev))));
            rateb= sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-obj.o.p_occ2b_0(:)).*(1.0e33.*obj.o.p_r2b(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-obj.o.p_occ2b_0(:)).*(1.0e33.*obj.o.p_r2b(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*obj.c_PA.*(1-obj.o.p_occ2b(:))).* (k_pab.*(1-obj.o.p_occ2b_0(:)).*(1.0e33.*obj.o.p_r2b(:)/(27*6.022e23))) .* r_PF_rev))));
            rate=ratea+rateb;
        end

        function rate = kdm(obj,k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            ratea = sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-obj.o.p_occ3a_0(:)).*(1.0e33.*obj.o.p_r3a(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-obj.o.p_occ3a_0(:)).*(1.0e33.*obj.o.p_r3a(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*obj.c_PA.*(1-obj.o.p_occ3a(:))).* (k_pab.*(1-obj.o.p_occ3a_0(:)).*(1.0e33.*obj.o.p_r3a(:)/(27*6.022e23))) .* r_PF_rev))));
            rateb= sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-obj.o.p_occ3b_0(:)).*(1.0e33.*obj.o.p_r3b(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-obj.o.p_occ3b_0(:)).*(1.0e33.*obj.o.p_r3b(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*obj.c_PA.*(1-obj.o.p_occ3b(:))).* (k_pab.*(1-obj.o.p_occ3b_0(:)).*(1.0e33.*obj.o.p_r3b(:)/(27*6.022e23))) .* r_PF_rev))));
            rate=ratea+rateb;
        end

        
    end
end