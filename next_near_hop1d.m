function H_ind=next_near_hop1d(t, k, rL, i_u, i_d, range, ii, Li, l_d, aux_s, flag_s, arg_s, M_s, M_aux)
%% Input
%   t: hopping term
%   k: twitst angle
%   r: site of lattice
%   i_u: configuration for the up sector
%   i_d: configuration for the down sector
%   range: conversion array for the desired sector
%   ii: site on the desired dimension
%   Li: size of the desired dimension
%   l_d: dimension of down sector
%   aux_s: auxiliary scalar to determine adjacent sites in the given
%   dimension
%   flag_s: logic value for determine indexation
%   arg_s: i_u or i_d depending on the sector
%   M_s: Size of lattice
    H_ind=zeros(2, 3);
    ind1=(i_u-1)*l_d+i_d;
    if ii-2<1
        if bitand(range(arg_s),2^(rL+(Li-2)*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL+(Li-2)*aux_s),1);
            N_fact=sum(M_aux(M_s-rL-(Li-2)*aux_s:M_s-rL-1));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(1,1)=ind1;
            H_ind(1, 2)=ind2;
            H_ind(1, 3)=-t*exp(k)*sgn;
        end
        
        if bitand(range(arg_s),2^(rL+2*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL+2*aux_s),1);
            N_fact=sum(M_aux(M_s-rL-2*aux_s:M_s-rL-1));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(2,1)=ind1;
            H_ind(2, 2)=ind2;
            H_ind(2, 3)=-t*sgn;
        end
    elseif ii+2>Li
        if bitand(range(arg_s),2^(rL-2*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL-2*aux_s),1);
            N_fact=sum(M_aux(M_s-rL+1:M_s-(rL-2*aux_s)));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(1,1)=ind1;
            H_ind(1, 2)=ind2;
            H_ind(1, 3)=-t*sgn;
        end
        if bitand(range(arg_s),2^(rL-aux_s*(Li-2)))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL-aux_s*(Li-2)),1);
            N_fact=sum(M_aux(M_s-rL+1:M_s-(rL-(Li-2)*aux_s)));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(2,1)=ind1;
            H_ind(2, 2)=ind2;
            H_ind(2, 3)=-t*exp(-k)*sgn;
        end
    else
        if bitand(range(arg_s),2^(rL+2*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL+2*aux_s),1);
            N_fact=sum(M_aux(M_s-rL-2*aux_s:M_s-rL-1));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(1,1)=ind1;
            H_ind(1, 2)=ind2;
            H_ind(1, 3)=-t*sgn;
        end
        
        if bitand(range(arg_s),2^(rL-2*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL-2*aux_s),1);
            N_fact=sum(M_aux(M_s-rL+1:M_s-rL+2*aux_s));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(2,1)=ind1;
            H_ind(2, 2)=ind2;
            H_ind(2, 3)=-t*sgn;
        end
    end
end