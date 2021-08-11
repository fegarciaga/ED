function [Ind, sgn]=Index_hop(ind1, ind2, flag, M , r, r_alt, M_aux, l_dn, i_up, i_dn)
    Ind=zeros(1,2);
    if min(M-r, M-r_alt)==M-r
        if r_alt+1<=r-1
            N_fact=sum(M_aux(M-r+1:M-r_alt));
        else
            N_fact=0;
        end
    else
        if r+1<=r_alt-1
            N_fact=sum(M_aux(M-r_alt:M-r-1));
        else
            N_fact=0;
        end
    end
    if flag==1
        ind=(ind2-1)*l_dn+i_dn;
    else
        ind=(i_up-1)*l_dn+ind2;
    end
    sgn=(-1)^(N_fact);
    Ind(1,1)=ind1;
    Ind(1, 2)=ind;
end