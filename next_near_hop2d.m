function H_ind=next_near_hop2d(t, k1, k2, r, i_up, i_dn, range, ii, jj, Li, Lj, l_dn, aux1, aux2, flag_s, arg_s, M_s, M_aux)
%%  Calculates next nearest neighbors hooping term for a higher dimensional lattice
%   Input:
%   t: hooping term
%   k1: twist angle for one of the axis
%   k2: twits angle for the remaining angle
%   r: site of lattice
%   i_up: configuration for the spin up sector
%   i_dn: configuration for the spin down sector
%   range: vector for indexing conversion
%   ii: position on first axis
%   jj: position on second axis
%   Li: lenght of first axis
%   Lj: lenght of second axis
%   l_dn: lenght of sector down
%   aux1: auxiliar scalar for indexation over the first axis
%   aux2: auxiliar scalar for indexation over the second axis
%   flag_s: boolean variable for indexation
%   arg_s: variable for indentifying the sector
%   M_s: size of lattice
%   Output: 
%   H_ind: hopping terms for state i_up or i_dn related to site (ii, jj)
    H_ind=zeros(4, 3);
    ind1=(i_up-1)*l_dn+i_dn;
    if ii==1
        if jj==1
            if bitand(range(arg_s),2^(r+(Li-1)*aux1+(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+(Li-1)*aux1+(Lj-1)*aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+(Li-1)*aux1+(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(k1)*exp(k2)*sgn;
            end
            
            if bitand(range(arg_s),2^(r+(Li-1)*aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+(Li-1)*aux1+aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+(Li-1)*aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(k1)*sgn;
            end
            
            if bitand(range(arg_s),2^(r+aux1+(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1+(Lj-1)*aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1+(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*exp(k2)*sgn;
            end
            
            if bitand(range(arg_s),2^(r+aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1+aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;  
            end
        
        elseif jj==Lj
            if bitand(range(arg_s),2^(r+(Li-1)*aux1-(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+(Li-1)*aux1-(Lj-1)*aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+(Li-1)*aux1-(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(k1)*exp(-k2)*sgn;  
            end
            if bitand(range(arg_s),2^(r+aux1-(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1-(Lj-1)*aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1-(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(-k2)*sgn;
            end
                
            if  bitand(range(arg_s),2^(r+(Li-1)*aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+(Li-1)*aux1-aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+(Li-1)*aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*exp(k1)*sgn;  
            end
            
            if  bitand(range(arg_s),2^(r+aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1-aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;  
            end
        else
            if bitand(range(arg_s),2^(r+(Li-1)*aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+(Li-1)*aux1+aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+(Li-1)*aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(k1)*sgn; 
                
            end
            if bitand(range(arg_s),2^(r+(Li-1)*aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+(Li-1)*aux1-aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+(Li-1)*aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(k1)*sgn; 
                
            end
            if bitand(range(arg_s),2^(r+aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1+aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*sgn; 
            end
            if bitand(range(arg_s),2^(r+aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1-aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;     
            end
        
        end   
    elseif ii==Li
        if jj==1
            if bitand(range(arg_s),2^(r-(Li-1)*aux1+(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-(Li-1)*aux1+(Lj-1)*aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-(Li-1)*aux1+(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(-k1)*exp(k2)*sgn;  
            end
                
            if bitand(range(arg_s),2^(r-(Li-1)*aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-(Li-1)*aux1+aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-(Li-1)*aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(-k1)*sgn;
            end
            if bitand(range(arg_s),2^(r-aux1+(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1+(Lj-1)*aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1+(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*exp(k2)*sgn;
            end
            if bitand(range(arg_s),2^(r-aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1+aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r,r-aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;    
            end
        elseif jj==Lj
            if bitand(range(arg_s),2^(r-(Li-1)*aux1-(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-(Li-1)*aux1-(Lj-1)*aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r,r-(Li-1)*aux1-(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(-k1)*exp(-k2)*sgn;
            end
            if bitand(range(arg_s),2^(r-(Li-1)*aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-(Li-1)*aux1-aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r,r-(Li-1)*aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(-k1)*sgn;
            end
            if bitand(range(arg_s),2^(r-aux1-(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1-(Lj-1)*aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1-(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*exp(-k2)*sgn;
            end
                    
            if bitand(range(arg_s),2^(r-aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1-aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;
            end
        else
            if bitand(range(arg_s),2^(r-(Li-1)*aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-(Li-1)*aux1+aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-(Li-1)*aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(-k1)*sgn;
            end
            if bitand(range(arg_s),2^(r-(Li-1)*aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-(Li-1)*aux1-aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-(Li-1)*aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(-k1)*sgn;
            end
            if bitand(range(arg_s),2^(r-aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1-aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*sgn;
            end
            if bitand(range(arg_s),2^(r-aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1+aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;
            end                    
        end
    elseif jj==1
        if mod(ii,Li)>1
            if bitand(range(arg_s),2^(r-aux1+(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1+(Lj-1)*aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1+(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(k2)*sgn;
            end
            if bitand(range(arg_s),2^(r+aux1+(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1+(Lj-1)*aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1+(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(k2)*sgn;
            end
            if bitand(range(arg_s),2^(r-aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1+aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*sgn;
            end
            if bitand(range(arg_s),2^(r+aux1+aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1+aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1+aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;       
            end     
        end
    elseif jj==Lj
        if mod(ii,Li)>1
            if bitand(range(arg_s),2^(r-aux1-(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1-(Lj-1)*aux2),1);
                [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1-(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(1, 3)=-t*exp(-k2)*sgn;    
            end
            if bitand(range(arg_s),2^(r+aux1-(Lj-1)*aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1-(Lj-1)*aux2),1);
                [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1-(Lj-1)*aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(2, 3)=-t*exp(-k2)*sgn;    
            end
            if bitand(range(arg_s),2^(r-aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r-aux1-aux2),1);
                [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(3, 3)=-t*sgn;   
            end
            if bitand(range(arg_s),2^(r+aux1-aux2))==0
                ind2=find(range==range(arg_s)-2^r+2^(r+aux1-aux2),1);
                [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1-aux2, M_aux, l_dn, i_up, i_dn);
                H_ind(4, 3)=-t*sgn;       
            end
        end                        
    else 
        if bitand(range(arg_s),2^(r+aux1+aux2))==0
            ind2=find(range==range(arg_s)-2^r+2^(r+aux1+aux2),1);
            [H_ind(1,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1+aux2, M_aux, l_dn, i_up, i_dn);
            H_ind(1, 3)=-t*sgn;  
        end
        if bitand(range(arg_s),2^(r-aux1+aux2))==0
            ind2=find(range==range(arg_s)-2^r+2^(r-aux1+aux2),1);
            [H_ind(2,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1+aux2, M_aux, l_dn, i_up, i_dn);
            H_ind(2, 3)=-t*sgn;  
        end
        if bitand(range(arg_s),2^(r+aux1-aux2))==0
            ind2=find(range==range(arg_s)-2^r+2^(r+aux1-aux2),1);
            [H_ind(3,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r+aux1-aux2, M_aux, l_dn, i_up, i_dn);
            H_ind(3, 3)=-t*sgn;
        end
        if bitand(range(arg_s),2^(r-aux1-aux2))==0
            ind2=find(range==range(arg_s)-2^r+2^(r-aux1-aux2),1);
            [H_ind(4,1:2), sgn]=Index_hop(ind1, ind2, flag_s, M_s , r, r-aux1-aux2, M_aux, l_dn, i_up, i_dn);
            H_ind(4, 3)=-t*sgn;
        end
    end
end
