function H_ind=near_hop(t, k, rL, i_u, i_d, range, ii, Li, l_d, aux_s, flag_s, arg_s, M_s, M_aux) 
%% Description
%   Calculates the index due to nearest neighbor hopping.
%   Input
%   t: hopping term
%   k: twist angle
%   r: position of the electron in lattice
%   i_u: index of state up
%   i_d: index of state dn
%   range: vector with the needed indexation
%   ii: possition of the electron in axis
%   Li: size of axis
%   l_d: size of the dn sector
%   aux_s: scalar to specify axis
%   flag_s: indicates the correspondent sector (usefull for indexation)
%   arg_s: i_u or i_d
%   M_s: size of lattice

%% implementation:
    H_ind=zeros(2, 3);
    % creates index corresponding to the given state
    % check for boundaries
    if ii==1
        % check if the hopping site is occupied, if not calculates the
        % index corresponding to the new state
        if bitand(range(arg_s),2^(rL+(Li-1)*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL+(Li-1)*aux_s));
            N_fact=sum(M_aux(M_s-(rL+(Li-1)*aux_s):M_s-rL-1));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(1,1)=(i_u-1)*l_d+i_d;
            H_ind(1, 2)=ind2;
            H_ind(1, 3)=-t*exp(k)*sgn;
       
        end
        if bitand(range(arg_s),2^(rL+aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL+aux_s));
            N_fact=0;
            if rL+1<=rL+aux_s-1
                N_fact=sum(M_aux(M_s-(rL+aux_s):M_s-rL-1));
            end
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(2,1)=(i_u-1)*l_d+i_d;
            H_ind(2, 2)=ind2;
            H_ind(2, 3)=-t*sgn;
        end
        
    elseif ii==Li
        if bitand(range(arg_s),2^(rL-aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL-aux_s));
            N_fact=0;
            if rL-aux_s+1<=rL-1
                N_fact=sum(M_aux(M_s-rL+1:M_s-(rL-aux_s)));
            end
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(1,1)=(i_u-1)*l_d+i_d;
            H_ind(1, 2)=ind2;
            H_ind(1, 3)=-t*sgn;
        end
        if bitand(range(arg_s),2^(rL-(Li-1)*aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL-(Li-1)*aux_s));
            N_fact=sum(M_aux(M_s-rL+1:M_s-(rL-(Li-1)*aux_s)));
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(2,1)=(i_u-1)*l_d+i_d;
            H_ind(2, 2)=ind2;
            H_ind(2, 3)=-t*exp(-k)*sgn;
        end
    else
        if bitand(range(arg_s),2^(rL-aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL-aux_s));
            N_fact=0;
            if rL-aux_s+1<=rL-1
                N_fact=sum(M_aux(M_s-rL+1:M_s-(rL-aux_s)));
            end
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(1,1)=(i_u-1)*l_d+i_d;
            H_ind(1, 2)=ind2;
            H_ind(1, 3)=-t*sgn;
        end
        
        if bitand(range(arg_s),2^(rL+aux_s))==0
            ind2=find(range==range(arg_s)-2^rL+2^(rL+aux_s));
            N_fact=0;
            if rL+1<=rL+aux_s-1
                N_fact=sum(M_aux(M_s-(rL+aux_s):M_s-rL-1));
            end
            if flag_s==1
                ind2=(ind2-1)*l_d+i_d;
            else
                ind2=(i_u-1)*l_d+ind2;
            end
            sgn=(-1)^(N_fact);
            H_ind(2,1)=(i_u-1)*l_d+i_d;
            H_ind(2, 2)=ind2;
            H_ind(2, 3)=-t*sgn;
        end
    end
end