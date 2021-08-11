function [range_up, range_dn, l_up, l_dn] = Index(N_u, N_d, lx, ly, lz)
%This function implements the indexation of hilbert space to meaningfull
% integers
%Input
%N_u: Number of spin up electrons
%N_d: Number of spin down electrons
%lx, ly, lz: size of the system
%Output
%range_up: array that relates index number to meaningfull integer for spin
% up
%range_dn: array that relates index number to meaningfull integer for spin
% down
%l_up: Dimension of Hilbert space for spin up sector
%l_dn: Dimension of Hilbert space for spin down sector
    N_sites = lx*ly*lz;
    l_up = nchoosek(N_sites, N_u);
    aux=zeros(N_u,1);
    range_up=zeros(l_up,1);
    aux_num=0;
    for i=1:N_u
        aux(i,1)=N_sites-i;
        aux_num=aux_num+2^(aux(i,1));
    end
    range_up(1,1)=aux_num;

    for i=1:l_up-1
        aux_num=0;
        i_aux=N_u;
        flag=1;
        while flag==1
           if i_aux==N_u
                if aux(i_aux)>0
                   flag=0; 
                else
                    i_aux=i_aux-1;
                end
           elseif aux(i_aux,1)-aux(i_aux+1,1)>1
               flag=0;
           else
                i_aux=i_aux-1;
           end    
        end
        aux(i_aux,1)=aux(i_aux,1)-1;
        if(i_aux~=N_u)
            for j=i_aux+1:N_u
                aux(j,1)=aux(j-1,1)-1;
            end
        end
    
        for j=1:N_u
            aux_num=aux_num+2^aux(j,1);
        end
        range_up(i+1,1)=aux_num;
    end

    l_dn = nchoosek(N_sites, N_d);
    aux=zeros(N_d,1);
    range_dn=zeros(l_dn,1);
    aux_num=0;

    for i=1:N_d
        aux(i,1)=N_sites-i;
        aux_num=aux_num+2^(aux(i,1));
    end
    range_dn(1,1)=aux_num;

    for i=1:l_dn-1
        aux_num=0;
        i_aux=N_d;
        flag=1;
        while flag==1
           if i_aux==N_d
                if aux(i_aux)>0
                   flag=0; 
                else
                    i_aux=i_aux-1;
                end
           elseif aux(i_aux,1)-aux(i_aux+1,1)>1
               flag=0;
           else
                i_aux=i_aux-1;
           end    
        end
        aux(i_aux,1)=aux(i_aux,1)-1;
        if(i_aux~=N_d)
            for j=i_aux+1:N_d
                aux(j,1)=aux(j-1,1)-1;
            end
        end
        for j=1:N_d
            aux_num=aux_num+2^aux(j,1);
        end
        range_dn(i+1,1)=aux_num;
    end
end