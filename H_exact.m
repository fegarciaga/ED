function E_val=H_exact(N_up, N_dn, Lx, Ly, Lz, tx, ty, tz, t2x, t2y, t2z, kx, ky, kz, U, car)
%% Makes indexation
    [range_u, range_d, l_u, l_d] = Index(N_up, N_dn, Lx, Ly, Lz);
%% Create array with index and values to insert in sparse matrix
    l=l_u*l_d;
    H_nonull=zeros(l*(car*(N_up+N_dn)+1), 3);
    M=Lx*Ly*Lz;
%% Fill matrix
    kx=sqrt(-1)*pi*kx;
    ky=sqrt(-1)*pi*ky;
    kz=sqrt(-1)*pi*kz;
    count=1;
    for i_up=1:l_u
        for i_dn=1:l_d
            % check for every possible basis state
            r=0;
            aux=0;
            for mz=1:Lz
                for iy=1:Ly
                    for jx=1:Lx
                        M_u=bitget(range_u(i_up),M:-1:1);
                        M_d=bitget(range_d(i_dn),M:-1:1);
                        %check in all possible sites
                        if bitand(range_u(i_up), 2^(r))~=0
                            %checks if there are electrons in the given
                            %place
                            if bitand(range_d(i_dn), 2^(r))~=0
                                % If there is a double occupancy it sums to
                                % the potential energy term
                                aux=aux+U;
                            end
                            if Lx~=1
                                % Search for hopping terms in the x axis
                                H_nonull(count:count+1,:)=near_hop(tx, kx, r, i_up, i_dn, range_u, jx, Lx, l_d, 1, 1, i_up, M, M_u);
                                count=count+2;
                            end

                            if Ly~=1
                                H_nonull(count:count+1,:)=near_hop(ty, ky, r, i_up, i_dn, range_u, iy, Ly, l_d, Lx, 1, i_up, M, M_u);
                                count=count+2;
                            end

                            if Lz~=1
                                H_nonull(count:count+1,:)=near_hop(tz, kz, r, i_up, i_dn, range_u, mz, Lz, l_d, Lx*Ly, 1, i_up, M, M_u);
                                count=count+2;
                            end

                            if Ly>3
                                if Lz==1
                                    if Lx==1
                                        %Search for next nearest neighbor
                                        %hopping in the 1D case
                                        H_nonull(count:count+1,:)=next_near_hop1d(t2y, ky, r, i_up, i_dn, range_u, iy, Ly, l_d, Lx, 1, i_up, M, M_u);
                                        count=count+2;
                                    end
                                end
                            end

                            if Lx>3
                                if Lz==1
                                    if Ly==1
                                        H_nonull(count:count+1,:)=next_near_hop1d(t2x, kx, r, i_up, i_dn, range_u, jx, Lx, l_d, 1, 1, i_up, M, M_u);
                                        count=count+2;
                                    end
                                end
                            end

                            if Lz>3
                                if Ly==1
                                    if Lx==1
                                        H_nonull(count:count+1,:)=next_near_hop1d(t2z, kz, r, i_up, i_dn, range_u, mz, Lz, l_d, Lx*Ly, 1, i_up, M, M_u);
                                        count=count+2;
                                    end
                                end
                            end

                            if Lx>1
                                if Ly>1
                                    % Search for next nearest neighbor
                                    % hopping in the 2D case
                                    H_nonull(count:count+3,:)=next_near_hop2d(t2z, kx, ky, r, i_up, i_dn, range_u, jx, iy, Lx, Ly, l_d, 1, Lx, 1, i_up, M, M_u);
                                    count=count+4;
                                end
                            end
                            
                            if Ly>1
                                if Lz>1
                                    H_nonull(count:count+3,:)=next_near_hop2d(t2x, ky, kz, r, i_up, i_dn, range_u, iy, mz, Ly, Lz, l_d, Lx, Lx*Ly, 1, i_up, M, M_u);
                                    count=count+4;
                                end
                            end

                            if Lx>1
                                if Lz>1
                                    H_nonull(count:count+3,:)=next_near_hop2d(t2y, kx, kz, r, i_up, i_dn, range_u, jx, mz, Lx, Lz, l_d, 1, Lx*Ly, 1, i_up, M, M_u);
                                    count=count+4;
                                end
                            end
                        end
                        % Does the same on the dn sector
                        if bitand(range_d(i_dn), 2^(r))~=0
                            if Lx~=1
                                H_nonull(count:count+1,:)=near_hop(tx, kx, r, i_up, i_dn, range_d, jx, Lx, l_d, 1, 0, i_dn, M, M_d);
                                count=count+2;
                            end

                            if Ly~=1
                                H_nonull(count:count+1,:)=near_hop(ty, ky, r, i_up, i_dn, range_d, iy, Ly, l_d, Lx, 0, i_dn, M, M_d);
                                count=count+2;
                            end

                            if Lz~=1
                                H_nonull(count:count+1,:)=near_hop(tz, kz, r, i_up, i_dn, range_d, mz, Lz, l_d, Lx*Ly, 0, i_dn, M, M_d);
                                count=count+2;
                            end

                            if Ly>3
                                if Lz==1
                                    if Lx==1
                                        H_nonull(count:count+1,:)=next_near_hop1d(t2y, ky, r, i_up, i_dn, range_d, iy, Ly, l_d, Lx, 0, i_dn, M, M_d);
                                        count=count+2;
                                    end
                                end
                            end

                            if Lx>3
                                if Lz==1
                                    if Ly==1
                                        H_nonull(count:count+1,:)=next_near_hop1d(t2x, kx, r, i_up, i_dn, range_d, jx, Lx, l_d, 1, 0, i_dn, M, M_d);
                                        count=count+2;
                                    end
                                end
                            end

                            if Lz>3
                                if Ly==1
                                    if Lx==1
                                        H_nonull(count:count+1,:)=next_near_hop1d(t2z, kz, r, i_up, i_dn, range_d, mz, Lz, l_d, Lx*Ly, 0, i_dn, M, M_d);
                                        count=count+2;
                                    end
                                end
                            end

                            if Lx>1
                                if Ly>1
                                    H_nonull(count:count+3,:)=next_near_hop2d(t2z, kx, ky, r, i_up, i_dn, range_d, jx, iy, Lx, Ly, l_d, 1, Lx, 0, i_dn, M, M_d);
                                    count=count+4;
                                end
                            end

                            if Ly>1
                                if Lz>1
                                    H_nonull(count:count+3,:)=next_near_hop2d(t2x, ky, kz, r, i_up, i_dn, range_d, iy, mz, Ly, Lz, l_d, Lx, Lx*Ly, 0, i_dn, M, M_d);
                                    count=count+4;
                                end
                            end

                            if Lx>1
                                if Lz>1
                                    H_nonull(count:count+3,:)=next_near_hop2d(t2x, ky, kz, r, i_up, i_dn, range_d, jx, mz, Lx, Lz, l_d, 1, Lx*Ly, 0, i_dn, M, M_d);
                                    count=count+4;
                                end
                            end
                        end
                        r=r+1;
                    end
                end
            end
            H_nonull(count, 1)=(i_up-1)*l_d+i_dn;
            H_nonull(count, 2)=H_nonull(count, 1);
            H_nonull(count, 3)=aux;
            count=count+1;
        end
    end
    %Selects non null terms
    b = H_nonull(any(H_nonull,2),:);
    % Creates sparse matrix
    H_ex=sparse(b(:,1), b(:,2), b(:,3), l, l);
    %display(H_ex);
    % Solve for lowest eigen-value
    E_val=eigs(H_ex,1, 'sa');
    %display(E_vec);
end