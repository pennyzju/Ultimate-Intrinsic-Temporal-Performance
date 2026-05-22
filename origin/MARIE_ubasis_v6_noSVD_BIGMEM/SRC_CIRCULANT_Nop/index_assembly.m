function [index_medium] = index_assembly(n,L,M,N)

size_index = L*M*N - n^3;

index_medium = zeros(size_index,3, 'int16');

ctr =0;
% x
for mx = n+1 : L
    for my = 1 : n
        for mz = 1 : n
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end
% y
for mx = 1 : n
    for my = n+1 : M
        for mz = 1 : n
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end
% xy
for mx = n+1 : L
    for my = n+1 : M
        for mz = 1 : n
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end
% z
for mx = 1 : n
    for my = 1 : n
        for mz = n+1 : N
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end
% xz
for mx = n+1 : L
    for my = 1 : n
        for mz = n+1 : N
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end
% yz
for mx = 1 : n
    for my = n+1 : M
        for mz = n+1 : N
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end
%xyz
for mx = n+1 : L
    for my = n+1 : M
        for mz = n+1 : N
            
            ctr = ctr + 1;
            
            index_medium(ctr,1) = mx;
            index_medium(ctr,2) = my;
            index_medium(ctr,3) = mz;
                        
        end
    end
end