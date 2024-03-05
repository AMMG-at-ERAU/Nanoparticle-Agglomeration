function [ ptr ] = A_clusters( dCNT )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
N = length(dCNT);
ptr = zeros(1,N);
for i=1:N
    r1 = i;
    b = dCNT(1:i-1,i)~=0;
    ptr(i) = -1;
        if ~isempty(find(b, 1))
            c = find(b);
            for j=1:length(c)
                [r2, ptr] = A_find_root(c(j), ptr); %find root of intersect CNT (follow A fast Monte Carlo algorithm for site or bond percolation)
                if r2 ~= r1 % CNT not in the three of CNT "i"
                    if ptr(r1) > ptr(r2)
                        ptr(r2) = ptr(r2) + ptr(r1); %merge number of CNT under the two threes
                        ptr(r1) = r2; %root of i changes to index of the intersect
                        r1 = r2;
                    else
                        ptr(r1) = ptr(r1) + ptr(r2);
                        ptr(r2) = r1;
                    end  
                end
            end  
        end
end

end

