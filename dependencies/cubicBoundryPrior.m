function mask = cubicBoundryPrior(N, Lx)
mask=zeros(N, N, N);
mid=ceil(N/2);
for i=1:N
    for j=1:N
        for k=1:N
            if(i>(mid-Lx) && i<(mid+Lx) && ...
                    j>(mid-Lx) && j<(mid+Lx) && ...
                    k>(mid-Lx) && k<(mid+Lx))
                mask(i, j, k)=1.0;
            end
        end
    end
end
end