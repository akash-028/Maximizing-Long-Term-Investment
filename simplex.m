function [phase,basicvar] = simplex(A,b,c,basicvar)
[m,n]          = size(A);
phase          = zeros(m+1,n+1);
phase(1:m,1:n) = A;
phase(m+1,1:n) = c(:);
phase(1:m,end) = b(:);
keep_running = true;
while keep_running
    if any(phase(end,1:n)<0)    % check if there is negative cost coeff.
        [~,J] = min(phase(end,1:n)); % yes, find the most negative
        % now check if corresponding column is unbounded
        if all(phase(1:m,J)<=0) 
          error('problem unbounded. All entries <= 0 in column %d',J);
        % do row operations to make all entries in the column 0 
        % except pivot
        else
            pivot_row = 0;
            min_found = inf;
            for i = 1:m
                if phase(i,J)>0
                    tmp = phase(i,end)/phase(i,J);
                    if tmp < min_found
                        min_found = tmp;
                        pivot_row = i;
                    end
                end
            end
            basicvar(pivot_row)=J;
            % normalize
            phase(pivot_row,:) = phase(pivot_row,:)/phase(pivot_row,J); 
            % now make all entries in J column zero.
            for i=1:m+1
                if i ~= pivot_row
                    phase(i,:)=phase(i,:)-sign(phase(i,J))*...
                        abs(phase(i,J))*phase(pivot_row,:);
                end
            end
        end
    else
        keep_running=false;       
    end
end