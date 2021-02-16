function [Y] = generalProjection(X, epsilon)

% this function projects each column x of a matrix X on the set defined by
%  {y in R^d s.t. y>=-epsilon, sum_i y_i=1}
% The solution is given, column-wise, by y_i=max(-epsilon, x_i-mu) for all
% i, where mu is given by solving sum_i max(-epsilon, x_i-mu)=1
%
% ****** Input ****** 
% X       : d-by-r matrix
% epsilon : scalar or r-by-1 vector, generally positive (though not compulsory)
% r       : factorization rank r 
% 
% ****** Output ****** 
% Y       : the projected matrix

% loop over each column of X
if nargin <= 1
    epsilon = 0;
end
if length(epsilon) == 1
    epsilon = epsilon*ones(size(X,2) ,1);
end
for i=1:size(X,2) 
    % sort each column of the input matrix
    x=X(:,i);
    x_bis=sort(x);
    
    len=length(x);
    index_min=1;
    index_max=len;

    % mu s.t. y_i < mu-epsilon, forall i
    mu_max=x_bis(len)+epsilon(i);
    % mu s.t. y_i > mu-epsilon, forall i
    mu_min=x_bis(1)+epsilon(i);

    somme_min=sum(x)-len*mu_min;
    somme_max=-len*epsilon;
   
    if (somme_min < 1)
        mu=(sum(x)-1)/len;
        y=max(repmat(-epsilon(i), len, 1), x-mu);
        Y(:,i)=y;
    else
        % apply a dichotomy to find the optimal value of mu 
        stop=0;
        while stop==0
            curr_ind=round((index_min+index_max)/2);
            mu=x_bis(curr_ind)+epsilon(i);
            y=max(repmat(-epsilon(i), len, 1), x-mu);
            somme=sum(y);
            if (somme < 1)
                index_max=curr_ind;
            elseif (somme > 1)
                index_min=curr_ind;
            else
               Y(:,i)=y;
               stop=1;
            end

            if index_max == index_min +1 
                stop=1;
            end

        end

        mu_inf=x_bis(index_min)+epsilon(i);
        mu_sup=x_bis(index_max)+epsilon(i);
        obj_inf=sum(max(repmat(-epsilon(i), len, 1), x-mu_inf));
        obj_sup=sum(max(repmat(-epsilon(i), len, 1), x-mu_sup));

        slope=(obj_sup-obj_inf)/(mu_sup-mu_inf);
        mu_opt=(1-obj_inf)/slope+mu_inf;
        
        % compute the corresponding column of Y
        y=max(repmat(-epsilon(i), len, 1), x-mu_opt);
        Y(:,i)=y;
    end

end

end

