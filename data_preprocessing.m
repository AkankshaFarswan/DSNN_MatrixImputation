function [X_process] = data_preprocessing(Y)
% load CLL_SUB_111.mat;
% X_CLL = X+1;
% X_CLL = log10(X_CLL);
% 
% 
% load ALLAML.mat;
% X_ALLAML=X;
% for j=1:size(X_ALLAML,2)
%     X_ALLAML(:,j) = X_ALLAML(:,j)-min(X_ALLAML(:,j));
% end
% 
% load GSE1159.mat;
% X_G =GSE1159gene';
% X_GSE = X_G +1;
% X_GSE = log10(X_GSE);

% for j=1:size(Y,2)
%     indx = find(Y(:,j)~=0);
%     min_val(j) = min(Y(indx,j));
%     X_process(:,j) = Y(:,j)-min_val;
% end
[Yr Yc] = size(Y);
for ir=1:1:Yr
    for ic=1:1:Yc
        if(Y(ir,ic)~=0)
          X_process(ir,ic) = log10(Y(ir,ic)+1);
        end
    end
end
%X_process = log10(Y+1);




