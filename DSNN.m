clc;
clear all;
close all;
load data.mat; % load data matrix where rows are samples and columns are genes
X_in = data;
[row,col]=size(X_in);
nr = row*col;
X_input = data_preprocessing(X_in); %log-transform data if the range is too high
parfor itr=1:1:10
    val_itr = zeros(1,9);
    val_itr2 = zeros(1,9);
    [row,col]=size(X_input);
    %Stage-1:Using compressive sensing to recover matrix from incomplete matrix
    for j = 1:1:9 
        m = floor((j/10)*nr);
        rng('shuffle');
        idx = randperm(nr,m);
        op_phi =  opPCISampling(nr, idx);      
        pts = spgSetParms('verbosity',0); 
        sigma = 0.000000001;
        tau = pi;
        H = op_DCT(op_phi,row,col);
        y = op_phi(X_input(:),1);
        z = spg_bpdn(H, y,sigma, opts);
        X_reconstructed = idct2(reshape(z,row,col)); 
        
        %Calculate NMSE after Stage-1
        temp1 = norm(data_reverseprocess(X_reconstructed)-X_in,'fro');
        val_itr(j) = (temp1/norm(X_in,'fro')).^2; 
        
        %Stage-2: Denoising
        X=X_reconstructed;
        Xobs = zeros(row,col);
        Xobs(idx)=X_input(idx);
        [row,col]=size(X);
        M1=zeros(row,col);
        M1((Xobs ~= 0))=1;  
        for c=1:col
            pivot=X(:,c);
            z=Xobs(Xobs(:,c)~=0,c);
            pivot(abs(pivot-mean(z))>=0.8*std(z))=0;
            X(:,c)=pivot;
        end
        X(M1==1)=Xobs(M1==1);
        M=zeros(row,col);
        M((X ~= 0))=1;
        Y=M.*X;
        
        %% parameters initialization
        lambda_2=0.01;
        lambda_4 = .00001;
        rho=1.1;
        W1=rand(row,col);
        B1=W1;
        max_iter=40;
        X= Y + ~Y.*rand(row,col);
        for iter=1:max_iter  
            fprintf('iteration=%d\n', iter); 
            svtmat=(X+B1);
            [U, S, V]=svd(svtmat);
            S=wthresh(S,'s',(lambda_2/lambda_4)); 
            W1=U*S*V';
            X= Y + ~Y.*W1;
            B1=X+B1-W1;
        end
        X_reconstructed2 = X;
        
        %Calculate NMSE after Stage-2
        temp2 = norm(X_in - data_reverseprocess(X_reconstructed2),'fro');
        val_itr2(j) = (temp2/norm(X_in,'fro')).^2;
    end
    nmse1(itr,2) = val_itr;
    nmse2(itr,:) = val_itr2;
end

