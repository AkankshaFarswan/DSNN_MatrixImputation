%clc;
%clear all;
%close all;
% load data matrix where rows are samples and columns are genes
load data.mat; 
X_in = data;
[row,col]=size(X_in);
nr = row*col;
%log-transform data using data_preprocessing if the range is too high
X_input = data_preprocessing(X_in);
itrn=30;
parfor itr=1:1:itrn        %(Using parfor for parallel computation)
    val_itr = zeros(1,9);
    val_itr2 = zeros(1,9);
    [row,col]=size(X_input);
    %Stage-1:Using compressive sensing to recover matrix from incomplete matrix
    for j = 1:1:9  % "j" determines the percentage of input observed 1 is for 10% observed input and 9 for 90% of observed input
        m = floor((j/10)*nr);
        rng('shuffle');
        idx = randperm(nr,m);  % randomly chosing the indices for introducing missing values
        op_phi =  opPCISampling(nr, idx);    %Creating phi object using the function opPCISampling  
        pts = spgSetParms('verbosity',0); 
        sigma = 0.000000001;
        tau = pi;
        H = op_DCT(op_phi,row,col);          %Creating DCT function handle using op_DCT function.
        y = op_phi(X_input(:),1);            %Input incomplete vector "y" with randomly introduced missing values 
        z = spg_bpdn(H, y,sigma, opts);      %Recovering complete matrix using "spgl" solver       
        X_reconstructed = idct2(reshape(z,row,col)); %Perform Inverse DCt transform on the output to get the imputed data matrix
        
        %Calculate NMSE after Stage-1
        %if data is log tranformed initially, use data_reversedataprocess() to convert back data to original values 
        X_rec1 = data_reverseprocess(X_reconstructed);
        temp1 = norm(X_rec1-X_in,'fro'); 
        val_itr(j) = (temp1/norm(X_in,'fro')).^2; 
        
        %Matrix from Stage-1 is considered a noisy version of original matrix.
        %Therefore, denoising is done to remove noise from the output matrix recovered from Stage-1
        %Stage-2: Denoising using nuclear norm minimization
        X=X_reconstructed;
        Xobs = zeros(row,col);
        Xobs(idx)=X_input(idx);
        [row,col]=size(X);
        M1=zeros(row,col);
        M1((Xobs ~= 0))=1;  
        %Checking for entries that are highly noisy and filling them with random entries.
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
        %if data is log tranformed initially, use data_reversedataprocess() to convert back data to original values 
        X_rec2 = data_reverseprocess(X_reconstructed2);
        temp2 = norm(X_rec2-X_in,'fro');
        val_itr2(j) = (temp2/norm(X_in,'fro')).^2;
    end
    nmse1(itr,:) = val_itr;
    nmse2(itr,:) = val_itr2;
end
avg_nmse1 = mean(nmse1); %(Averaged NMSE after stage-1 over 30 iterations)
avg_nmse2 = mean(nmse2); %(Averaged NMSE after stage-1 + stage-2 over 30 iterations)
%Generate semilogy plot
x=10:10:90
semilogy(x,avg_nmse1,x,avg_nmse2)

