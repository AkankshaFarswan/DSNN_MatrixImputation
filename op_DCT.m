function op = op_DCT(phi,n1,n2)

op = @(z,mode) opRWT_CSMRI_2D_intrnl(phi,n1,n2,z,mode);


function y = opRWT_CSMRI_2D_intrnl(phi,n1,n2,z,mode)

if mode==1
    %y = kron(dctmtx(n2)',dctmtx(n1)')*z;
    y = idct2(reshape(z,n1,n2));
    y = phi(y(:),1);
else
    temp = phi(z,2);
    %y = kron(dctmtx(n2),dctmtx(n1))*temp;
    y = dct2(reshape(temp,n1,n2));
    y = y(:);
end
