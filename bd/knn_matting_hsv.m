function alpha = knn_matting_hsv(im, scrib)

    lambda=1000;
    level=2;
    factor=1;
    
    scrib=reshape(scrib(:,:,1),[],1);
    [m, n, d]=size(im);
    val=scrib>0.9;
    map=(scrib<0.1)+val;
    nn=[10;2];
    [a, b]=ind2sub([m n],1:m*n);
    feature=[reshape(cos(im(:,:,1)*2*pi),m*n,1)'*factor;reshape(sin(im(:,:,1)*2*pi),m*n,1)'*factor;reshape(im(:,:,2:3),m*n,2)'/2;[a;b]/sqrt(m*m+n*n)*level+rand(2,m*n)*1e-6];
    now=0;
    for i=1:size(nn,1)
        ind=vl_kdtreequery(vl_kdtreebuild(feature),feature,feature,'NUMNEIGHBORS',nn(i),'MAXNUMCOMPARISONS',nn(i)*2);
        a=reshape(repmat(uint32(1:m*n),nn(i),1),[],1);
        b=reshape(ind,[],1);
        row(now+1:now+m*n*nn(i),:)=[min(a,b) max(a,b)];
        feature(5:6,:)=feature(5:6,:)/100;
        now=now+m*n*nn(i);
    end
    row=unique(row,'rows');
    value=max(1-sum(abs(feature(:,row(:,1))-feature(:,row(:,2))))/6,0);
    A=sparse(double(row(:,1)),double(row(:,2)),value,m*n,m*n);
    A=A+A';
    D=spdiags(sum(A,2),0,n*m,n*m);
    M=D-A+lambda*spdiags(map,0,m*n,m*n);
    L=ichol(M);
    x=pcg(M,lambda*val,[],2000,L,L');
    alpha=reshape(x,m,n);
end