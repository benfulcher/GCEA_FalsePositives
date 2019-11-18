function sumSqResid = SpatialLagForm(y,distMat,rho,d0,doPairwise)

W = exp(-distMat/d0);
W(logical(eye(size(W)))) = 0;

if doPairwise
    Y = y*y';
    residual = Y - rho*W;
    residual = residual(logical(triu(size(residual),+1)));
    sumSqResid = sum(residual.^2);
else
    residual = y - rho*W*y;
    sumSqResid = sum(residual.^2);
end

end
