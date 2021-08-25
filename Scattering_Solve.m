function Solution = Scattering_Solve(omega,Dim,h,BC,Matx,Maty,Matz, Source)
    %be careful about reshape [N M] and Matrix(:) operation
    N = round(Dim(1)/h);%num of x dim grid points
    M = round(Dim(2)/h);%num of y dim grid points
    temp = Source.';
    Source = temp(:);
    [pmlx, pmly, pmlz] = PML(Dim,h,BC);
    Maxwell = Scatt_Maxwell_Operator_Construct(omega,Dim,h,BC,(Matx.*pmlx).',(Maty.*pmly).',(Matz.*pmlz).');
    Solution = reshape(Maxwell\Source,[N,M]).';
    %Solution = reshape(pcg(Maxwell,Source),[N,M]).';
end

