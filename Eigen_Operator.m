function Eigen_Maxwell = Eigen_Operator(Dim,h,BC,Matx,Maty,Matz)
    [pmlx, pmly, pmlz] = PML(Dim,h,BC);
    Eigen_Maxwell = Eigen_Maxwell_Operator_Construct(Dim,h,BC,(Matx.*pmlx).',(Maty.*pmly).',(Matz.*pmlz).');
end

