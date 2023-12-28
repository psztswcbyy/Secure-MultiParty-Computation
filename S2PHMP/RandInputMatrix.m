function M = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp)
Ori_M = rand(N,'double') + randi([FirstNumMin,FirstNumMax],N);
Exp_M = 10.^(randi([minEp, maxEp],N));
Sign_M = sign(2*rand(N)-1);
M = Sign_M.*Ori_M.*Exp_M;
end