function boot = be_perfrom_boot(XB,XF)

rng(1)
[tmp] = bootstrp(1e3,@(x) [mean(x) var(x)],XB);
boot.G1 = tmp(:,1);
boot.V1 = tmp(:,2);
rng(2)
[tmp] = bootstrp(1e3,@(x) [mean(x) var(x)],XF);
boot.G2 = tmp(:,1);
boot.V2 = tmp(:,2);

boot.n1 = repmat(length(XB),1,1e3);
boot.n2 = repmat(length(XF),1,1e3);
% 
boot.CI1 = prctile(boot.G1,[100*.05/2,100*(1-.05/2)]);
boot.CI2 = prctile(boot.G2,[100*.05/2,100*(1-.05/2)]);

end
