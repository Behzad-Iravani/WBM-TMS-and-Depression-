function plot_colorbars2(mdl)


index = ~cellfun(@(x) contains('zage',x) | contains('sex_F',x) | contains('sex_M',x) | contains('(Intercept)',x), mdl.name{1});
sig  = sign(mdl.CI(:,index));
sig  = sig(1,:).*sig(2,:);

values = nanmean(mdl.coeff(:,index));


col = [[0 77 129;237 137 93]./255;ones(1,3)./2];
col =col([2 3 1],:);
hold on
b = zeros(1,length(values));
for i = 1:length(values)
if values(i) >0 && sig(i)>0
    b(i) = 1;
  
elseif values(i)<0  && sig(i)>0
     b(i) = -1;
% else
%     b.FaceColor = col(3,:);  
% 
end

end
imagesc(b)
c = 0;
for i =find(index)
    c = c + 1;
    if b(c) ==1
        text(c,1,sprintf('CI = [%1.2f %1.2f]', mdl.CI(1,i),mdl.CI(2,i)),...
    HorizontalAlignment="center", Rotation= 90 , Color='w')

    else
        text(c,1,sprintf('CI = [%1.2f %1.2f]', mdl.CI(1,i),mdl.CI(2,i)),...
    HorizontalAlignment="center", Rotation= 90 )
    end

end
axis tight
line([[1.5:(sum(index)+.5)]'*ones(1,2)]', repmat(ylim(),sum(index),1)' , 'color','k')
set(gca, 'ytick','', 'xtick', 1:sum(index))
colormap(col)
caxis([-1 1])
end