function AS=asDistances(AS)
% AS=asDdistances(AS);
Dist=zeros(AS.nv);
for k=1:numel(AS.Crit)
    Dist(AS.Cpaires(k,1),AS.Cpaires(k,2))=AS.Crit(k);
end
Dist(AS.orphelines,:)=[];
Dist(:,AS.orphelines)=[];
AS.Dist=triU(Dist)';