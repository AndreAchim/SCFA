function sat=asSatPaire(AS,var)
% sat=asSatPaire(AS,var);
% var(1,2)
% sat(1,2) 
% retourne les saturations des variables v1 et v2 (v2>v1) sur leur facteur commun
% % % crit, le critère associé à la paire, peut servir pour en pondérer la moyenne 
% la pondération ne s'est pas avérée une bonne amélioration
inverse=var(1)>var(2);
if inverse
    var=[var(2) var(1)];
end
r=AS.R(var(1),var(2));
pa=find(ismember(AS.Cpaires,var,'rows'));
% crit=AS.Crit(pa);
p=AS.Ppaires(pa);
if r*p<0   % j ou k est orpheline si signe(poids)~=signe(corr.observée)
    sat=[0 0];
    return
end
s1=sqrt(r/p);
s2=r/s1;
if inverse
    sat=[s2,s1];
else
    sat=[s1,s2];
end