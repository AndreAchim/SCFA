function AS=asSaturations(AS)
% AS=asSaturations(AS);
ng=numel(AS.Gr);
AS.Var=zeros(ng,1);  % variable qui pourra représenter chaque facteur
AS.GrDe=zeros(AS.nv,1);  % groupe de chaque variable unifactorielle
for gr=1:ng
    var=sort(AS.Gr{gr});
    AS.Gr{gr}=var;  % s'assurer que les variables sont en ordre croissant pour les retrouver dans AS.Cpaires
    n=numel(var);
    satur=zeros(n-1,n);   % pour pouvoir examiner l'homogénéité des n-1 estimations de chaque variable
    % crit=zeros(n-1,n);
    rg=zeros(n,1);  % les rangs où écrire les saturations de chaque variable
    for j=1:n-1
        for k=j+1:n  % pour chaque paire des variables du facteur
            if var(j)<0 || var(k)<0
                break;   % ignorer une variable de rang rendu négatif (rouvée orpheline)
            end
            % [sat,cr]=asSatPaire(AS,var([j k]));
            sat=asSatPaire(AS,var([j k]));
            if all(sat==0)   % asSatPaire retourne [0,0] si signe(poids)~=signe(corr.observée) (j ou k est orpheline)
                keyboard   % la suite sous cette condition reste à valider
                f=[j k];
                s=sum(AS.R(:,f).^2)-1;
                f=f(1+(s(2)<s(1)));
                AS=declareOrpheline(AS,f,gr);
                var(var==f)=-var(var==f);
                % satur(:,f)=[];   % <<<<< pas une bonne idée d'enlever une colonne
            else
                if satur(1,j)~=0 && sign(sat(1))~=sign(satur(1,j))  % assurer polarités constantes
                    sat=-sat;
                end
                rg([j k])=rg([j k])+1;
                satur(rg(j),j)=sat(1);
                satur(rg(k),k)=sat(2);
                % crit(rg(j),j)=cr;
                % crit(rg(k),k)=cr;
            end
        end
    end
    if any(var<0)
        AS.pertinent(AS.pertinent==-var(var<0))=[];
    end
    AS.satur{gr}=satur;
    if size(satur,1)>1
        % po=1./crit;
        % po=po./sum(po);
        % satur=satur.*po;
        % satur=sum(satur);
        satur=mean(satur);
    end
    if sum(satur)<0
        satur=-satur;
    end
    AS.Fct(AS.Gr{gr},gr)=satur;
    f=find(satur==max(abs(satur)),1);
    f=var(f);
    AS.Var(gr)=f;
    AS.GrDe(f)=gr;
end