function AS=asCoplanaire(AS)
% AS=asCoplanaire(AS);
% si trois groupes n'occupent qu'un plan, en enlever un
% Mais comment faire avant d'avoir estimé les saturations?
% En estimant ces saturations (sans moyenner sur toutes les paires du groupe)
AS.coplan=[];
AS.GrCoplan=[];
ote=[];
ng=numel(AS.Gr);
if ng<3
    return
end
options=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-6,'TolX',1e-6);
trios=nchoosek(1:ng,3);
cible=1:numel(AS.pertinent);
for tri=trios'  % pour chaque trio de groupes
    nt=0;
    pire=0;
    for i=AS.Gr{tri(1)}
        for j=AS.Gr{tri(2)}
            for k=AS.Gr{tri(3)}
                nt=nt+1;
                melange=[i j k];
                if any(melange<0) % ignorer un groupe à exclure
                    pire=-1; 
                else
                    % si les deux plus grandes corrélations ne sont pas essentiellement égales,
                    % essayer avec nouvelles initialisations
                    mul=[1 .5 2];
                    for e=1:numel(mul)
                        P=sign(AS.GS(:,melange(1:end-1))'*AS.GS(:,melange(end)));  % initialiser à +ou - 1 selon signe de corrélations
                        P=mul(e)*P./sqrt(numel(P));
                        [P,cr]=fminsearch(@(P) asCrit(P,AS.GS,melange,cible),P,options);
                        [~,cor]=asCrit(P,AS.GS,melange,cible);
                        cor=cor';
                        c=sort(abs(cor)); 
                        if c(end)-c(end-1)<.0001 && all(abs(log10(P))<2)
                            break;
                        end
                    end
                    if c(end)-c(end-1)>.0001
                        keyboard  % l'optimisation n'a pas trouvé un vrai creux
                    end
                    if cr>pire
                        pire=cr;
                    end
                end
            end
            if pire>0
                pire=pire*(AS.N-1);
                if pire<chi2inv(.95.^(1/nt),1)  % les groupes dans tri semblent coplanaires
                    % caculer les saturations à partir des deux premières variables de groupe
                    s=zeros(3,2);
                    v=zeros(3,2);
                    for g=1:3
                        v(g,:)=AS.Gr{tri(g)}(1:2);  % les deux premières variabes de chaque groupe
                        s(g,:)=asSatPaire(AS,v(g,:));  % leurs saturations
                    end
                    C=zeros(3,3);
                    for g1=1:2
                        for g2=(g1+1):3  % pour les trois paires degroupes
                            co=0;
                            for j=1:2 % pour les deux variables du premier groupe
                                for k=1:2  % deux variables du deuxième groupe corrélées à celles du prem
                                    co=co+AS.GS(:,v(g1,j))'*AS.GS(:,v(g2,k))/(s(g1,j)*s(g2,k));
                                end
                            end
                            C(g1,g2)=abs(co)/4;
                        end
                    end
                    C=triU(C);
                    mi=min(C); % corrélation de la paire la plus orthogonale
                    % les positions dans C sont 12, 13 et 23, les rangs manquant sont 3, 2 et 1 d'où:
                    f=4-find(C==mi,1);
                    f=tri(f);
                    ote=[ote f];
                    AS.GrCoplan{end+1}=AS.Gr{f};
                    AS.Gr{f}=-AS.Gr{f};
                    % en garder trace dans AS.coplan
                    AS.coplan=[AS.coplan;tri'];
                end
            end
        end
    end
end
for k=1:numel(ote)
AS.reste=[AS.reste -AS.Gr{ote(k)}];
end
AS.Gr(ote)=[];
