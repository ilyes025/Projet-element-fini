
clear all

clc


E=200*10^9; I=29*10^-6; h=0.15; L=1; F0=60*10^3; q0=2.4*10^3; nh=6;

%E=1; I=1; h=0.15; L=1; F0=60; q0=100; nh=6;
%****************************************************************** TP2 **************************************************************************************************
%Réponse a l'exercice 1

  %Derivé seconde de la fonction psi
    function y=psisec(i,xb,he) 
      switch (i)
        case(1)
          y=-6/he^2+12*xb./he^3;
        case(2)
          y=-4/he+6*xb./he^2;
        case(3)
          y=6/he^2-12*xb./he^3;
        case(4)
          y=-2/he+6*xb./he^2;
        endswitch
    endfunction
  
  
  %La matrice de rigidité élémentaire
    function [ke] =kele(x1,x2,E,I)
      he=x2-x1;
      %he=x;
      k=zeros(4,4);
      xb=0:1/100:he; 
      for i=1:4
        for j=1:4
          k(i,j)=E.*I.*trapz( xb,psisec(i,xb,he).*psisec(j,xb,he) );
        endfor
      endfor
      ke=k;
    endfunction
  
  
%Réponse a l'exercice 2

  %La fonction psi
    function y=psi(i,xb,he) 
      switch (i)
        case(1)
          y=1-3*(xb./he).^2+2*(xb./he).^3;
        case(2)
          y=xb.-2/he*xb.^2+xb.^3/he^2;
        case(3)
          y=3*(xb./he).^2-2*(xb./he).^3;
        case(4)
          y=-xb.^2/he+xb.^3/he^2;
        endswitch
    endfunction
    
  %Vecteur force élémentaire
    function [fe]=fele(x1,x2,q0,L)
      % L=1; variable globale !!
      he=x2-x1;
      %he=x;
      f=[];
      xb=0:1/100:he;
      for i=1:4
        
      q=(xb.*q0./L); % charge linéaire
       %=q0 ; % charge constante
       f(i)= trapz( xb,q.*psi(i,xb,he));
       
       endfor
      fe=f';
    endfunction
%L'exécution de programme avec x1=0 et x2=1

ke=kele(0,0.5,E,I)
fe=fele(0,0.5,q0,L)

%********************************************************************** TP3 ****************************************************************************************

%Réponse a l'exercice 1
    %La matrice globale de rigidité 
    %1er maillage
    function k= kele_gene(n,l,E,I)
      k=zeros(2*n+2,2*n+2);
      x=0:l/n:l;
      l=0;
      for i=1:2:(2*n)
        l=l+1;
        bloc=i:i+3;
        mat=k(bloc,bloc)+kele(x(l),x(l+1),E,I);
        k(bloc,bloc)=mat;
      endfor
    endfunction
     %2eme maillage
    function k= kele_gene1(n,l,E,I)
      k=zeros(2*n+2,2*n+2);
      x=zeros(n+1,1);
      for i=1:n+1
        x(i)=l/2*(1-cos((i-1)*pi/n));
      endfor
      l=0; 
      for i=1:2:(2*n)
        l=l+1;
        bloc=i:i+3;
        mat=k(bloc,bloc)+kele(x(l),x(l+1),E,I);
        k(bloc,bloc)=mat;
      endfor
    endfunction
    

ka=kele_gene(nh,1,E,I)
    %Le vecteur force global
    %1er maillage
    function f= fele_gene(n,l,q0,L)
      f=zeros(2*n+2,1);
      x=0:l/n:l;
      l=0;
      for i=1:2:2*n
        l=l+1;
        bloc=i:i+3;
        vect=f(bloc,1)+fele(x(l),x(l+1),q0,L);
        f(bloc,1)=vect;
      endfor
    endfunction
    %2eme maillage
    function f= fele_gene1(n,l,q0,L)
      f=zeros(2*n+2,1);
      x=zeros(n+1,1);
      for i=1:n+1
        x(i)=L/2*(1-cos((i-1)*pi/n));
      endfor
      l=0;
      for i=1:2:2*n
        l=l+1;
        bloc=i:i+3;
        vect=f(bloc,1)+fele(x(l),x(l+1),q0,L);
        f(bloc,1)=vect;
      endfor
    endfunction
fa=fele_gene(nh,1,q0,L)

 %Réponse a l'exercice 2
  U=zeros(2*nh+2,2);  % Les vecteurs U et Q ont deux colonnes
  Q=zeros(2*nh+2,2);  % la première définit l'état (connu / inconnu) (1/0) et la seconde contient les valeurs
  %pour U
  U(1,1)=1;
  U(2,1)=1;
  %pour Q
  for i=3:2*nh+2
    Q(i,1)= 1;
  endfor
  Q(2*nh+1,2)=F0;
  % affichage des deux vecteurs
  Q
  U
  % Pour l'ex 2 on peut utiliser une seule colonne dans les vecteurs si notre valeur est inconnue on peut mettre NaN.

 
%******************************************************************** TP4*******************************************************************************************

  %1er maillage
  x=0:1/nh:1;
  % fonction de la 2eme maillage
  function x1=maill_2(nh,L)
    x1=zeros(nh+1,1);
    for i=1:nh+1
      x1(i)=L/2*(1-cos((i-1)*pi/nh));
    endfor
  endfunction
  
  
  %deplacement dans les deux maillage
    function [UT,U]=deplacement_maill_1(nh,q0,L,E,I,Q)
      F=Q(3:2*nh+2,2).+fele_gene(nh,1,q0,L)(3:2*nh+2);
      K1=kele_gene(nh,1,E,I)(3:2*nh+2,3:2*nh+2);
      U1=inv(K1)*F;
      UT=zeros(2*nh+2,1);
      UT(3:2*nh+2)=U1;
      U=UT(1:2:2*nh+2);
    endfunction
    function [UT,U]=deplacement_maill_2(nh,q0,L,E,I,Q)
      F=Q(3:2*nh+2,2).+fele_gene1(nh,1,q0,L)(3:2*nh+2);
      K1=kele_gene1(nh,1,E,I)(3:2*nh+2,3:2*nh+2);
      U1=inv(K1)*F;
      UT=zeros(2*nh+2,1);
      UT(3:2*nh+2)=U1;
      U=UT(1:2:2*nh+2);
    endfunction
    [UT1,U1]=deplacement_maill_1(nh,q0,L,E,I,Q);
    [UT2,U2]=deplacement_maill_2(nh,q0,L,E,I,Q);
  
    % déformation théorique
    function v =deformation_theo(x,nh,q0,L,E,I,F0)
      v =(F0/(6*E*I)*x.^2).*(3*L-x)+q0*L^4/(120*E*I)*(20*x.^2/L^2-10*x.^3/L^3+x.^5/L^5);
    endfunction

    x=0:1/nh:1;
    v=deformation_theo(x,nh,q0,L,E,I,F0);
    v1=deformation_theo(maill_2(nh,L),nh,q0,L,E,I,F0);
    subplot(2,2,2)
    plot(x,U2,x,v1);
    title ("déformation avec 2eme maillage");
    subplot(2,2,1)
    plot(x,U1,x,v);
    title ("déformation avec 1er maillage");
  %détermination de Q pour un element donné
    function q=Q(i,nh,U,E,I,L,q0)      % i cad le ieme element
      x=0:1/nh:1;
      n=zeros(nh,1);
      j=1;
      for k=1:nh
        n(k)=j;
        j=j+2;
      endfor
      q=(kele(x(i),x(i+1),E,I)*U(n(i):n(i)+3,1)).-fele(x(i),x(i+1),q0,L);
    endfunction
  % Mfz_pratique
    function m=Mfz_pratique(nh,U,E,I,L,q0)
      m=zeros(nh+1,1);
      m(1)=-Q(1,nh,U,E,I,L,q0)(2);
      for i=2:nh+1
      m(i)=Q(i-1,nh,U,E,I,L,q0)(4);
      endfor
    endfunction
  % Ty_partique
    function t=Ty_partique(nh,U,E,I,L,q0)
      t=zeros(nh+1,1);
      t(1)=-Q(1,nh,U,E,I,L,q0)(1);
      for i=2:nh+1
      t(i)=Q(i-1,nh,U,E,I,L,q0)(3);
      endfor
    endfunction

  % fonctions théoriques
    function t=Ty_theo(x,I,L,q0,F0)
      t=F0+q0*L/2*(1-(x.^2)./L^2);
    endfunction
    
    function m=Mfz_thep(x,I,L,q0,F0)
      m= F0*(L-x).+q0*L^2/6*(2-3*x./L+x.^3/L^3);
    endfunction
    
  

  tt=Ty_theo(x,I,L,q0,F0);
  tp=Ty_partique(nh,UT1,E,I,L,q0);
  mt=Mfz_thep(x,I,L,q0,F0);
  mp=Mfz_pratique(nh,UT1,E,I,L,q0);
  subplot(2,2,3);
  plot(x,tp,x,tt)
  title ("Ty");
  subplot(2,2,4);
  %axis([0,1,0,70000]);
  plot(x,mp,x,mt)
  title ("Mfz");