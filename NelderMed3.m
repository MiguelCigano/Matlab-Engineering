clc;
clear all;
close all;

X=[1 0; 0 0.25; 0 0]; %Vertices del simplex inicial, en filas.

n=length(X);
alpha=1;
betha=0.5;
gamma=1.5;
v=2;
Error=0.001
fi=zeros(3,1);

%Generación contornos.
figure; 
[Xplot,Yplot]=meshgrid(-6:0.2:6,-6:0.2:6);
Z=(Xplot.^2)-4*Xplot+(Yplot.^2)-Yplot-(Xplot.*Yplot);
contour(Xplot,Yplot,Z, 'LineWidth', 1.5);
hold on;

Q1=30; %Condicion inicial del error.

while Error<=Q1

    %Generación gráfica de los vertices.    
    for i=1:v+1
        plot(X(i,1),X(i,2),'d','LineWidth',2,...
        'MarkerEdgeColor','b',... 
        'MarkerFaceColor', 'b',...
        'MarkerSize',6)
        hold on;
    end 

    %Generación gráfica de las lineas que unen los vertices.
    line([X(1,1),X(2,1)],[X(1,2),X(2,2)],'Color','red','LineStyle','--');
    line([X(1,1),X(3,1)],[X(1,2),X(3,2)],'Color','red','LineStyle','--');
    line([X(2,1),X(3,1)],[X(2,2),X(3,2)],'Color','red','LineStyle','--');

    %Evaluación de cada vertice del simplex en la función.
    %Se concidera x1=x y x2=y.
    for i=1:n 
        syms x y;
        f=subs(((x.^2)-4*x+(y.^2)-y-(x*y)),{x,y},{X(i,:)});
        fi(i,:)=f; %fi -> vector columna donde se guardan las func. evaluadas.
    end  

    %Busqueda en "fi" del mayor valor de las funciones evaluadas.
    [Max,IM] = max(fi);
    [I_row, I_col] = ind2sub(size(fi),IM);

    %Vector fila correspondiente al vertice que da el mayor valor evaluado.
    Xh=X(I_row,:);  

    %Busqueda en "fi" del menor valor de las funciones evaluadas.
    [Min,Im] = min(fi);
    [I_row, I_col] = ind2sub(size(fi),Im);

    %Vector fila correspondiente al vertice que da el menor valor evaluado.
    Xl=X(I_row,:);



    %Calculo del valor medio y su vector fila correspondiente
    for i=1:n
         fimed=subs((x),{x},{fi(i,:)});
         if (fimed<Max) && (fimed>Min)
             VVmed=fimed; 
             Vmed=double(VVmed); 
             Imed=i;
             %Vector fila correspondiente al vertice de valor medio.
             Xmed=X(Imed,:);  
         end
    end

    %Calculo de x0=Xo.
    Xo=(1/2)*((Xl)+(Xmed));
    Xr=(1+alpha)*Xo-alpha.*Xh;
    ffXo=subs(((x.^2)-4*x+(y.^2)-y-(x*y)),{x,y},{Xo(:,:)});
    fXo=double(ffXo);
    ffXr=subs(((x.^2)-4*x+(y.^2)-y-(x*y)),{x,y},{Xr(:,:)});
    fXr=double(ffXr);

    if Min<=fXr && Max>fXr
        b=0; 
    end

    %Probamos la expansión
    %Calculo del vector Xe y el valor de la función evaluada con el.
    Xe=gamma*Xr+(1-gamma)*Xo;
    fXeSyms=subs(((x.^2)-4*x+(y.^2)-y-(x*y)),{x,y},{Xe(:,:)});
    fXe=double(fXeSyms);

    if fXe<Min 
        X(IM,:)=Xe;
        for i=1:v+1
            plot(X(i,1),X(i,2),'d','LineWidth',1,...
            'MarkerEdgeColor','b',... 
            'MarkerFaceColor', 'b',...
            'MarkerSize',6)
            hold on
        end 

        line([X(1,1),X(2,1)],[X(1,2),X(2,2)],'Color','red','LineStyle','--');
        line([X(1,1),X(3,1)],[X(1,2),X(3,2)],'Color','red','LineStyle','--');
        line([X(2,1),X(3,1)],[X(2,2),X(3,2)],'Color','red','LineStyle','--');
    end 


    if fXe>Min
        X(IM,:)=Xr;
        for i=1:v+1
            plot(X(i,1),X(i,2),'d','LineWidth',1,...
            'MarkerEdgeColor','b',... 
            'MarkerFaceColor', 'b',...
            'MarkerSize',6)
            hold on
        end 
  
        line([X(1,1),X(2,1)],[X(1,2),X(2,2)],'Color','red','LineStyle','--');
        line([X(1,1),X(3,1)],[X(1,2),X(3,2)],'Color','red','LineStyle','--');
        line([X(2,1),X(3,1)],[X(2,2),X(3,2)],'Color','red','LineStyle','--');
    
    end
    
    if Min<fXr && Max<fXr
        %Contracción
        Xc=(betha*(Xh))+(1-betha)*Xo;
        ffXc=subs(((x.^2)-4*x+(y.^2)-y-(x*y)),{x,y},{Xc(:,:)});
        fXc=double(ffXc);
        X(IM,:)=Xc;
        if fXc<Max && fXc<fXr
            for i=1:v+1
                plot(X(i,1),X(i,2),'d','LineWidth',1,...
                'MarkerEdgeColor','b',... 
                'MarkerFaceColor', 'b',...
                'MarkerSize',6)
                hold on
            end 
            
            line([X(1,1),X(2,1)],[X(1,2),X(2,2)],'Color','red','LineStyle','--');
            line([X(1,1),X(3,1)],[X(1,2),X(3,2)],'Color','red','LineStyle','--');
            line([X(2,1),X(3,1)],[X(2,2),X(3,2)],'Color','red','LineStyle','--');
                  
        end
    end
     
    Q=((((fi(1,:)-fXo).^2 + (fi(2,:)-fXo).^2 +(fi(3,:)-fXo).^2)/3).^0.5);
    Q1 = double(Q);

end

Xpf=(1/3)*((X(1,:))+(X(2,:))+(X(3,:)));
 plot(Xpf(:,1),Xpf(:,2),'d','LineWidth',1,...
'MarkerEdgeColor','r',... 
'MarkerFaceColor', 'r',...
'MarkerSize',6)

