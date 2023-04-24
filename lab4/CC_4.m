g = 'circleg';
bc = 'circleb1';
c=1; a=0; f=1;

[p,e,t] = initmesh(g, 'Hmax', 0.2);

u = -(p(1,:).^2+p(2,:).^2-1)/4;
uh = assempde(bc,p,e,t,c,a,f);
uh1 = -(  p(1,:).^2 - p(1,:)./2 + p(2,:).^2 - p(2,:)./2  - 1/2 )/2.* (p(1,:).^2 + p(2,:).^2 - 1);

t_ = t; t_(4,:)=[];

%ищу самостоятельно производную
points = length(p);
traingles = length(t_);
edges_of_traingle = length(t_(:,1));

normgraddifferror = zeros(1,points);
%normgraddiffuh = zeros(1,points);
graduh = zeros([2 points]);

error = u-uh1';
%uh = uh1';

%в каждой точке
for i =1:1:points
    
    %ищу соседние точки
    for col = 1:1:traingles
    
      %найти текущую точку в треугольнике
      if ismember(i, t_(:,col))

        for k = 1:1:edges_of_traingle
            j = t_(k,col); %сосед

            if i ~= j
                eij = [p(1,j)-p(1,i) p(2,j)-p(2,i)];
                tga = (error(j)-error(i))/norm(eij);
                %tga = (uh(j)-uh(i))/norm(eij);
                eij = eij./norm(eij);
            
                grad = eij.*tga;
                norm_grad = norm(grad);
            
                %запомнить
                 if normgraddifferror(i) < norm_grad
                    normgraddifferror(i) = norm_grad;
                 %if normgraddiffuh(i) < norm_grad
                    %normgraddiffuh(i) = norm_grad;
                    graduh(1, i) = grad(1);
                    graduh(2, i) = grad(2);
                
                end
            
            end
        end
      end
    
    end

end

%Нашли градиент ошибки. Начинаем его осреднять

%Найду площадь всех треугольников
traingles_measures = zeros(1,traingles);

for tr = 1:1:traingles
    
    %tr = [p1,p2,p3] - три точки треугольника. Площадь по формуле Герона

    a = sqrt((p(1, t_(1, tr)) - p(1, t_(2, tr)))*(p(1, t_(1, tr)) - p(1, t_(2, tr))) + (p(2, t_(1, tr)) - p(2, t_(2, tr)))*(p(2, t_(1, tr)) - p(2, t_(2, tr))));
    b = sqrt((p(1, t_(2, tr)) - p(1, t_(3, tr)))*(p(1, t_(2, tr)) - p(1, t_(3, tr))) + (p(2, t_(2, tr)) - p(2, t_(3, tr)))*(p(2, t_(2, tr)) - p(2, t_(3, tr))));
    c = sqrt((p(1, t_(3, tr)) - p(1, t_(1, tr)))*(p(1, t_(3, tr)) - p(1, t_(1, tr))) + (p(2, t_(3, tr)) - p(2, t_(1, tr)))*(p(2, t_(3, tr)) - p(2, t_(1, tr))));

    half_perimeter = (a+b+c)/2;
    
    traingles_measures(tr) = sqrt(half_perimeter*(half_perimeter-a)*(half_perimeter-b)*(half_perimeter-c));
end

%осредняю градиенты
%new_normgraddiffuh = zeros([1 points]);
new_normgraddifferror = zeros([1 points]);
neighbour_traingles_i = [];
Tij = [];
Omegai = 0;
gi = zeros([2 points]);

size = 1;
for i=1:1:points
    
    neighbour_traingles_i = [];
    Tij = [];
    Omegai = 0;
    
    %ищу соседние треугольники
    size = 1;
    for j = 1:1:traingles
        
      %найти текущую точку в треугольнике
      if ismember(i, t_(:,j))
          neighbour_traingles_i(1, size) =  t_(1,j);
          neighbour_traingles_i(2, size) =  t_(2,j);
          neighbour_traingles_i(3, size) =  t_(3,j);
          Tij(size) = traingles_measures(j);
          
          Omegai = Omegai +  Tij(size);
          size = size + 1;
      end
      
    end
    
    %осредняю градиенты в этой точке
    for j = 1:1:size-1
        gi(1, i) = Tij(j)/Omegai * graduh(1, neighbour_traingles_i(1, j));
        gi(2, i) = Tij(j)/Omegai * graduh(1, neighbour_traingles_i(1, j));
    end
    
    %осредненный градиент
    %new_normgraddiffuh(i) = norm(gi);
    new_normgraddifferror(i) = norm(gi);
end


ei = zeros([2 points]);
%индикатор!
for i = 1:1:points
    ei(1,i) = norm(gi(:,i) - graduh(:,i));
    ei(2,i) = i; %номер точки
end

%Сортировка по максимальной ошибке
for i = 1:1:points-1
    for j = 1:1:points-1
        if ei(1,i)>ei(1, i+1)
            temp = ei(1,i);
            ei(1,i) = ei(1,i+1);
            ei(1,i+1) = temp;
            
            temp = ei(2,i);
            ei(2,i) = ei(2,i+1);
            ei(2,i+1) = temp;
        end
    end
end

% 30% наибольшей ошибки
lim = round(0.3*points);

figure
pdemesh(p,e,t);
hold all;
pdesurf(p,t, normgraddifferror');
colormap('hsv');
xlabel('x'), ylabel('y'), zlabel('||grad(u - uh)||');
colorbar;


%Поверхность разницы численного градиента и осреднения
figure
pdemesh(p,e,t);
hold all;
colormap('hsv');
xlabel('x'), ylabel('y'), zlabel('||grad(uh) - G(grad(uh))||');
colorbar;

for i = 1:1:lim
    plot3( p(1, ei(2, points-i+1)), p(2, ei(2, points-i+1)), 1 ,'ro');
end






%Чтобы отмечать н
figure
subplot(2,1,1)
pdemesh(p,e,t);
hold all;
colormap('hsv');
xlabel('x'), ylabel('y');
title('||G(grad(uh)) - grad(uh)||L2(Tij)');
%colorbar;

subplot(2,1,2)
pdemesh(p,e,t);
hold all;
colormap('hsv');
xlabel('x'), ylabel('y');
title('||grad(u) - grad(uh)||L2(Tij)');
%colorbar;


%Поверхность реальной ошибки
figure
pdemesh(p,e,t);
hold all;
pdesurf(p,t, abs(u'-uh));
colormap('hsv');
xlabel('x'), ylabel('y'), zlabel('||u - uh||');
colorbar;


