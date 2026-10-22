%% Tp Final PDI
T = readtable('hemorrhage_diagnosis.csv');   % tabla: un renglon por corte
head(T)    % muestra SOLO los primeros renglones

conH = T(T.No_Hemorrhage == 0, :);   % solo los renglones donde SOLO hay hemorragia
conH(1:30 , {'PatientNumber', 'SliceNumber'}) ;  % paciente y corte de los primeros 30 con hemorragia

%% Primer corte con hemorragia de la lista
p = conH.PatientNumber(15);       % numero de paciente
s = conH.SliceNumber(15);         % numero de corte

%% Cargo imagenes
[Ibrain, M] = cargar(p, s);       % ventana de cerebro y mascara de los medicos
Ihueso = imread(fullfile('Patients_CT', sprintf('%03d', p), 'bone', sprintf('%d.jpg', s)));

%% Muestro Imagenes
figure
subplot(1,3,1);
imshow(Ibrain);
title(sprintf('Paciente %d, corte %d', p, s));

subplot(1,3,2);
imshow(Ihueso);
title('Ventana de hueso');

subplot(1,3,3);
imshow(Ibrain)
hold on
visboundaries(M, 'Color', 'r');  
title('Hemorragia marcada');

size(Ibrain)                      
figure
imshow(Ibrain)
hold on
visboundaries(M, 'Color', 'r');    
impixelinfo         

%% segmentacion de hueso y cerebro y hemorragia
[interior, cf, ref, craneo, disco] = segmentar_cerebro(Ibrain);
cerebro = Ibrain;  
cerebro(~interior) = 0;

figure
subplot(2,2,1);
imshow(craneo);
title('1) Hueso (>240)');

subplot(2,2,2);
imshow(disco);
title('2) Relleno');

subplot(2,2,3);
imshow(interior);
title('3) Interior')

subplot(2,2,4);
imshow(cerebro);
hold on
visboundaries(M, 'Color', 'r');
title('4) Cerebro sin cráneo + hemorragia')

% veo que porcentaje de la hemorragia quedo adentro de la mascara del cerebro:
porc = 100 * sum(M(:) & interior(:)) / sum(M(:))

%% Histograma: sangre vs cerebro sano
figure
histogram(cerebro(M & interior), 0:255, 'Normalization','pdf')
hold on
histogram(cerebro(~M & interior), 0:255, 'Normalization','pdf')
legend('Sangre','Cerebro sano');
xlabel('Nivel de gris');

%% Segmento sangre (método umbralizado directo)
cf = medfilt2(cerebro, [5 5]); 
sangre = (cf > 150) & (cf < 215) & interior; % rango de la sangre, solo adentro del cerebro
sangre = bwareaopen(sangre, 50);           % saco puntitos chicos

% Comparo con la máscara de los médicos (Dice: 1 = perfecto, 0 = nada)
dice = 2*sum(sangre(:) & M(:)) / max(sum(sangre(:)) + sum(M(:)), 1)

figure, imshow(Ibrain), hold on
visboundaries(M, 'Color', 'r')        % rojo = médicos
visboundaries(sangre, 'Color', 'g')   % verde = nuestro
title(sprintf('Paciente %d, corte %d | Dice = %.2f', p, s, dice));

%% Veo porque dice me dan tan chico:
nuestro = sum(sangre(:))
medicos = sum(M(:))
acierto = sum(sangre(:) & M(:))

%% Segmento sangre con DELTA y R_OPEN
% DELTA define qué tan más clara debe ser la sangre en comparación con el cerebro sano del mismo corte
% Ropen es el radio de un disco para una apertura morfologica

DELTA  = 35;      
R_OPEN = 2;       
sangre_opt = umbralizar(cf, ref, interior, DELTA, R_OPEN);
dice_opt = 2*sum(sangre_opt(:) & M(:)) / max(sum(sangre_opt(:)) + sum(M(:)), 1);

figure, imshow(Ibrain), hold on
visboundaries(M, 'Color','r');           % rojo = medicos
visboundaries(sangre_opt, 'Color','g');  % verde = nuestro
title(sprintf('Optimizado | Paciente %d, corte %d | Dice = %.2f', p, s, dice_opt));
drawnow

%% Mejores DELTA y R_OPEN para cada tipo 
tipos_nombres = {'Intraventricular','Intraparenchymal','Subarachnoid','Epidural','Subdural'};
N_por_tipo = 20;                 % cortes por tipo (inf = todos)
DELTAs = 0:5:80;
ROPENs = 0:5;
nT = numel(tipos_nombres);

% me creo un array con solo cortes con un unico tipo de hemorragia
sumaT = zeros(height(conH),1);
for i = 1:nT, sumaT = sumaT + conH.(tipos_nombres{i}); end
unico = conH(sumaT == 1, :);

%pongo vectores vacios para ir asignandoles los mejores valores
DELTA_opt = nan(nT,1);  
ROPEN_opt = nan(nT,1);  
Dice_max = nan(nT,1);  
nCortes = zeros(nT,1);

for i = 1:nT

    sub = unico(unico.(tipos_nombres{i}) == 1, :);% esto es para no agarrar toooodos los cortes 
    if height(sub) > N_por_tipo
        sub = sub(round(linspace(1, height(sub), N_por_tipo)), :);
    end
    C = cargar_set(sub);                         % cada corte se carga una sola vez
    nCortes(i) = numel(C);
    if isempty(C), continue; 
    end
    Dm = nan(numel(DELTAs), numel(ROPENs));      % Dice para cada combinacion
    for a = 1:numel(DELTAs)
        for b = 1:numel(ROPENs)
            Dm(a,b) = dice_set(C, DELTAs(a), ROPENs(b));
        end
    end % hasta aca lo que hago es ver la combinación de los Ropen y Delta para obtener el mejor dice
    [Dice_max(i), idx] = max(Dm(:)); %aca obtengo el mas alto
    [a, b] = ind2sub(size(Dm), idx);
    DELTA_opt(i) = DELTAs(a);  
    ROPEN_opt(i) = ROPENs(b);
    fprintf('%-18s DELTA=%d apertura=%d -> Dice=%.3f (%d cortes)\n', ...
        tipos_nombres{i}, DELTA_opt(i), ROPEN_opt(i), Dice_max(i), nCortes(i));
end

resumen = table(tipos_nombres', nCortes, DELTA_opt, ROPEN_opt, Dice_max, ...
    'VariableNames', {'Tipo','Cortes','DELTA','Apertura','Dice'});
disp(resumen)

%% Funciones auxiliares
function [I, M] = cargar(p, s)
    % lee la imagen en ventana de cerebro y la mascara de los medicos
    carpeta = fullfile('Patients_CT', sprintf('%03d', p), 'brain');
    I = imread(fullfile(carpeta, sprintf('%d.jpg', s)));
    M = imread(fullfile(carpeta, sprintf('%d_HGE_Seg.jpg', s)));
    if size(I,3) > 1, I = rgb2gray(I); end
    if size(M,3) > 1, M = M(:,:,1); end
    M = M > 128;
end

function [interior, cf, ref, craneo, disco] = segmentar_cerebro(I)
    craneo  = I > 240;
    cerrado = imclose(craneo, strel('disk', 15));  % tapa los huecos del anillo
    disco   = imfill(cerrado, 'holes');
    
    interior0 = disco & ~craneo;
    interior0 = bwareafilt(interior0, 1);
    
    if ~any(interior0(:))               % no se encontro el cerebro
        interior = false(size(I));
        cf = [];
        ref = NaN;
        return; 
    end
    
    interior = imerode(interior0, strel('disk', 2));
    if ~any(interior(:))
        cf = [];
        ref = NaN;
        return; 
    end
    
    cf  = medfilt2(I, [5 5]);
    ref = median(double(cf(interior)));   % gris tipico del cerebro de ESTE corte
end

function sangre = umbralizar(cf, ref, interior, DELTA, rOpen)
    sangre = false(size(cf));
    if isnan(ref) || isempty(cf), return; end
    sangre = (cf > ref + DELTA) & (cf < 245) & interior;   % mas claro que el cerebro
    if rOpen > 0, sangre = imopen(sangre, strel('disk', rOpen)); end
    sangre = imclose(sangre, strel('disk', 3));
    sangre = imfill(sangre, 'holes');
    sangre = bwareaopen(sangre, 30);
end

function C = cargar_set(sub)
    % carga una vez todos los cortes de la lista (los que fallan se saltean)
    C = {};
    for k = 1:height(sub)
        try
            [I, M] = cargar(sub.PatientNumber(k), sub.SliceNumber(k));
            [int, cf, ref] = segmentar_cerebro(I);
            if any(M(:)) && ~isnan(ref)
                C{end+1} = struct('M', M, 'int', int, 'cf', cf, 'ref', ref);
            end
        catch
        end
    end
end

function dm = dice_set(C, DELTA, rOpen)
    % Dice promedio sobre los cortes ya cargados
    d = nan(numel(C),1);
    for k = 1:numel(C)
        sk = umbralizar(C{k}.cf, C{k}.ref, C{k}.int, DELTA, rOpen);
        d(k) = 2*sum(sk(:) & C{k}.M(:)) / max(sum(sk(:)) + sum(C{k}.M(:)), 1);
    end
    dm = mean(d, 'omitnan');
end