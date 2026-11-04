%% TP Final PDI segmentacion de hemorragia en un paciente
% weno solo en un paciente donde anda bien:
%  1. los pasos del metodo en un corte (CT -> cerebro -> candidatos -> final)
%  2. varios cortes con la hemorragia (rojo = medicos, verde = nuestro)
%  3. un visor con slider para recorrer todos los cortes
%  4. la reconstruccion 3D (medicos vs nuestro)
%
% metodo final (lo elegimos con PruebaTPFINAL.m):
%  - saco el craneo y me quedo con la region cerebral
%  - umbral relativo: sangre = mas clara que el cerebro de ese corte + DELTA,
%    con histeresis (umbral fuerte para las semillas, debil para crecer un poco)
%  - armo los objetos 3D y me quedo con los que:
%       son gruesos (radio >= RADIO_MIN)             -> filtro de forma
%       son brillantes en promedio (contraste medio)  -> filtro de densidad
%       no tienen picos de calcio/hueso (pico < TECHO)
%
% chequeo: el Dice tiene que dar igual que en PruebaTPFINAL (ej: CT-225 -> 0.927)

clear; clc; close all

%% 1) Que paciente muestro
DATASET = fullfile(pwd, 'Seg-CQ500');   % carpeta del dataset 
CASO    = 'CQ500-CT-225';               % cambiar por cualquiera de la lista de abajo

% casos donde el metodo anda bien (Dice 3D con el metodo final):
%   grandes : CQ500-CT-75  (94 mL, 0.79) | CQ500-CT-488 (69 mL, 0.80) | CQ500-CT-225 (54 mL, 0.93)
%   medianas: CQ500-CT-226 (21 mL, 0.84) | CQ500-CT-249 (15 mL, 0.89)
%   chicas  : CQ500-CT-22  (6 mL, 0.90)  | CQ500-CT-10  (3 mL, 0.84)  | CQ500-CT-196 (0.5 mL, 0.89)

VER_3D           = true;    % reconstruccion 3D: medicos vs nuestro
VER_3D_COMPANERA = false;   % true = ademas la reconstruccion de reconstruir3D_suave (craneo + cerebro + hemorragia de los medicos)

%% 2) Parametros del metodo final (los mismos que eligio PruebaTPFINAL.m)
P.HUESO       = 100;   % HU: arriba de esto es hueso
P.EROSION     = 8;     % radio (px) con el que me alejo del borde del craneo
P.MIN_CEREBRO = 5000;  % px: si un corte tiene menos cerebro que esto, lo ignoro (base y vertice)
P.HISTERESIS  = 8;     % HU: el umbral debil es el umbral fuerte menos esto
P.TECHO       = 90;    % HU maximo de la sangre (mas que esto es calcio o hueso)
P.R_OPEN      = 2;     % radio de la apertura (la hago sobre las semillas)
P.AREA_MIN    = 100;   % px: semillas mas chicas que esto se descartan
P.CRECER      = 5;     % voxeles que puede crecer cada semilla
P.VOL_MIN     = 1000;  % voxeles: objetos 3D mas chicos que esto se descartan
P.DELTA       = 15;    % HU sobre el cerebro (umbral fuerte)
P.RADIO_MIN   = 4;     % voxeles: radio minimo del objeto (filtro de forma)
P.CONT_MIN    = 16;    % HU: contraste medio minimo del objeto (filtro de densidad)

%% 3) Busco y leo el paciente
mascaras = dir(fullfile(DATASET, '**', '*mask*.nii*'));   % una mascara por paciente
if isempty(mascaras), error('No encontre mascaras dentro de %s.', DATASET); end
[~, carpetas] = cellfun(@fileparts, {mascaras.folder}, 'UniformOutput', false);   % nombre de cada paciente (= su carpeta)
k = find(strcmp(carpetas, CASO), 1);
if isempty(k), error('No encontre el paciente %s en %s.', CASO, DATASET); end

fprintf('Paciente: %s\n', CASO);
[CT, mask, info, archCT, archMask] = leerPaciente(mascaras(k));
vox    = info.PixelDimensions(1:3);                % tamano del voxel (mm)
volVox = prod(vox) / 1000;                         % mL por voxel
fprintf('Volumen %d x %d x %d | voxel %.2f x %.2f x %.2f mm\n', size(CT), vox);

%% 4) Segmento con el metodo final
tic
interior  = calcularInterior(CT, P);               % region cerebral (sin craneo)
[cf, ref] = filtrarYReferencia(CT, interior);      % CT con mediana y gris tipico del cerebro de cada corte
cand      = candidatos(cf, ref, interior, P.DELTA, P);   % umbral + histeresis + limpieza
sangre    = filtrarObjetos(cand, cf, ref, P);      % filtros de forma, densidad y pico
fprintf('Segmentacion lista (%.0f s)\n', toc);

CTv = im2uint8(mat2gray(CT, [0 80]));              % ventana de cerebro (0 a 80 HU) para mostrar; uint8 gasta menos memoria
clear CT                                           % ya no lo necesito

%% 5) Que tan bien anduvo
TP   = nnz(sangre & mask);                         % voxeles acertados
dice = 2 * TP / max(nnz(sangre) + nnz(mask), 1);   % 1 = perfecto, 0 = nada en comun
sens = TP / max(nnz(mask), 1);                     % de la sangre real, cuanta encontre
prec = TP / max(nnz(sangre), 1);                   % de lo que marque, cuanto era sangre
fprintf('\nDice 3D: %.3f | Sensibilidad: %.3f | Precision: %.3f\n', dice, sens, prec);
fprintf('Volumen -> nuestro: %.1f mL | medicos: %.1f mL\n', nnz(sangre) * volVox, nnz(mask) * volVox);

area      = squeeze(sum(sum(mask, 1), 2));         % pixeles de hemorragia en cada corte
[~, zMax] = max(area);                             % el corte con mas hemorragia

%% 6) Figura 1: los pasos del metodo en el corte con mas hemorragia
c = rot90(CTv(:,:,zMax));                          % rot90 es solo para que se vea derecho
figure('Name', 'Pasos del metodo')
subplot(2, 2, 1); imshow(c); title('1. CT original (ventana 0-80 HU)')
subplot(2, 2, 2); imshow(c); hold on
if any(interior(:,:,zMax), 'all'), visboundaries(rot90(interior(:,:,zMax)), 'Color', 'c'); end
title('2. Region cerebral (saco el craneo)')
subplot(2, 2, 3); imshow(labeloverlay(c, rot90(cand(:,:,zMax)), 'Colormap', [1 1 0], 'Transparency', 0.4))
title('3. Candidatos (umbral + histeresis)')
subplot(2, 2, 4); imshow(c); hold on
if any(mask(:,:,zMax), 'all'),   visboundaries(rot90(mask(:,:,zMax)),   'Color', 'r'); end
if any(sangre(:,:,zMax), 'all'), visboundaries(rot90(sangre(:,:,zMax)), 'Color', 'g'); end
title(sprintf('4. Final (verde) vs medicos (rojo) | Dice corte %.2f', diceCorte(sangre, mask, zMax)))
sgtitle(sprintf('%s - corte %d', CASO, zMax))

%% 7) Figura 2: seis cortes repartidos a lo largo de la hemorragia
zs      = find(area > 0);                          % cortes donde hay hemorragia
mostrar = unique(zs(round(linspace(1, numel(zs), min(6, numel(zs))))));   % 6 cortes repartidos
figure('Name', 'Cortes con hemorragia')
for i = 1:numel(mostrar)
    z = mostrar(i);
    subplot(2, 3, i)
    imshow(rot90(CTv(:,:,z))); hold on
    if any(mask(:,:,z), 'all'),   visboundaries(rot90(mask(:,:,z)),   'Color', 'r'); end
    if any(sangre(:,:,z), 'all'), visboundaries(rot90(sangre(:,:,z)), 'Color', 'g'); end
    title(sprintf('Corte %d | Dice %.2f', z, diceCorte(sangre, mask, z)))
end
sgtitle(sprintf('%s | Dice 3D = %.2f (rojo = medicos, verde = nuestro)', CASO, dice))

%% 8) Figura 3: visor con slider para recorrer todos los cortes
fig = figure('Name', 'Visor de cortes');
ax  = axes(fig, 'Position', [0.05 0.12 0.9 0.82]);
nZ  = size(CTv, 3);
sl  = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.15 0.03 0.7 0.04], ...
    'Min', 1, 'Max', nZ, 'Value', zMax, 'SliderStep', [1 10] / max(nZ - 1, 1));
mover = @(src, ~) dibujarCorte(ax, round(src.Value), CTv, mask, sangre);   % que hago cuando muevo el slider
sl.Callback = mover;
addlistener(sl, 'ContinuousValueChange', mover);   % asi se actualiza mientras arrastro
dibujarCorte(ax, zMax, CTv, mask, sangre);         % arranco en el corte con mas hemorragia

%% 9) Figura 4: reconstruccion 3D
if VER_3D
    ver3D(interior, mask, sangre, vox, volVox, sprintf('%s | Dice 3D = %.2f', CASO, dice));
end
if VER_3D_COMPANERA
    if exist('reconstruir3D_suave', 'file') == 2
        reconstruir3D_suave('todo', archCT, archMask, interior)
    else
        warning('No encontre reconstruir3D_suave.m en esta carpeta.');
    end
end

%% ============ FUNCIONES ============

function [CT, mask, info, archCT, archMask] = leerPaciente(m)
    % leo el CT (ya en HU) y la mascara de los medicos de un paciente
    carpeta  = m.folder;
    archMask = fullfile(carpeta, m.name);

    todos = dir(fullfile(carpeta, '*.nii*'));      % todos los .nii de la carpeta del paciente
    nom   = {todos.name};
    esCT  = ~contains(nom, 'mask', 'IgnoreCase', true) & ...
            ~contains(nom, 'brain', 'IgnoreCase', true);   % el CT es el que no es mascara ni "brain"
    if ~any(esCT), error('No encontre el CT'); end
    archCT = fullfile(carpeta, nom{find(esCT, 1)});

    info = niftiinfo(archCT);
    CT   = double(niftiread(archCT));              % paso a double para poder filtrar y hacer cuentas
    pend = info.MultiplicativeScaling;  if pend == 0, pend = 1; end   % 0 significa "sin escala"
    CT   = CT * pend + info.AdditiveOffset;        % a unidades HU (niftiread no lo hace solo)

    mask = niftiread(archMask) > 0;                % mascara binaria: true = hemorragia
    if ~isequal(size(CT), size(mask)), error('CT y mascara con distinto tamano'); end
end

function sangre = filtrarObjetos(cand, cf, ref, P)
    % para cada objeto 3D candidato decido si queda (los 3 filtros del metodo final)
    % es lo mismo que hace PruebaTPFINAL con las tablas de objetos, pero aca armo la mascara
    sangre = false(size(cand));
    if ~any(cand(:)), return; end
    D      = bwdist(~cand);                        % distancia al borde = radio de la esfera que entra
    cc     = bwconncomp(cand, 26);                 % separo los objetos 3D
    refPac = median(ref, 'omitnan');               % HU tipico del cerebro del paciente
    n      = zeros(1, 4);                          % cuantos: [saco por radio, por contraste, por pico, quedan]
    fprintf('\nObjetos candidatos: %d (muestro los que pasan el filtro de forma)\n', cc.NumObjects);
    for k = 1:cc.NumObjects
        idx       = cc.PixelIdxList{k};
        [~, ~, z] = ind2sub(size(cand), idx);      % en que corte esta cada voxel
        contraste = double(cf(idx)) - ref(z);      % HU por encima del cerebro de su corte
        radio     = max(D(idx));
        cMed      = mean(contraste);
        pico      = max(contraste) + refPac;       % pico en HU (lo calculo igual que en PruebaTPFINAL)

        if radio < P.RADIO_MIN
            n(1) = n(1) + 1;  continue             % muy finito: lo saco (son muchos, no los muestro)
        elseif cMed < P.CONT_MIN
            n(2) = n(2) + 1;  motivo = 'lo saco: poco contraste';
        elseif pico >= P.TECHO
            n(3) = n(3) + 1;  motivo = 'lo saco: pico de calcio/hueso';
        else
            n(4) = n(4) + 1;  motivo = 'QUEDA';
            sangre(idx) = true;
        end
        fprintf('   radio %5.1f | %7d vox | contraste medio %5.1f HU | pico %5.0f HU -> %s\n', ...
            radio, numel(idx), cMed, pico, motivo);
    end
    fprintf('Saco %d por radio, %d por contraste y %d por pico. Quedan %d objetos.\n', n);
end

function d = diceCorte(sangre, mask, z)
    % dice de un solo corte
    s = sangre(:,:,z);  m = mask(:,:,z);
    d = 2 * nnz(s & m) / max(nnz(s) + nnz(m), 1);
end

function dibujarCorte(ax, z, CTv, mask, sangre)
    % dibujo el corte z en el visor (lo llama el slider)
    cla(ax);
    imshow(rot90(CTv(:,:,z)), 'Parent', ax); hold(ax, 'on');
    if any(mask(:,:,z), 'all'),   visboundaries(ax, rot90(mask(:,:,z)),   'Color', 'r'); end
    if any(sangre(:,:,z), 'all'), visboundaries(ax, rot90(sangre(:,:,z)), 'Color', 'g'); end
    hold(ax, 'off');
    title(ax, sprintf('Corte %d de %d   (rojo = medicos, verde = nuestro)', z, size(CTv, 3)));
end

function ver3D(cerebro, mask, sangre, vox, volVox, titulo)
    % reconstruccion 3D: a la izquierda los medicos, a la derecha nuestro.
    % el cerebro va transparente para ubicarse. Las dos vistas giran juntas.
    figure('Name', 'Reconstruccion 3D', 'Color', 'w');
    hem     = {mask, sangre};
    nombres = {'Medicos', 'Nuestro'};
    colores = {[0.85 0.10 0.10], [0.10 0.70 0.20]};
    ax      = gobjects(1, 2);
    for i = 1:2
        ax(i) = subplot(1, 2, i);
        superficie3D(cerebro, 2, vox, [0.75 0.75 0.75], 0.08, 0.1);   % cerebro: 1 de cada 2 pixeles y pocas caras (es enorme)
        superficie3D(hem{i},  1, vox, colores{i}, 1, 0.5);            % hemorragia: resolucion completa
        daspect([1 1 1]); view(3); axis tight vis3d off                % ya esta en mm, asi que 1:1:1
        camlight headlight; lighting gouraud
        title(sprintf('%s (%.1f mL)', nombres{i}, nnz(hem{i}) * volVox))
    end
    setappdata(gcf, 'enlace', linkprop(ax, {'CameraPosition', 'CameraUpVector', 'CameraTarget', 'CameraViewAngle'}));   % giran juntas
    rotate3d on                                    % para girarlo con el mouse
    sgtitle(titulo)
end

function superficie3D(BW, paso, vox, color, alfa, fraccion)
    % dibujo la superficie de un volumen binario (isosurface), suavizada y en mm
    if ~any(BW(:)), return; end
    % recorto la caja donde esta el objeto (asi no trabajo con todo el volumen)
    f = find(any(any(BW, 2), 3));  c = find(any(any(BW, 1), 3));  s = find(any(any(BW, 1), 2));
    f1 = max(f(1) - 2, 1);  f2 = min(f(end) + 2, size(BW, 1));
    c1 = max(c(1) - 2, 1);  c2 = min(c(end) + 2, size(BW, 2));
    s1 = max(s(1) - 2, 1);  s2 = min(s(end) + 2, size(BW, 3));
    V  = smooth3(double(BW(f1:paso:f2, c1:paso:c2, s1:s2)), 'box', 3);   % suavizo para que no se vean escalones
    fv = isosurface(V, 0.5);                       % superficie donde pasa de 0 a 1
    if isempty(fv.vertices), return; end
    fv = reducepatch(fv, fraccion);                % menos triangulos = gira mas rapido
    % paso de indices de la caja a mm (x = columna, y = fila, z = corte)
    v = fv.vertices;
    fv.vertices = [(c1 - 1 + (v(:,1) - 1) * paso) * vox(2), ...
                   (f1 - 1 + (v(:,2) - 1) * paso) * vox(1), ...
                   (s1 - 1 + (v(:,3) - 1))        * vox(3)];
    patch(fv, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alfa);
end

function interior = calcularInterior(CT, P)
    % saco el craneo: hueso -> cierro el anillo -> relleno -> resto el hueso -> erosiono
    interior = false(size(CT));
    for z = 1:size(CT, 3)
        corte  = CT(:,:,z);
        craneo = corte > P.HUESO;                                   % 1. todo lo que es hueso
        disco  = imfill(imclose(craneo, strel('disk', 15)), 'holes');   % 2. cierro el anillo del craneo y lo relleno
        int    = disco & ~craneo;                                   % 3. lo de adentro, sin el hueso
        int    = mayorObjeto(int);                                  % 4. me quedo con el objeto mas grande (el cerebro)
        int    = imerode(int, strel('disk', P.EROSION));            % 5. achico para alejarme del hueso
        if nnz(int) < P.MIN_CEREBRO, int(:) = false; end            % corte de base o vertice: lo anulo
        interior(:,:,z) = int;
    end
end

function [cf, ref] = filtrarYReferencia(CT, interior)
    % filtro de mediana y gris tipico del cerebro de cada corte
    nZ  = size(CT, 3);
    cf  = zeros(size(CT), 'single');               % single para gastar menos memoria
    ref = nan(nZ, 1);
    for z = 1:nZ
        cf(:,:,z) = medfilt2(CT(:,:,z), [5 5]);    % mediana 5x5: saca ruido sin borrar bordes
        int = interior(:,:,z);
        if any(int(:))
            tmp    = cf(:,:,z);
            ref(z) = median(tmp(int));             % mediana del cerebro de este corte
        end
    end
end

function sangre = candidatos(cf, ref, interior, DELTA, P)
    % umbral con histeresis limitada + limpieza 3D (todavia sin el filtro de forma)
    fuerte = false(size(cf));                      % semillas: seguro es sangre
    debil  = false(size(cf));                      % candidatos: puede ser sangre
    for z = 1:size(cf, 3)
        if isnan(ref(z)), continue; end            % corte sin cerebro: lo salteo
        c = cf(:,:,z);
        f = (c > ref(z) + DELTA)                 & (c < P.TECHO) & interior(:,:,z);   % umbral fuerte
        d = (c > ref(z) + DELTA - P.HISTERESIS)  & (c < P.TECHO) & interior(:,:,z);   % umbral debil
        if P.R_OPEN > 0, f = imopen(f, strel('disk', P.R_OPEN)); end   % apertura: borra puntitos e hilos finos
        f = bwareaopen(f, P.AREA_MIN);             % saco semillas muy chicas
        fuerte(:,:,z) = f;
        debil(:,:,z)  = d | f;                     % el debil incluye al fuerte
    end

    sangre  = fuerte;
    vecinos = true(3, 3, 3);                       % los 26 vecinos en 3D
    for k = 1:P.CRECER                             % cada vuelta crece 1 voxel, pero solo sobre los candidatos
        sangre = imdilate(sangre, vecinos) & debil;
    end

    for z = 1:size(sangre, 3)
        s = sangre(:,:,z);
        if ~any(s(:)), continue; end
        sangre(:,:,z) = imfill(imclose(s, strel('disk', 3)), 'holes');   % cierre une zonas cercanas, imfill tapa agujeros
    end
    sangre = bwareaopen(sangre, P.VOL_MIN, 26);    % limpieza 3D: saco objetos chicos
end

function BW = mayorObjeto(BW)
    % dejo solo el objeto mas grande
    cc = bwconncomp(BW);
    if cc.NumObjects == 0, return; end
    [~, k] = max(cellfun(@numel, cc.PixelIdxList));
    BW = false(size(BW));
    BW(cc.PixelIdxList{k}) = true;
end