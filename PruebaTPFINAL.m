%% ============================================================
%  SEGMENTACIÓN DE HEMORRAGIA EN VOLUMEN NIfTI (método relativo)
%  Idea: en cada corte, la sangre es lo que es más claro que el
%  cerebro "típico" de ESE corte (mediana) + un margen DELTA.
%  Después se limpia con morfología y se compara con la máscara
%  de los médicos usando el coeficiente de Dice.
%% ============================================================
% clear; clc; close all      % limpio variables, consola y figuras

%% ---------- 1) RUTAS DE LOS ARCHIVOS ----------
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta), carpeta = pwd; 
    end   % si se pegó en la consola, uso la carpeta actual

archCT    = fullfile(carpeta, 'CT.nii');                    % tomografía original
archBrain = fullfile(carpeta, 'CQ500-CT-4_brain.nii');      % volumen sin cráneo
archMask  = fullfile(carpeta, 'ICH_mask.nii');              % máscara de los médicos

% Verifico que existan los tres; si falta uno, aviso cuál
for a = {archCT, archBrain, archMask}
    if ~isfile(a{1}), error('No encuentro el archivo: %s', a{1}); end
end

%% ---------- 2) LECTURA DE LOS VOLÚMENES ----------
CT    = double(niftiread(archCT));      % pasa a double para poder filtrar y comparar
brain = double(niftiread(archBrain));   % lo uso solo para saber dónde está el cerebro
mask  = niftiread(archMask) > 0;        % máscara binaria (true = hemorragia)
nZ    = size(CT, 3);                    % cantidad de cortes axiales


%% ---------- 3) PARÁMETROS (lo único que hay que tocar) ----------
P.usarBrain = true;   % true: el interior sale del volumen sin cráneo; false: lo calculo desde el CT
P.huesoHU   = 100;    % umbral de hueso en HU (solo se usa si usarBrain = false)
P.erosion   = 3;      % radio de erosión: me alejo del borde del cráneo (evita falsos positivos)
P.delta     = 15;     % HU por encima del cerebro para considerar que es sangre
P.techo     = 100;    % tope en HU: por encima ya es hueso o calcificación, no sangre
P.rOpen     = 2;      % radio de la apertura morfológica (borra puntitos y hilos finos)
P.minArea   = 30;     % tamaño mínimo (en píxeles) de un objeto para conservarlo

%% ---------- 4) PRECÁLCULO POR CORTE ----------
interior = false(size(CT));   % máscara del interior del cerebro, corte a corte
cf       = zeros(size(CT));   % CT filtrado con mediana (reduce ruido)
ref      = nan(nZ, 1);        % gris típico del cerebro de cada corte
se       = strel('disk', P.erosion);   % elemento estructural para erosionar

for z = 1:nZ
    Cz = CT(:,:,z);   % corte z del CT

    % --- Defino el interior del cerebro ---
    if P.usarBrain
        % Todo lo que no es cero en el volumen "brain"; imfill tapa huecos internos
        int = imfill(brain(:,:,z) > 0, 'holes');
    else
        % Calculo el interior desde el CT: cráneo -> cierro anillo -> relleno -> resto el hueso
        craneo = Cz > P.huesoHU;
        disco  = imfill(imclose(craneo, strel('disk', 15)), 'holes');
        int    = disco & ~craneo;
    end
    if any(int(:)), int = bwareafilt(int, 1); end  % me quedo con el objeto más grande (el cerebro)
    int = imerode(int, se);                        % achico un poco para alejarme del hueso
    interior(:,:,z) = int;

    % --- Filtro de mediana 5x5: suaviza el ruido sin borrar bordes ---
    cf(:,:,z) = medfilt2(Cz, [5 5]);

    % --- Gris típico del cerebro en este corte (mediana dentro del interior) ---
    if any(int(:))
        tmp = cf(:,:,z);
        ref(z) = median(tmp(int));
    end
end

%% ---------- 5) SEGMENTACIÓN DE TODO EL VOLUMEN ----------
sangre = false(size(CT));   % resultado: máscara de sangre detectada
for z = 1:nZ
    sangre(:,:,z) = segmentar(cf(:,:,z), ref(z), interior(:,:,z), P);
end

% ---------- 6) EVALUACIÓN CON DICE ----------
conMask = squeeze(any(any(mask, 1), 2));   % cortes donde los médicos marcaron hemorragia

dice3D = 2*sum(sangre(:) & mask(:)) / max(sum(sangre(:)) + sum(mask(:)), 1);

% Dice por corte (solo cortes con hemorragia)
diceZ = zeros(nZ, 1);
for z = find(conMask)'
    S = sangre(:,:,z);  M = mask(:,:,z);
    diceZ(z) = 2*sum(S(:) & M(:)) / max(sum(S(:)) + sum(M(:)), 1);
end

fprintf('\nDice 3D (volumen completo):          %.3f\n', dice3D);
fprintf('Dice promedio (cortes con hemorr.):  %.3f\n', mean(diceZ(conMask)));
fprintf('Voxeles: nuestro=%d, medicos=%d, acierto=%d\n', ...
    nnz(sangre), nnz(mask), nnz(sangre & mask));

figure
plot(find(conMask), diceZ(conMask), 'o-');
xlabel('Corte'); 
ylabel('Dice'); 
title('Dice por corte'); 
grid on

%% ---------- 7) BARRIDO RÁPIDO DE DELTA ----------
barrer = true;   % poner en false para saltear esta sección
if barrer
    deltas = 5:5:40;
    dd = zeros(size(deltas));
    for i = 1:numel(deltas)
        Q = P;  Q.delta = deltas(i);
        S = false(size(CT));
        for z = find(conMask)'
            S(:,:,z) = segmentar(cf(:,:,z), ref(z), interior(:,:,z), Q);
        end
        A = S(:,:,conMask);  B = mask(:,:,conMask);
        dd(i) = 2*sum(A(:) & B(:)) / max(sum(A(:)) + sum(B(:)), 1);
        fprintf('delta = %2d -> Dice = %.3f\n', deltas(i), dd(i));
    end
    figure, plot(deltas, dd, 'o-'); grid on
    xlabel('DELTA (HU)'); ylabel('Dice'); title('Dice según DELTA');
    [~, k] = max(dd);
    fprintf('Mejor DELTA: %d (Dice %.3f). Cargalo en P.delta y volvé a correr.\n', deltas(k), dd(k));
end

%% ---------- 8) VISOR: rojo = médicos, verde = nuestro ----------
CTv = mat2gray(CT, [0 80]);        % ventana de cerebro (0 a 80 HU) escalada a 0..1 para mostrar
z0  = find(conMask, 1);            % arranco en el primer corte con hemorragia
if isempty(z0), z0 = round(nZ/2); end

f  = figure('Name','Segmentación vs máscara','Position',[100 100 700 650]);
ax = axes(f);
% repmat a 3 canales para que visboundaries dibuje bien sobre la imagen
h  = imshow(repmat(rot90(CTv(:,:,z0)), [1 1 3]), 'Parent', ax); hold(ax, 'on');

% Guardo en la figura todo lo que necesita el callback del slider
d.CTv = CTv; d.mask = mask; d.sangre = sangre; d.diceZ = diceZ;
d.ax = ax; d.h = h; d.nZ = nZ;
f.UserData = d;

% Slider para recorrer los cortes
uicontrol('Style','slider','Min',1,'Max',nZ,'Value',z0, ...
    'SliderStep',[1/(nZ-1) 5/(nZ-1)], 'Units','normalized', ...
    'Position',[0.2 0.02 0.6 0.04], ...
    'Callback', @(s,~) actualizar(f, round(s.Value)));
actualizar(f, z0);   % dibujo el corte inicial
reconstruir3D_suave('craneo')       % craneo, desde CT.nii
reconstruir3D_suave('cerebro')      % cerebro, desde el volumen sin craneo
reconstruir3D_suave('hemorragia')   % lo que marco el medico (ICH_mask.nii)
reconstruir3D_suave('todo')         % las tres juntas
%% ============ FUNCIONES LOCALES ============

function sangre = segmentar(cf, ref, interior, P)
    % Segmenta la sangre de UN corte con el método relativo.
    sangre = false(size(cf));                          % por defecto, nada
    if isnan(ref) || ~any(interior(:)), return; end    % corte sin cerebro válido: salgo

    % Más claro que el cerebro del corte (+delta), menos que el techo, y solo dentro del cerebro
    sangre = (cf > ref + P.delta) & (cf < P.techo) & interior;

    % Apertura: erosiona y dilata, borra puntitos y hilos finos
    if P.rOpen > 0, sangre = imopen(sangre, strel('disk', P.rOpen)); end
    % Cierre: une zonas de sangre cercanas y tapa pequeños cortes
    sangre = imclose(sangre, strel('disk', 3));
    % Relleno: tapa agujeros dentro de cada región de sangre
    sangre = imfill(sangre, 'holes');
    % Descarto objetos más chicos que minArea
    sangre = bwareaopen(sangre, P.minArea);
end

function actualizar(f, z)
    % Callback del slider: redibuja el corte z con los contornos.
    d = f.UserData;
    d.h.CData = repmat(rot90(d.CTv(:,:,z)), [1 1 3]);   % imagen del corte
    delete(findobj(d.ax, 'Type', 'line'));              % borro contornos del corte anterior
    m = rot90(d.mask(:,:,z));                           % máscara de médicos (misma rotación)
    s = rot90(d.sangre(:,:,z));                         % nuestra segmentación
    if any(m(:)), visboundaries(d.ax, m, 'Color', 'r'); end   % rojo = médicos
    if any(s(:)), visboundaries(d.ax, s, 'Color', 'g'); end   % verde = nuestro
    title(d.ax, sprintf('Corte %d/%d | Dice = %.3f', z, d.nZ, d.diceZ(z)));
end
