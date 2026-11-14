%% TP Final PDI - Segmentacion de hemorragia en CT (Seg-CQ500, NIfTI)
% En cada corte, la sangre es lo que es mas claro que el cerebro tipico de ESE
% corte (mediana) + DELTA. Pasos:
%   1. Cerebro: adentro del craneo, a mas de MARGEN mm del hueso (3D), objeto 3D mas grande
%   2. Sangre: mas claro que el gris del corte + DELTA y menos que TECHO
%   3. Limpieza 2D (apertura, cierre, relleno) y 3D (espesor minimo GROSOR_MIN mm)
%   4. Histeresis: cada objeto necesita algun punto mas denso que el cerebro + NUCLEO
% Se compara con la mascara de los medicos (Dice), se guardan las dos mascaras por
% separado y se abre un visor. Solo para evaluar: naranja = falso positivo (objeto
% nuestro que no toca ninguno de los medicos), azul = falso negativo (hemorragia de
% los medicos que no toca ninguno de los nuestros).

clear; clc; close all

%% 1) Parametros
DATASET    = fullfile(pwd, 'Seg-CQ500');    % carpeta del dataset (al lado de este script)
RESULTADOS = fullfile(pwd, 'resultados');   % aca se guardan las mascaras separadas
CTS        = [55  249 90 227 95];           % numero de CT del CQ500 (carpeta CQ500CT227, ...)
HUESO      = 100;    % HU: por encima es hueso
MARGEN     = 6;      % mm: distancia minima al hueso, en 3D
DELTA      = 15;     % HU por encima del cerebro para considerar sangre
TECHO      = 100;    % HU: por encima ya no es sangre
R_OPEN     = 2;      % pixeles: radio de la apertura en cada corte
GROSOR_MIN = 4;      % mm: espesor minimo en 3D (0 = sin filtro 3D)
NUCLEO     = 25;     % HU sobre el cerebro que necesita algun punto de cada objeto
guardar    = true;   % guarda las dos mascaras (.nii.gz) y la tabla de Dice
barrer     = false;  % prueba varios DELTA en todas las CT (tarda mas)
deltas     = 5:5:40;
%El dice que mejor queda es el de margen 10 y grosor 4 
%% 2) Buscar las mascaras y el numero de CT de cada una
mascaras = dir(fullfile(DATASET, '**', '*mask*.nii*'));
if isempty(mascaras), error('No encontre mascaras dentro de %s.', DATASET); end
nums = cellfun(@(t) str2double(char(t)), regexpi(strcat({mascaras.folder}, filesep, ...
       {mascaras.name}), 'CQ500[-_ ]*CT[-_ ]*(\d+)', 'tokens', 'once'));
if guardar && ~exist(RESULTADOS, 'dir'), mkdir(RESULTADOS); end

%% 3) Procesar cada CT
res = {};  tabla = [];  dd = [];
for ct = CTS
    k = find(nums == ct, 1);
    if isempty(k), warning('No encontre la CT %d.', ct); continue; end
    t0 = tic;  carpeta = mascaras(k).folder;  archMask = fullfile(carpeta, mascaras(k).name);

    % Leer el CT en HU (el .nii de la carpeta que no es mascara ni "brain") y la mascara
    f = dir(fullfile(carpeta, '*.nii*'));
    f = f(~contains({f.name}, {'mask', 'brain'}, 'IgnoreCase', true));
    info = niftiinfo(fullfile(carpeta, f(1).name));
    m = info.MultiplicativeScaling;  if m == 0, m = 1; end       % 0 = "sin escala"
    CT   = double(niftiread(info)) * m + info.AdditiveOffset;
    mask = niftiread(archMask) > 0;                               % medicos
    esp  = reshape(double(info.PixelDimensions(1:3)), 1, 3);      % voxel en mm

    % Nuestra mascara
    [interior, cf, ref] = prepararVolumen(CT, esp, HUESO, MARGEN);
    sangre = segmentar(cf, ref, interior, esp, DELTA, TECHO, R_OPEN, GROSOR_MIN, NUCLEO);

    % Evaluar: Dice 3D, Dice medio en los cortes con hemorragia, falsos positivos y negativos
    diceZ = arrayfun(@(z) calcDice(sangre(:,:,z), mask(:,:,z)), find(squeeze(any(any(mask, 1), 2)))');
    fp = sangre & ~mask;
    fn = mask & ~sangre;             % de los medicos que no encontramos
    mL = prod(esp) / 1000;                                        % mL de un voxel
    tabla(end+1, :) = [ct, calcDice(sangre, mask), mean(diceZ), nnz(mask) * mL, nnz(sangre) * mL, ...
                       nnz(fp) * mL, nnz(fn) * mL];
    fprintf(['CT %d: Dice 3D %.3f | por corte %.3f | medicos %.1f mL | nuestra %.1f mL | ' ...
             'falso positivo %.1f mL | falso negativo %.1f mL (%.0f s)\n'], tabla(end, :), toc(t0));

    % Guardar las dos mascaras por separado, con el encabezado de la de los medicos
    if guardar
        try
            im = niftiinfo(archMask);
            im.Datatype = 'uint8';  im.BitsPerPixel = 8;  im.ImageSize = size(mask);
            im.PixelDimensions = im.PixelDimensions(1:3);
            niftiwrite(uint8(mask),   fullfile(RESULTADOS, sprintf('CT%d_medicos', ct)), im, 'Compressed', true);
            niftiwrite(uint8(sangre), fullfile(RESULTADOS, sprintf('CT%d_nuestra', ct)), im, 'Compressed', true);
        catch err, warning('No pude guardar las mascaras de la CT %d: %s', ct, err.message); end
    end
    if barrer       % Dice 3D para cada DELTA
        dd(end+1, :) = arrayfun(@(D) calcDice(segmentar(cf, ref, interior, esp, D, TECHO, ...
                                R_OPEN, GROSOR_MIN, NUCLEO), mask), deltas);
    end

    res{end+1} = prepararVisor(tabla(end, :), CT, mask, sangre, fp, fn, interior, esp, info, HUESO);
    clear CT cf interior sangre mask fp fn                        % libero memoria
end
if isempty(res), error('No se pudo procesar ninguna CT.'); end

%% 4) Resumen (y barrido de DELTA si barrer = true)
resumen = array2table(tabla, 'VariableNames', ...
    {'CT', 'Dice3D', 'DiceCortes', 'VolMedicos_mL', 'VolNuestra_mL', 'FalsoPositivo_mL', 'FalsoNegativo_mL'});
disp(resumen)
fprintf('Promedio: Dice 3D %.3f | Dice por corte %.3f\n', mean(tabla(:, 2)), mean(tabla(:, 3), 'omitnan'));
if guardar, writetable(resumen, fullfile(RESULTADOS, 'resumen_dice.csv')); end
if barrer
    figure('Name', 'Dice segun DELTA');  grid on;  hold on
    plot(deltas, dd', 'o-', deltas, mean(dd, 1), 'k-', 'LineWidth', 1.5);
    xlabel('DELTA (HU)');  ylabel('Dice 3D');
    legend([arrayfun(@(c) sprintf('CT %d', c), tabla(:, 1), 'UniformOutput', false); {'Promedio'}]);
    [mejor, i] = max(mean(dd, 1));
    fprintf('Mejor DELTA en promedio: %d (Dice %.3f)\n', deltas(i), mejor);
end

%% 5) Visor: una CT por vez (3D y cortes), Medicos | Nuestra
visor(res);

%% ============ FUNCIONES ============

function [interior, cf, ref] = prepararVolumen(CT, esp, HUESO, MARGEN)
    % Cerebro (interior), CT filtrado con mediana (cf) y gris tipico de cada corte (ref)
    nZ = size(CT, 3);
    hueso = mayor3D(imreconstruct(CT > 300, CT > HUESO));    % > HUESO pegado a hueso denso (un coagulo no cuenta)
    cerca = imdilate(hueso, strel('cuboid', 2 * round(MARGEN ./ esp) + 1));   % a menos de MARGEN mm, en 3D
    interior = false(size(CT));  cf = CT;
    for z = 1:nZ
        disco = imfill(imclose(hueso(:,:,z), strel('disk', 15)), 'holes');   % craneo cerrado y relleno
        interior(:,:,z) = disco & ~cerca(:,:,z);                             % adentro y lejos del hueso
        cf(:,:,z) = medfilt2(CT(:,:,z), [5 5]);                              % saca ruido sin borrar bordes
    end
    interior = mayor3D(interior);                   % objeto 3D mas grande: saca los huecos de la base
    ref = nan(nZ, 1);
    for z = 1:nZ, c = cf(:,:,z);  ref(z) = median(c(interior(:,:,z))); end   % NaN si no hay cerebro
end

function sangre = segmentar(cf, ref, interior, esp, DELTA, TECHO, R_OPEN, GROSOR_MIN, NUCLEO)
    R = reshape(ref, 1, 1, []);                                        % gris del cerebro de cada corte
    s = (cf > R + DELTA) & (cf < TECHO) & interior;                    % 1. umbral relativo
    s = imclose(imopen(s, strel('disk', R_OPEN)), strel('disk', 3));   % 2. limpieza en cada corte
    for z = 1:size(s, 3), s(:,:,z) = imfill(s(:,:,z), 'holes'); end
    s = s & zonaGruesa(s, esp, GROSOR_MIN);                            % 3. espesor minimo en 3D
    sangre = imreconstruct(s & (cf > R + NUCLEO), s);                  % 4. histeresis: nucleo denso
end

function zona = zonaGruesa(S, esp, G)
    % Zona alrededor de donde entra una esfera de G mm de diametro (en una grilla de
    % 1 mm, para que sea redonda). Venas, hoz y rayas son finas: ahi no entra.
    if G <= 0 || ~any(S(:)), zona = true(size(S)); return; end
    S1     = imresize3(uint8(S) * 255, round(size(S) .* esp), 'linear') > 127;
    nucleo = bwdist(~S1) > G / 2;              % centros donde entra la esfera
    zona1  = bwdist(nucleo) <= G;              % alrededor del nucleo (recupera el borde)
    zona   = imresize3(uint8(zona1), size(S), 'nearest') > 0;
end

function BW = mayor3D(BW)
    % Se queda con el objeto 3D mas grande
    cc = bwconncomp(BW);
    if cc.NumObjects == 0, return; end
    [~, k] = max(cellfun(@numel, cc.PixelIdxList));
    BW = false(size(BW));  BW(cc.PixelIdxList{k}) = true;
end

function d = calcDice(A, B)
    d = 2 * nnz(A & B) / max(nnz(A) + nnz(B), 1);   % 1 = coincidencia perfecta, 0 = nada en comun
end

function r = prepararVisor(datos, CT, mask, sangre, fp, fn, interior, esp, info, HUESO)
    % Superficies 3D (en mm) y cortes 2D que usa el visor
    r.datos = datos;                     % [CT Dice3D DiceCortes VolMed VolNuestra FalsoPos FalsoNeg]
    p = max(1, round(size(CT, 1) / 128));  pz = max(1, round(p * esp(1) / esp(3)));
    baja = @(V) V(1:p:end, 1:p:end, 1:pz:end);       % craneo y cerebro: alcanza con menos resolucion
    r.fv = [superficie(mask & ~fn, esp), superficie(fn, esp), ...         % 1 medicos, 2 falso negativo
            superficie(sangre & ~fp, esp), superficie(fp, esp), ...       % 3 nuestra, 4 falso positivo
            superficie(imfill(baja(interior), 'holes'), esp .* [p p pz]), ...   % 5 cerebro
            superficie(mayor3D(baja(CT) > HUESO), esp .* [p p pz])];          % 6 craneo
    V = vertcat(r.fv.vertices);
    r.lim = [min(V, [], 1); max(V, [], 1)]' + [-5 5];                    % limites de los ejes
    r.invertirZ = false;  try, r.invertirZ = info.Transform.T(3, 3) < 0; catch, end

    % Cortes con hemorragia (medicos o nuestra), recortados alrededor del cerebro
    r.zz = find(squeeze(any(any(mask | sangre, 1), 2)));
    if isempty(r.zz), r.zz = round(size(CT, 3) / 2); end
    [f, c] = find(any(interior | mask, 3));
    f = max(min(f) - 20, 1):min(max(f) + 20, size(CT, 1));
    c = max(min(c) - 20, 1):min(max(c) + 20, size(CT, 2));
    r.ct2D  = im2uint8(mat2gray(CT(f, c, r.zz), [0 80]));               % ventana de cerebro
    r.med2D = uint8(mask(f, c, r.zz)) + uint8(fn(f, c, r.zz));          % 1 encontrada, 2 falso negativo
    r.nue2D = uint8(sangre(f, c, r.zz)) + uint8(fp(f, c, r.zz));        % 1 coincide, 2 falso positivo
    r.diceZ = arrayfun(@(z) calcDice(sangre(:,:,z), mask(:,:,z)), r.zz);
    [~, r.kIni] = max(squeeze(sum(sum(r.med2D > 0, 1), 2)));           % arranca en el corte con mas sangre
end

function fv = superficie(M, esp)
    % Superficie 3D suavizada de una mascara M, con vertices en mm
    fv = struct('faces', zeros(0, 3), 'vertices', zeros(0, 3));
    if ~any(M(:)), return; end
    [f, c, s] = ind2sub(size(M), find(M));                              % recorto a la caja del objeto
    lo = max([min(f) min(c) min(s)] - 2, 1);  hi = min([max(f) max(c) max(s)] + 2, size(M));
    sub = padarray(M(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), [1 1 1]);
    sup = isosurface(smooth3(double(sub), 'gaussian'), 0.5);
    if isempty(sup.vertices), return; end
    if size(sup.faces, 1) > 20000, sup = reducepatch(sup, 20000); end   % rota mas fluido
    fv.faces    = sup.faces;
    fv.vertices = (sup.vertices + lo([2 1 3]) - 3) .* esp([2 1 3]);     % indices -> mm (x = columna)
end

function visor(res)
    % Ventana: 3D Medicos | 3D Nuestra | cortes, con botones Anterior / Siguiente
    fig = figure('Name', 'Medicos vs Nuestra', 'NumberTitle', 'off', 'Color', 'w', ...
                 'Units', 'normalized', 'OuterPosition', [0.03 0.05 0.94 0.9]);
    d = struct('res', {res}, 'k', 1, 'ax3D', gobjects(1, 2), 'link', []);
    d.ax2D = [axes('Parent', fig, 'Position', [0.68 0.53 0.30 0.36]), ...
              axes('Parent', fig, 'Position', [0.68 0.11 0.30 0.36])];
    d.titulo = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.02 0.94 0.96 0.04], ...
                         'BackgroundColor', 'w', 'FontSize', 13, 'FontWeight', 'bold');
    uicontrol(fig, 'String', '< Anterior', 'Units', 'normalized', 'Position', [0.18 0.02 0.12 0.05], ...
              'FontSize', 11, 'Callback', @(s, ~) dibujarPaciente(ancestor(s, 'figure'), -1));
    uicontrol(fig, 'String', 'Siguiente >', 'Units', 'normalized', 'Position', [0.36 0.02 0.12 0.05], ...
              'FontSize', 11, 'Callback', @(s, ~) dibujarPaciente(ancestor(s, 'figure'), 1));
    d.slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.70 0.035 0.26 0.03], ...
                         'Min', 1, 'Max', 2, 'Value', 1, 'Callback', @(s, ~) dibujarCortes(ancestor(s, 'figure')));
    d.rot = rotate3d(fig);  d.rot.Enable = 'on';     % rotar con el mouse (solo los 3D)
    guidata(fig, d);
    dibujarPaciente(fig, 0);
end

function dibujarPaciente(fig, paso)
    % Pasa a otra CT (paso = -1 / +1) y dibuja los dos 3D y la barra de cortes
    d = guidata(fig);  d.k = mod(d.k - 1 + paso, numel(d.res)) + 1;  r = d.res{d.k};
    d.titulo.String = sprintf('CT %d  (%d de %d)      Dice 3D = %.3f      Dice medio por corte = %.3f', ...
                              r.datos(1), d.k, numel(d.res), r.datos(2), r.datos(3));
    vista = [-37.5 30];                              % mantengo el angulo de la CT anterior
    if isgraphics(d.ax3D(1)), [az, el] = view(d.ax3D(1));  vista = [az el];  delete(d.ax3D); end

    % Superficies: 1 medicos (rojo), 2 falso negativo (azul), 3 nuestra (verde),
    %              4 falso positivo (naranja), 5 cerebro, 6 craneo
    colores = [0.85 0.10 0.10; 0.20 0.40 1.00; 0.10 0.65 0.20; 1.00 0.55 0.00; ...
               0.50 0.50 0.65; 0.70 0.65 0.55];
    alfas   = [1 1 1 1 0.12 0.10];                   % opacidad
    capas   = {[6 5 1 2], [6 5 3 4]};                % que se dibuja en cada panel
    titulos = {{sprintf('Medicos: %.1f mL', r.datos(4)), sprintf('azul = falso negativo: %.1f mL', r.datos(7))}, ...
               {sprintf('Nuestra: %.1f mL', r.datos(5)), sprintf('naranja = falso positivo: %.1f mL', r.datos(6))}};
    for j = 1:2
        ax = axes('Parent', fig, 'Position', [0.01 + 0.33 * (j - 1), 0.11, 0.32, 0.78]);
        pintar(ax, r.fv(capas{j}), colores(capas{j}, :), alfas(capas{j}));
        light('Parent', ax, 'Position', [1 1 1], 'Style', 'infinite');
        light('Parent', ax, 'Position', [-1 -1 -1], 'Style', 'infinite');
        set(ax, 'XLim', r.lim(1, :), 'YLim', r.lim(2, :), 'ZLim', r.lim(3, :), ...
                'DataAspectRatio', [1 1 1], 'Visible', 'off');
        if r.invertirZ, ax.ZDir = 'reverse'; end
        view(ax, vista);  axis(ax, 'vis3d');  title(ax, titulos{j}, 'FontSize', 12);
        d.ax3D(j) = ax;
    end
    d.link = linkprop(d.ax3D, {'CameraPosition', 'CameraUpVector'});    % giran juntos

    n = numel(r.zz);  d.slider.Value = 1;            % la barra recorre los cortes con hemorragia
    if n > 1, set(d.slider, 'Max', n, 'Value', r.kIni, 'SliderStep', min([1 5] / (n - 1), 1), 'Enable', 'on');
    else,     set(d.slider, 'Max', 2, 'Enable', 'off'); end
    guidata(fig, d);
    dibujarCortes(fig);
end

function pintar(ax, fv, colores, alfas)
    % Dibuja cada superficie de fv que no este vacia, con su color y opacidad
    for i = find(arrayfun(@(s) ~isempty(s.faces), fv))
        patch('Parent', ax, 'Faces', fv(i).faces, 'Vertices', fv(i).vertices, 'FaceColor', colores(i, :), ...
              'FaceAlpha', alfas(i), 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
              'AmbientStrength', 0.5, 'SpecularStrength', 0.2);
    end
end

function dibujarCortes(fig)
    % Corte elegido con la barra: Medicos arriba (azul = falso negativo),
    % Nuestra abajo (naranja = falso positivo)
    d = guidata(fig);  r = d.res{d.k};
    i = min(max(round(d.slider.Value), 1), numel(r.zz));
    if strcmp(d.slider.Enable, 'on'), d.slider.Value = i; end
    img = rot90(r.ct2D(:,:,i));
    imshow(labeloverlay(img, rot90(r.med2D(:,:,i)), 'Colormap', [1 0 0; 0.2 0.4 1], 'Transparency', 0.5), ...
           'Parent', d.ax2D(1));
    imshow(labeloverlay(img, rot90(r.nue2D(:,:,i)), 'Colormap', [0 0.8 0; 1 0.55 0], 'Transparency', 0.5), ...
           'Parent', d.ax2D(2));
    title(d.ax2D(1), sprintf('Medicos | corte %d  (%d de %d)', r.zz(i), i, numel(r.zz)));
    title(d.ax2D(2), sprintf('Nuestra | Dice del corte = %.2f', r.diceZ(i)));
    setAllowAxesRotate(d.rot, d.ax2D, false);        % los cortes no se rotan
end