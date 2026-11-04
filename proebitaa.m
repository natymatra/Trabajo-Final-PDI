%% TP Final PDI - Pasada periferica (recomendacion del profe)
% la idea (del profe): segmentar distinto segun el tipo de hemorragia.
%  - las de AFUERA (epidural, subdural) estan pegadas al craneo. Ahi el gris de
%    la sangre (50-80 HU) si se distingue del cerebro (~30 HU); el problema es
%    que nuestro metodo borra esa zona (me alejo 8 px del hueso).
%  - entonces hago una segunda pasada SOLO en una banda pegada al craneo, con su
%    propia regla, y la sumo a la pasada central (la de siempre).
% lo pruebo en los pacientes donde el metodo central no se podia aplicar
% (menos del 90% de la hemorragia dentro de la region cerebral).
% los parto en AJUSTE (elijo la regla) y EVALUACION (mido), como siempre.
%
% necesita elegibilidad_51.csv (lo genera PruebaTPFINAL.m con SOLO_DENTRO = true).
% el progreso queda en stats_periferica.mat: si lo corto, retoma donde quedo.

clear; clc; close all

%% 1) Parametros
DATASET       = fullfile(pwd, 'Seg-CQ500');   % carpeta del dataset (al lado de este script)
UMBRAL_DENTRO = 90;    % los pacientes con menos de esto dentro son los "perifericos"
SEMILLA       = 1;     % para partir en ajuste y evaluacion siempre igual

% pasada central: el metodo final de PruebaTPFINAL (no lo toco)
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

% pasada periferica
P.ANCHO_BANDA = 35;    % px: ancho de la banda pegada al craneo donde busco (~17 mm)
P.SEPARACION  = 2;     % px: no miro lo pegadisimo al hueso (el filtro de mediana lo mezcla con hueso)
P.CONTACTO    = 6;     % px: "toca el craneo" = esta a menos de esto del hueso
P.VOL_MIN_P   = 500;   % voxeles: objetos perifericos mas chicos se descartan
P.DELTAS_P    = [15 20 25];   % HU sobre el cerebro: los umbrales que pruebo en la banda
P.VERSION     = 1;     % si subo este numero, obligo a recalcular el cache

% reglas que pruebo para quedarme con un objeto periferico (no hace falta recalcular si las cambio)
CONTS_P     = [16 20 24];   % HU: contraste medio minimo
RADIOS_P    = [0 2 3];      % voxeles: radio minimo (0 = sin filtro de forma)
CONTACTOS_P = [0 30];       % %: cuanto del objeto tiene que tocar el craneo (0 = no lo pido)

%% 2) Que pacientes son "perifericos"
mascaras = dir(fullfile(DATASET, '**', '*mask*.nii*'));   % una mascara por paciente
if isempty(mascaras), error('No encontre mascaras dentro de %s.', DATASET); end
[~, orden] = sort({mascaras.folder});             % mismo orden que PruebaTPFINAL
mascaras = mascaras(orden);
nTot = numel(mascaras);
todosNombres = cell(nTot, 1);
for i = 1:nTot
    [~, todosNombres{i}] = fileparts(mascaras(i).folder);   % el nombre del paciente es el de su carpeta
end

if ~isfile('elegibilidad_51.csv')
    error('Falta elegibilidad_51.csv: corre primero PruebaTPFINAL.m con SOLO_DENTRO = true.');
end
E = readtable('elegibilidad_51.csv', 'TextType', 'string');   % % de la hemorragia dentro, de los 51
[esta, pos] = ismember(string(todosNombres), E.Caso);
pct = nan(nTot, 1);
pct(esta) = E.PctDentro(pos(esta));

perifericos = find(pct < UMBRAL_DENTRO);          % donde el metodo central no alcanzaba
n = numel(perifericos);
fprintf('Pacientes perifericos (menos de %d%% dentro): %d\n', UMBRAL_DENTRO, n);
if n < 2, error('Hay muy pocos pacientes perifericos.'); end

rng(SEMILLA);
esAjuste = false(n, 1);
esAjuste(randperm(n, floor(n / 2))) = true;       % la mitad ajusta, la otra mitad evalua

%% 3) Proceso los pacientes (esto es lo lento)
archCache = 'stats_periferica.mat';
cfg = P;
S = struct();
if isfile(archCache)
    L = load(archCache);
    if isequal(L.cfg, cfg)                        % mismos parametros: reuso lo guardado
        S = L.S;
        fprintf('Se retoma desde %s\n', archCache);
    else
        fprintf('Cambiaron los parametros: se recalcula todo.\n');
    end
end

nombres = cell(n, 1);  R = cell(n, 1);  ok = false(n, 1);
for i = 1:n
    m = mascaras(perifericos(i));
    nombres{i} = todosNombres{perifericos(i)};
    clave = matlab.lang.makeValidName(nombres{i});
    if ~isfield(S, clave)
        fprintf('(%d/%d) %s (%.0f%% dentro) ... ', i, n, nombres{i}, pct(perifericos(i)));
        tic
        try
            S.(clave) = procesarPeriferico(m, P);
            fprintf('listo (%.0f s)\n', toc);
        catch ME
            warning('Fallo %s: %s', nombres{i}, ME.message);   % si falla uno, sigo con los demas
            S.(clave) = [];
        end
        save(archCache, 'S', 'cfg');              % guardo cada vez que termino un paciente
    else
        fprintf('(%d/%d) %s: ya procesado\n', i, n, nombres{i});
    end
    R{i}  = S.(clave);
    ok(i) = ~isempty(R{i});
end

%% 4) Elijo la regla periferica usando SOLO los de ajuste
aj = find(esAjuste & ok);
if isempty(aj), error('Ningun paciente de ajuste se pudo procesar.'); end
nD = numel(P.DELTAS_P);  nC = numel(CONTS_P);  nR = numel(RADIOS_P);  nK = numel(CONTACTOS_P);
M = zeros(nD, nC, nR, nK);                        % dice medio de ajuste para cada combinacion
for a = 1:nD
    for c = 1:nC
        for b = 1:nR
            for q = 1:nK
                d = zeros(1, numel(aj));
                for j = 1:numel(aj)
                    d(j) = diceUnion(R{aj(j)}, a, CONTS_P(c), RADIOS_P(b), CONTACTOS_P(q), P.TECHO);
                end
                M(a, c, b, q) = mean(d);
            end
        end
    end
end

dCentralAj = zeros(1, numel(aj));                 % dice de ajuste con la pasada central sola
for j = 1:numel(aj)
    r = R{aj(j)};
    dCentralAj(j) = 2 * r.tpC / max(r.nC + r.nMask, 1);
end

fprintf('\nDice 3D medio de AJUSTE:\n');
fprintf('   solo pasada central          -> %.3f\n', mean(dCentralAj));
for a = 1:nD
    Ma = M(a, :, :, :);
    [m1, i1] = max(Ma(:));
    [c1, b1, q1] = ind2sub([nC nR nK], i1);
    fprintf('   central + periferica DELTA_P = %d -> %.3f (CONT_MIN = %d, RADIO_MIN = %d, CONTACTO >= %d%%)\n', ...
        P.DELTAS_P(a), m1, CONTS_P(c1), RADIOS_P(b1), CONTACTOS_P(q1));
end
[mejor, iB] = max(M(:));
[aB, cB, bB, qB] = ind2sub(size(M), iB);
fprintf('Elijo: DELTA_P = %d | CONT_MIN = %d | RADIO_MIN = %d | CONTACTO >= %d%% (Dice de ajuste %.3f)\n\n', ...
    P.DELTAS_P(aB), CONTS_P(cB), RADIOS_P(bB), CONTACTOS_P(qB), mejor);

%% 5) Resultados de todos los perifericos con esa regla
diceC = nan(n, 1);  sensC = nan(n, 1);  precC = nan(n, 1);   % solo central
diceP = nan(n, 1);  sensP = nan(n, 1);  precP = nan(n, 1);   % central + periferica
pctB  = nan(n, 1);  pctF  = nan(n, 1);  volM  = nan(n, 1);
for i = 1:n
    if ~ok(i), continue; end
    r = R{i};
    diceC(i) = 2 * r.tpC / max(r.nC + r.nMask, 1);
    sensC(i) = r.tpC / max(r.nMask, 1);
    precC(i) = r.tpC / max(r.nC, 1);
    [diceP(i), sensP(i), precP(i)] = diceUnion(r, aB, CONTS_P(cB), RADIOS_P(bB), CONTACTOS_P(qB), P.TECHO);
    pctB(i) = r.pctBanda;                         % cuanto de la hemorragia cae en la banda
    pctF(i) = r.pctFuera;                         % cuanto no cae ni en la region central ni en la banda
    volM(i) = r.nMask * r.volVox;
end

grupo = repmat("Evaluacion", n, 1);
grupo(esAjuste) = "Ajuste";
T = table(string(nombres), grupo, pct(perifericos), pctB, pctF, volM, diceC, diceP, sensC, sensP, precC, precP, ...
    'VariableNames', {'Caso', 'Grupo', 'PctDentro', 'PctBanda', 'PctFuera', 'Vol_medicos_mL', ...
                      'Dice_central', 'Dice_con_periferica', 'Sens_central', 'Sens_con', 'Prec_central', 'Prec_con'});
T = sortrows(T, 'Grupo');                         % primero los de ajuste
disp(T)
writetable(T, 'resultados_periferica.csv');      % lo guardo para el informe

okAj = esAjuste & ok;  okEv = ~esAjuste & ok;
fprintf('Dice 3D medio -> AJUSTE:     central %.3f | central + periferica %.3f\n', mean(diceC(okAj)), mean(diceP(okAj)));
fprintf('Dice 3D medio -> EVALUACION: central %.3f | central + periferica %.3f\n', mean(diceC(okEv)), mean(diceP(okEv)));
fprintf('Mediana (evaluacion):        central %.3f | central + periferica %.3f\n', median(diceC(okEv)), median(diceP(okEv)));
fprintf('Sensibilidad / precision medias (evaluacion): central %.3f / %.3f | con periferica %.3f / %.3f\n', ...
    mean(sensC(okEv)), mean(precC(okEv)), mean(sensP(okEv)), mean(precP(okEv)));
dif = diceP(okEv) - diceC(okEv);
fprintf('Pacientes de evaluacion donde mejora / empeora / queda igual: %d / %d / %d\n', ...
    nnz(dif > 0.01), nnz(dif < -0.01), nnz(abs(dif) <= 0.01));
fprintf('(PctFuera alto = la hemorragia cae en cortes de base o vertice que ignoro; ahi no hay pasada que la encuentre)\n');

figure('Name', 'Pasada periferica')
bar([diceC, diceP]); grid on                      % dos barras por paciente
xticks(1:n); xticklabels(nombres); xtickangle(45)
ylabel('Dice 3D'); ylim([0 1])
legend('solo central', 'central + periferica', 'Location', 'northoutside', 'Orientation', 'horizontal')
title('Pacientes con hemorragia periferica')

%% ============ FUNCIONES ============

function R = procesarPeriferico(m, P)
    % leo un paciente, hago la pasada central y guardo las estadisticas de los
    % objetos de la banda periferica para cada DELTA_P
    [CT, mask, info] = leerPaciente(m);
    [interior, intra, hueso] = calcularRegiones(CT, P);   % region central, todo lo de adentro del craneo, hueso
    [cf, ref] = filtrarYReferencia(CT, interior);  % CT filtrado y gris tipico del cerebro de cada corte
    clear CT                                       % libero memoria

    R.nMask  = nnz(mask);                          % voxeles de hemorragia segun los medicos
    R.volVox = prod(info.PixelDimensions(1:3)) / 1000;   % mL por voxel

    % 1) pasada central (la de siempre)
    cand    = candidatos(cf, ref, interior, P.DELTA, P);
    central = metodoCentral(cand, cf, ref, P);
    R.nC    = nnz(central);                        % voxeles que marca la central
    R.tpC   = nnz(central & mask);                 % y cuantos acierta

    % 2) banda pegada al craneo
    [banda, cercaHueso] = calcularBanda(intra, hueso, P);
    R.pctBanda = 100 * nnz(mask & banda) / max(R.nMask, 1);
    R.pctFuera = 100 * nnz(mask & ~interior & ~banda) / max(R.nMask, 1);

    % 3) objetos de la banda, para cada umbral
    Pp = P;  Pp.VOL_MIN = P.VOL_MIN_P;             % en la banda acepto objetos mas chicos
    R.objP = cell(1, numel(P.DELTAS_P));
    for d = 1:numel(P.DELTAS_P)
        candP     = candidatos(cf, ref, banda, P.DELTAS_P(d), Pp);   % mismo umbral con histeresis, pero en la banda
        R.objP{d} = estadisticasPerifericas(candP, mask, central, cf, ref, cercaHueso);
    end
end

function [interior, intra, hueso] = calcularRegiones(CT, P)
    % igual que calcularInterior, pero ademas devuelvo todo lo de adentro del
    % craneo (sin erosionar) y el hueso, que los necesito para la banda
    interior = false(size(CT));  intra = false(size(CT));  hueso = false(size(CT));
    for z = 1:size(CT, 3)
        corte  = CT(:,:,z);
        craneo = corte > P.HUESO;                                   % 1. todo lo que es hueso
        disco  = imfill(imclose(craneo, strel('disk', 15)), 'holes');   % 2. cierro el anillo del craneo y lo relleno
        dentro = mayorObjeto(disco & ~craneo);                      % 3 y 4. lo de adentro, sin hueso, el objeto mas grande
        int    = imerode(dentro, strel('disk', P.EROSION));         % 5. achico para alejarme del hueso
        if nnz(int) < P.MIN_CEREBRO                                 % corte de base o vertice: lo anulo
            int(:) = false;  dentro(:) = false;
        end
        interior(:,:,z) = int;  intra(:,:,z) = dentro;  hueso(:,:,z) = craneo;
    end
end

function [banda, cercaHueso] = calcularBanda(intra, hueso, P)
    % banda pegada al craneo: lo de adentro del craneo que esta a menos de
    % ANCHO_BANDA px del borde, sin lo pegadisimo al hueso
    banda = false(size(intra));  cercaHueso = false(size(intra));
    discoBanda = strel('disk', P.ANCHO_BANDA);
    discoSep   = strel('disk', P.SEPARACION);
    discoCont  = strel('disk', P.CONTACTO);
    for z = 1:size(intra, 3)
        dentro = intra(:,:,z);
        if ~any(dentro(:)), continue; end
        h = hueso(:,:,z);
        b = dentro & ~imerode(dentro, discoBanda); % anillo de ANCHO_BANDA px pegado al borde
        b = b & ~imdilate(h, discoSep);            % saco lo pegadisimo al hueso
        banda(:,:,z)      = b;
        cercaHueso(:,:,z) = imdilate(h, discoCont);   % zona a menos de CONTACTO px del hueso
    end
end

function central = metodoCentral(cand, cf, ref, P)
    % el metodo final de siempre: objetos gruesos, brillantes y sin picos de calcio
    central = false(size(cand));
    if ~any(cand(:)), return; end
    D      = bwdist(~cand);                        % radio en cada voxel
    cc     = bwconncomp(cand, 26);                 % objetos 3D
    refPac = median(ref, 'omitnan');               % HU tipico del cerebro del paciente
    for k = 1:cc.NumObjects
        idx       = cc.PixelIdxList{k};
        [~, ~, z] = ind2sub(size(cand), idx);
        contraste = double(cf(idx)) - ref(z);
        if max(D(idx)) >= P.RADIO_MIN && mean(contraste) >= P.CONT_MIN && max(contraste) + refPac < P.TECHO
            central(idx) = true;
        end
    end
end

function T = estadisticasPerifericas(sangre, mask, central, cf, ref, cercaHueso)
    % para cada objeto 3D de la banda guardo una fila:
    % [radio, voxeles nuevos, aciertos nuevos, contraste medio, pico HU, % que toca el craneo]
    % "nuevos" = los que la pasada central no marco (asi no los cuento dos veces)
    T = zeros(0, 6);
    if ~any(sangre(:)), return; end
    D  = bwdist(~sangre);
    cc = bwconncomp(sangre, 26);
    T  = zeros(cc.NumObjects, 6);
    for k = 1:cc.NumObjects
        idx       = cc.PixelIdxList{k};
        [~, ~, z] = ind2sub(size(sangre), idx);
        nuevo     = ~central(idx);
        contraste = double(cf(idx)) - ref(z);
        T(k,:) = [max(D(idx)), nnz(nuevo), nnz(mask(idx) & nuevo), mean(contraste), ...
                  max(double(cf(idx))), 100 * mean(cercaHueso(idx))];
    end
end

function [dice, sens, prec, vol] = diceUnion(r, a, contMin, radioMin, contactoMin, TECHO)
    % metricas de la pasada central + los objetos perifericos que pasan la regla
    T    = r.objP{a};
    sel  = T(:,1) >= radioMin & T(:,4) >= contMin & T(:,5) < TECHO & T(:,6) >= contactoMin;
    nS   = r.nC  + sum(T(sel, 2));                 % voxeles que marco en total
    TP   = r.tpC + sum(T(sel, 3));                 % voxeles acertados en total
    dice = 2 * TP / max(nS + r.nMask, 1);
    sens = TP / max(r.nMask, 1);
    prec = TP / max(nS, 1);
    vol  = nS * r.volVox;
end

function [CT, mask, info] = leerPaciente(m)
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