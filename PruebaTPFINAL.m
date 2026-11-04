%% TP Final PDI - Segmentacion de hemorragia en CT (Seg-CQ500, NIfTI)
% Metodo relativo: en cada corte, la sangre es lo que es mas claro que el
% cerebro tipico de ESE corte (mediana) + un margen DELTA.
% Despues se limpia con morfologia y se compara con la mascara de los
% medicos usando el coeficiente de Dice.

clear; clc; close all

%% 1) Buscar los archivos del paciente dentro de Seg-CQ500
DATASET  = fullfile(pwd, 'Seg-CQ500');   % carpeta del dataset (al lado de este script)
PACIENTE = 1;                            % que paciente usar (1, 2, 3...)

% Busco todas las mascaras (una por paciente) en todas las subcarpetas
mascaras = dir(fullfile(DATASET, '**', '*mask*.nii*'));
if isempty(mascaras)
    error('No encontre mascaras dentro de %s. Revisar la carpeta.', DATASET);
end

% Muestro la lista de pacientes encontrados
fprintf('Pacientes encontrados: %d\n', numel(mascaras));
for k = 1:numel(mascaras)
    fprintf('  %2d) %s\n', k, mascaras(k).folder);
end

% Carpeta del paciente elegido
carpeta  = mascaras(PACIENTE).folder;
archMask = fullfile(carpeta, mascaras(PACIENTE).name);

% El CT es el .nii de la misma carpeta que no es mascara ni "brain"
todos = dir(fullfile(carpeta, '*.nii*'));
nombres = {todos.name};
esCT = ~contains(nombres, 'mask', 'IgnoreCase', true) & ...
       ~contains(nombres, 'brain', 'IgnoreCase', true);
if ~any(esCT)
    error('No encontre el CT en %s', carpeta);
end
archCT = fullfile(carpeta, nombres{find(esCT, 1)});

% Volumen sin craneo (solo se usa si usarBrain = true)
b = dir(fullfile(carpeta, '*brain*.nii*'));
archBrain = '';
if ~isempty(b), archBrain = fullfile(carpeta, b(1).name); end

fprintf('\nUsando paciente %d:\n  CT:      %s\n  Mascara: %s\n', ...
    PACIENTE, archCT, archMask);

%% 2) Parametros (lo unico que hay que tocar)
usarBrain = false;  % false = sacamos el craneo nosotros desde el CT (todo propio)
HUESO     = 100;    % HU: por encima es hueso
EROSION   = 3;      % radio para alejarse del borde del craneo
DELTA     = 15;     % HU por encima del cerebro para considerar sangre
TECHO     = 100;    % HU: por encima ya no es sangre (hueso, calcificacion)
R_OPEN    = 2;      % radio de la apertura morfologica
AREA_MIN  = 30;     % pixeles: objetos mas chicos se descartan

barrer    = true;   % true = probar varios DELTA
mostrar3D = false;  % true = reconstruccion 3D (es lenta)

%% 3) Leer los volumenes
info = niftiinfo(archCT);
CT   = double(niftiread(archCT));

% Paso a HU con la escala del encabezado (niftiread no la aplica sola)
pendiente = info.MultiplicativeScaling;
if pendiente == 0, pendiente = 1; end        % 0 significa "sin escala"
CT = CT * pendiente + info.AdditiveOffset;

mask = niftiread(archMask) > 0;              % true = hemorragia
if usarBrain
    if isempty(archBrain), error('No hay volumen "brain" para este paciente.'); end
    brain = double(niftiread(archBrain));
end

nZ = size(CT, 3);
fprintf('Volumen %d x %d x %d | HU min = %.0f, max = %.0f\n', size(CT), min(CT(:)), max(CT(:)));
% Control: el aire deberia estar cerca de -1000 HU y el hueso por encima de 300

if ~isequal(size(CT), size(mask))
    error('El CT y la mascara no tienen el mismo tamano.');
end

%% 4) Sacar el craneo y calcular el gris tipico del cerebro (corte por corte)
interior = false(size(CT));   % donde esta el cerebro
cf       = zeros(size(CT));   % CT filtrado con mediana
ref      = nan(nZ, 1);        % gris tipico del cerebro en cada corte

for z = 1:nZ
    corte = CT(:,:,z);

    % Interior del cerebro
    if usarBrain
        int = imfill(brain(:,:,z) > 0, 'holes');
    else
        craneo = corte > HUESO;                                  % 1. hueso
        disco  = imfill(imclose(craneo, strel('disk', 15)), 'holes'); % 2. cierro y relleno
        int    = disco & ~craneo;                                % 3. adentro sin hueso
    end
    if any(int(:))
        int = bwareafilt(int, 1);                % me quedo con el objeto mas grande
    end
    int = imerode(int, strel('disk', EROSION));  % me alejo del hueso
    interior(:,:,z) = int;

    % Filtro de mediana: saca ruido sin borrar bordes
    cf(:,:,z) = medfilt2(corte, [5 5]);

    % Gris tipico del cerebro de este corte
    if any(int(:))
        tmp    = cf(:,:,z);
        ref(z) = median(tmp(int));
    end
end

%% 5) Segmentar la sangre en todo el volumen
sangre = segmentarVolumen(cf, ref, interior, DELTA, TECHO, R_OPEN, AREA_MIN);

%% 6) Evaluar con Dice
conMask = squeeze(any(any(mask, 1), 2));     % cortes donde los medicos marcaron sangre

dice3D = calcDice(sangre, mask);             % volumen completo (cuenta TODOS los errores)

diceZ = nan(nZ, 1);                          % Dice por corte, solo cortes con hemorragia
for z = find(conMask)'
    diceZ(z) = calcDice(sangre(:,:,z), mask(:,:,z));
end

fprintf('\nDice 3D (volumen completo):        %.3f\n', dice3D);
fprintf('Dice medio (cortes con hemorragia): %.3f\n', mean(diceZ, 'omitnan'));
fprintf('Voxeles -> nuestro: %d | medicos: %d | acierto: %d\n', ...
    nnz(sangre), nnz(mask), nnz(sangre & mask));

figure
plot(find(conMask), diceZ(conMask), 'o-'); grid on
xlabel('Corte'); ylabel('Dice'); title('Dice por corte');

%% 7) Barrido de DELTA (evaluando el volumen completo)
if barrer
    deltas = 5:5:40;
    dd = zeros(size(deltas));
    for i = 1:numel(deltas)
        S = segmentarVolumen(cf, ref, interior, deltas(i), TECHO, R_OPEN, AREA_MIN);
        dd(i) = calcDice(S, mask);
        fprintf('DELTA = %2d -> Dice = %.3f\n', deltas(i), dd(i));
    end
    figure
    plot(deltas, dd, 'o-'); grid on
    xlabel('DELTA (HU)'); ylabel('Dice 3D'); title('Dice segun DELTA');
    [~, k] = max(dd);
    fprintf('Mejor DELTA: %d (Dice %.3f)\n', deltas(k), dd(k));
end

%% 8) Ver los cortes con mas hemorragia (rojo = medicos, verde = nuestro)
CTv = mat2gray(CT, [0 80]);                  % ventana de cerebro para mostrar

area = squeeze(sum(sum(mask, 1), 2));        % cuanta hemorragia hay en cada corte
[~, orden] = sort(area, 'descend');
mostrar = orden(1:min(6, nnz(area)));        % hasta 6 cortes con mas sangre

figure('Name', 'Segmentacion vs mascara')
for i = 1:numel(mostrar)
    z = mostrar(i);
    subplot(2, 3, i)
    imshow(rot90(CTv(:,:,z))); hold on
    visboundaries(rot90(mask(:,:,z)),   'Color', 'r');
    if any(sangre(:,:,z), 'all')
        visboundaries(rot90(sangre(:,:,z)), 'Color', 'g');
    end
    title(sprintf('Corte %d | Dice = %.2f', z, diceZ(z)));
end

%% 9) Reconstruccion 3D (opcional)
if mostrar3D
    reconstruir3D_suave('todo')
end

%% ============ FUNCIONES ============

function sangre = segmentarVolumen(cf, ref, interior, DELTA, TECHO, R_OPEN, AREA_MIN)
    % Aplica la segmentacion corte por corte a todo el volumen
    sangre = false(size(cf));
    for z = 1:size(cf, 3)
        if isnan(ref(z)), continue; end      % corte sin cerebro: lo salteo

        % 1. Mas claro que el cerebro (+DELTA), menos que el techo, solo adentro
        s = (cf(:,:,z) > ref(z) + DELTA) & (cf(:,:,z) < TECHO) & interior(:,:,z);

        % 2. Limpieza morfologica
        if R_OPEN > 0
            s = imopen(s, strel('disk', R_OPEN));   % borra puntitos e hilos
        end
        s = imclose(s, strel('disk', 3));           % une zonas cercanas
        s = imfill(s, 'holes');                      % tapa agujeros
        s = bwareaopen(s, AREA_MIN);                 % saca objetos chicos

        sangre(:,:,z) = s;
    end
end

function d = calcDice(A, B)
    % Dice: 1 = coincidencia perfecta, 0 = nada en comun
    d = 2 * nnz(A & B) / max(nnz(A) + nnz(B), 1);
end