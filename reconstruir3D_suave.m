function reconstruir3D_suave(tipo, carpeta)
% RECONSTRUIR3D_SUAVE  Reconstruccion 3D lisa (sin efecto "capas") con contornos Canny.
%
%   reconstruir3D_suave('craneo')       -> craneo   (CT.nii)
%   reconstruir3D_suave('cerebro')      -> cerebro  (CQ500-CT-4_brain.nii)
%   reconstruir3D_suave('hemorragia')   -> mascara del medico (ICH_mask.nii)
%   reconstruir3D_suave('todo')         -> las tres juntas
%   reconstruir3D_suave(tipo, carpeta)  -> si los .nii estan en otra carpeta
%
% Requiere Image Processing Toolbox (imresize3, imgaussfilt3, edge3: R2018b o mas nuevo).

    %% ---------- PARAMETROS ----------
    P.huesoHU     = 150;        % umbral de hueso (HU)
    P.paso        = 2;          % 1 = maxima calidad (mas lento), 2 = recomendado, 3 = rapido
    P.sigma       = 1.5;        % suavizado de la superficie (en voxeles nuevos). Mas alto = mas liso
    P.sigmaHem    = 1.0;        % suavizado de la hemorragia (menor para no borrarla)
    P.caras       = 0.5;        % fraccion de caras que se conservan

    P.canny       = true;       % true = resalta contornos con Canny 3D
    P.cannyUmbral = [0.10 0.25];% umbrales [bajo alto] de Canny (0..1). Subilos si hay mucho ruido
    P.cannySigma  = 1.5;        % suavizado previo de Canny
    P.fuerza      = 0.85;       % que tan oscuros se ven los contornos (0..1)
    P.profund     = 2;          % en el cerebro, mira Canny "unos voxeles adentro" de la superficie

    %% ---------- RUTAS ----------
    if nargin < 2
        carpeta = fileparts(mfilename('fullpath'));
        if isempty(carpeta), carpeta = pwd; end
    end
    archCT    = fullfile(carpeta, 'CT.nii');
    archBrain = fullfile(carpeta, 'CQ500-CT-4_brain.nii');
    archMask  = fullfile(carpeta, 'ICH_mask.nii');

    tipo = lower(tipo);
    if ~ismember(tipo, {'craneo','cerebro','hemorragia','todo'})
        error('tipo debe ser ''craneo'', ''cerebro'', ''hemorragia'' o ''todo''.');
    end

    %% ---------- GRILLA CUBICA ----------
    % Trabajo con dimensiones [Y X Z] (lo que espera isosurface) y re-muestreo a voxel cubico.
    info    = niftiinfo(archCT);
    vox     = info.PixelDimensions(1:3);   vox = vox([2 1 3]);
    sz      = info.ImageSize(1:3);         sz  = sz([2 1 3]);
    h       = mean(vox(1:2)) * P.paso;                 % lado del voxel nuevo (mm)
    newSize = max(round(sz .* vox / h), 4);

    figure('Color','k','Name',['3D suave: ' tipo]);
    hold on

    %% ---------- CRANEO ----------
    if ismember(tipo, {'craneo','todo'})
        CT = leer(archCT);
        BW = mayorObjeto(CT > P.huesoHU);
        if strcmp(tipo,'craneo'), alfa = 1; else, alfa = 0.08; end
        dibujar(BW, CT, [0 1000], [0.92 0.92 0.85], alfa, P.sigma, 0, newSize, P);
        clear CT BW
    end

    %% ---------- CEREBRO ----------
    if ismember(tipo, {'cerebro','todo'})
        B  = leer(archBrain);
        BW = false(size(B));
        for z = 1:size(B,3)
            BW(:,:,z) = imfill(B(:,:,z) > 0, 'holes');
        end
        BW = mayorObjeto(BW);
        if P.canny, CTc = leer(archCT); else, CTc = []; end   % Canny sobre el CT en ventana de cerebro
        if strcmp(tipo,'cerebro'), alfa = 1; else, alfa = 0.25; end
        dibujar(BW, CTc, [0 80], [0.90 0.60 0.65], alfa, P.sigma, P.profund, newSize, P);
        clear B BW CTc
    end

    %% ---------- HEMORRAGIA ----------
    if ismember(tipo, {'hemorragia','todo'})
        M = permute(niftiread(archMask) > 0, [2 1 3]);
        dibujar(M, [], [], [0.85 0.10 0.10], 1, P.sigmaHem, 0, newSize, P);
    end

    %% ---------- ESCENA ----------
    daspect([1 1 1]);                      % la grilla ya es cubica
    view(3); axis tight vis3d off
    camlight headlight; camlight('right');
    lighting gouraud; material([0.55 0.65 0.25 12]);
    rotate3d on
    title(['Reconstruccion 3D: ' tipo], 'Color','w');
end

%% ============ FUNCIONES LOCALES ============

function V = leer(arch)
    % Lee el .nii y lo deja en orden [Y X Z]
    V = permute(single(niftiread(arch)), [2 1 3]);
end

function BW = mayorObjeto(BW)
    % Deja solo la componente conexa mas grande
    cc = bwconncomp(BW, 26);
    if cc.NumObjects == 0, return; end
    [~, k] = max(cellfun(@numel, cc.PixelIdxList));
    BW = false(size(BW));
    BW(cc.PixelIdxList{k}) = true;
end

function S = distanciaConSigno(BW)
    % Mapa de distancias con signo por corte: positivo adentro, negativo afuera.
    % Interpolar ESTO entre cortes (y no la mascara) es lo que elimina los escalones.
    S = zeros(size(BW), 'single');
    for z = 1:size(BW,3)
        b = BW(:,:,z);
        if any(b(:))
            S(:,:,z) = bwdist(~b) - bwdist(b);
        else
            S(:,:,z) = -30;                  % corte vacio: "lejos" del objeto
        end
    end
end

function Es = bordesCanny(V, ventana, newSize, P)
    % Canny 3D sobre el CT (ventana de grises) ya llevado a la grilla cubica.
    V = imresize3(single(mat2gray(V, ventana)), newSize, 'linear');
    try
        E = edge3(V, 'approxcanny', P.cannyUmbral, P.cannySigma);
    catch
        % Plan B si tu MATLAB no tiene edge3: gradiente 3D con umbral
        G = imgradient3(imgaussfilt3(V, P.cannySigma));
        E = G > quantile(G(:), 0.92);
    end
    Es = imgaussfilt3(single(E), 1);         % engrosa un poco las lineas
end

function dibujar(BW, Vint, ventana, color, alfa, sigma, profund, newSize, P)
    if ~any(BW(:)), warning('Volumen vacio: no hay nada para dibujar.'); return; end

    % 1) Interpolo entre cortes y suavizo el mapa de distancias
    S = imresize3(distanciaConSigno(BW), newSize, 'linear');
    S = imgaussfilt3(S, sigma);

    % 2) Superficie en S = 0 (lisa, sin escalones)
    fv = isosurface(S, 0);
    fv = reducepatch(fv, P.caras);
    v  = fv.vertices;

    % 3) Contornos Canny pintados sobre la superficie
    if P.canny && ~isempty(Vint)
        Es = bordesCanny(Vint, ventana, newSize, P);

        % Direccion hacia adentro (S crece hacia adentro) para mirar un poco debajo de la superficie
        [gx, gy, gz] = gradient(S);
        d = [interp3(gx, v(:,1), v(:,2), v(:,3), 'linear', 0), ...
             interp3(gy, v(:,1), v(:,2), v(:,3), 'linear', 0), ...
             interp3(gz, v(:,1), v(:,2), v(:,3), 'linear', 0)];
        d = d ./ max(sqrt(sum(d.^2, 2)), eps);
        q = v + profund * d;

        Ev = interp3(Es, q(:,1), q(:,2), q(:,3), 'linear', 0);
        t  = P.fuerza * min(Ev / (quantile(Ev, 0.98) + eps), 1);   % 0 = sin borde, 1 = borde fuerte
        colorBorde = color * 0.15;
        C  = (1 - t) * color + t * colorBorde;                     % mezcla color base / color de borde

        p = patch(fv, 'FaceVertexCData', C, 'FaceColor', 'interp');
    else
        p = patch(fv, 'FaceColor', color);
    end
    p.EdgeColor = 'none';
    p.FaceAlpha = alfa;
end
