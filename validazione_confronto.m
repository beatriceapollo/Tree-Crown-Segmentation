%% SCRIPT DI VALIDAZIONE DOPPIA: SPORCA (Baseline) vs PULITA (CHM Filtered)

clear; clc; close all;

fprintf('=== AVVIO VALIDAZIONE DOPPIA (OTTIMIZZATA) ===\n');

% 1. CARICAMENTO DATI

load('dati_validazione.mat');

fprintf('Estrazione dei %d centroidi previsti...\n', max(finalLabelMap(:)));

props = regionprops(finalLabelMap, 'Centroid', 'Area');

alberiValidi = [props.Area] >= 30; % Rimuove micro-frammenti

predCentroids = cat(1, props(alberiValidi).Centroid);

numPred = size(predCentroids, 1);

predX = predCentroids(:,1);

predY = predCentroids(:,2);

[xWorldExt, yWorldExt] = intrinsicToWorld(R, [cMin cMax], [rMin rMax]);

minX_W = min(xWorldExt); maxX_W = max(xWorldExt);

minY_W = min(yWorldExt); maxY_W = max(yWorldExt);

% 2. COSTRUZIONE DEI DUE GROUND TRUTH (SPORCO E PULITO)

dirtyGTPolygons = [];

cleanGTPolygons = [];

numTotalGT = length(truthShapes);

hWait1 = waitbar(0, 'Analisi del Ground Truth in corso...');

for i = 1:numTotalGT

if mod(i, 100) == 0, waitbar(i/numTotalGT, hWait1); end


polyX_W = truthShapes(i).X; polyY_W = truthShapes(i).Y;


% Se il poligono è dentro la nostra immagine

if any(polyX_W >= minX_W & polyX_W <= maxX_W) && any(polyY_W >= minY_W & polyY_W <= maxY_W)

[xPix, yPix] = worldToIntrinsic(R, polyX_W, polyY_W);

xRel = xPix - cMin + 1; yRel = yPix - rMin + 1;

xRel(isnan(xRel)) = []; yRel(isnan(yRel)) = [];


poligono = polyshape(xRel, yRel);


% Aggiungiamo SEMPRE al Ground Truth "Sporco"

dirtyGTPolygons = [dirtyGTPolygons; poligono];


% Controllo per il Ground Truth "Pulito" (Filtro CHM)

boxMinX = max(1, floor(min(xRel))); boxMaxX = min(size(chmCrop, 2), ceil(max(xRel)));

boxMinY = max(1, floor(min(yRel))); boxMaxY = min(size(chmCrop, 1), ceil(max(yRel)));

chmPatch = chmCrop(boxMinY:boxMaxY, boxMinX:boxMaxX);


% Se è alto almeno 2 metri, lo aggiungiamo ANCHE a quello "Pulito"

if max(chmPatch(:)) >= 2.0

cleanGTPolygons = [cleanGTPolygons; poligono];

end

end

end

close(hWait1);

numDirtyGT = length(dirtyGTPolygons);

numCleanGT = length(cleanGTPolygons);

fprintf('Etichette totali (Sporche): %d\n', numDirtyGT);

fprintf('Etichette valide (Pulite) : %d\n', numCleanGT);

fprintf('Etichette scartate (Cespugli/Basse): %d\n\n', numDirtyGT - numCleanGT);

% =========================================================================

% FUNZIONE LOCALE PER IL MATCHING VELOCE

% =========================================================================

function [TP, FP, FN, Precision, Recall, F1] = calcolaMetriche(predX, predY, gtPolygons)

numPred = length(predX);

numGT = length(gtPolygons);

TP = 0;

matchedGT = false(numGT, 1);

isTP = false(numPred, 1);


hWait = waitbar(0, 'Calcolo incroci spaziali...');

for g = 1:numGT

if mod(g, 100) == 0, waitbar(g/numGT, hWait); end


vX = gtPolygons(g).Vertices(:,1); vY = gtPolygons(g).Vertices(:,2);

minVx = min(vX); maxVx = max(vX); minVy = min(vY); maxVy = max(vY);


% Filtro veloce: cerca solo centroidi vicini al box

cands = find(predX >= minVx & predX <= maxVx & predY >= minVy & predY <= maxVy);


if ~isempty(cands)

in = inpolygon(predX(cands), predY(cands), vX, vY);

matchedIdx = cands(in);

if ~isempty(matchedIdx)

for m = 1:length(matchedIdx)

pidx = matchedIdx(m);

if ~isTP(pidx)

isTP(pidx) = true; matchedGT(g) = true; TP = TP + 1;

break;

end

end

end

end

end

close(hWait);


FP = numPred - TP;