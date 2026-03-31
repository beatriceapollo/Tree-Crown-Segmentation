%% SCRIPT: CROP INTERATTIVO CON GROUND TRUTH SOVRAPPOSTO
clear; clc; close all;

% --- 1. PATH DEI FILE (Usa quelli del tassello 552000 o quello che preferisci) ---
imgPath = 'ABBY_rgb/2019_ABBY_3_555000_5067000_image.tif';
shpPath = 'ABBY_labels/2019_ABBY_3_555000_5067000_image.shp';

fprintf('Caricamento dati...\n');
[img, R] = readgeoraster(imgPath);
truthShapes = shaperead(shpPath);

% Normalizzazione veloce per la visualizzazione
imgDouble = double(img);
for b = 1:size(imgDouble, 3)
    band = imgDouble(:,:,b);
    low_in = prctile(band(band>0), 1);   
    high_in = prctile(band(band>0), 99); 
    band(band < low_in) = low_in;
    band(band > high_in) = high_in;
    imgDouble(:,:,b) = (band - low_in) / (high_in - low_in);
end

% --- 2. SELEZIONE INTERATTIVA DEL CROP (IN PIXEL) ---
figCrop = figure('Name', 'Seleziona Area', 'WindowState', 'maximized');
imshow(imgDouble);
title('DISEGNA UN RETTANGOLO su un gruppo di alberi (poi doppio click per confermare)');

% getrect restituisce [xMin_pixel, yMin_pixel, larghezza, altezza]
rect = getrect(figCrop); 
close(figCrop);

% Calcoliamo i limiti riga/colonna esatti
cMin = max(1, round(rect(1)));
rMin = max(1, round(rect(2)));
cMax = min(size(img,2), round(rect(1) + rect(3)));
rMax = min(size(img,1), round(rect(2) + rect(4)));

% Ritagliamo l'immagine
imgCrop = imgDouble(rMin:rMax, cMin:cMax, :);
fprintf('Area ritagliata: %d x %d pixel.\n', size(imgCrop,1), size(imgCrop,2));

% --- 3. CONVERSIONE E FILTRAGGIO DELLO SHAPEFILE ---
% Troviamo i limiti geografici (in metri) del nostro rettangolo di pixel
[xWorldExt, yWorldExt] = intrinsicToWorld(R, [cMin cMax], [rMin rMax]);
minX_W = min(xWorldExt); maxX_W = max(xWorldExt);
minY_W = min(yWorldExt); maxY_W = max(yWorldExt);

figure('Name', 'Ground Truth di Dettaglio', 'WindowState', 'maximized');
imshow(imgCrop);
hold on;

alberiTrovati = 0;

% Cicliamo su tutti i poligoni/box dello shapefile
for i = 1:length(truthShapes)
    polyX_World = truthShapes(i).X;
    polyY_World = truthShapes(i).Y;
    
    % Controlliamo se la bounding box cade (almeno in parte) nell'area selezionata
    if any(polyX_World >= minX_W & polyX_World <= maxX_W) && ...
       any(polyY_World >= minY_W & polyY_World <= maxY_W)
        
        % MAGIA: Trasformiamo le coordinate Mappa (metri) in coordinate Pixel
        [xPixelAssoluti, yPixelAssoluti] = worldToIntrinsic(R, polyX_World, polyY_World);
        
        % Trasliamo i pixel in base a dove abbiamo fatto il ritaglio
        xPixelRelativi = xPixelAssoluti - cMin + 1;
        yPixelRelativi = yPixelAssoluti - rMin + 1;
        
        % Disegniamo il poligono sull'immagine croppata!
        plot(xPixelRelativi, yPixelRelativi, 'y-', 'LineWidth', 1.5);
        
        alberiTrovati = alberiTrovati + 1;
    end
end

title(sprintf('Ground Truth sull''area selezionata (%d alberi)', alberiTrovati), 'FontSize', 14);
hold off;
fprintf('Completato! Mostrati %d alberi nell''area di interesse.\n', alberiTrovati);