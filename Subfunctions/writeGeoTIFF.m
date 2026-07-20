function writeGeoTIFF(filename, data, epsg, corners, cmap, clim)
    % created with free version auf claude.ai and checked by T. Wunderlich,
    % CAU Kiel, 2026)

    % corners: 4x2 Matrix with UTM-Koordinaten of four corners of the area
    % corners = [xUpperLeft yUpperLeft; xUpperRight yUpperRight; xLowerLeft yLowerLeft; xLowerRight yLowerRight]
    %
    % clim: [min max] for colorscale
    % cmap: nx3 Colormap name (e.g. 'jet')
    
    if nargin < 4; cmap = 'jet'; end
    if nargin < 5; clim = [min(data(:)), max(data(:))]; end
    
    [rows, cols] = size(data);
    
    % Affin-Transformation aus Ecken berechnen
    % oben-links → oben-rechts gibt uns X-Richtung (pro Pixel)
    xOL = corners(1,1);  yOL = corners(1,2);
    xOR = corners(2,1);  yOR = corners(2,2);
    xUL = corners(3,1);  yUL = corners(3,2);
    
    % Pixel-Schrittweite in X- und Y-Richtung (als Vektor, wegen Rotation)
    a = (xOR - xOL) / cols;  % dX pro Pixel in Spaltenrichtung
    b = (yOR - yOL) / cols;  % dY pro Pixel in Spaltenrichtung
    c = (xUL - xOL) / rows;  % dX pro Pixel in Zeilenrichtung
    d = (yUL - yOL) / rows;  % dY pro Pixel in Zeilenrichtung
    
    % Transformationsmatrix (zeilenweise für TIFF-Tag)
    % [a b c xOL d e f yOL 0 0 0 0 0 0 0 1]
    transform = [a, b, 0, xOL, ...
                 c, d, 0, yOL, ...
                 0, 0, 0, 0,   ...
                 0, 0, 0, 1];
    
    % TIFF schreiben
    t = Tiff(filename, 'w');
    t.setTag('ImageWidth',          cols);
    t.setTag('ImageLength',         rows);
    t.setTag('Photometric',         Tiff.Photometric.MinIsBlack);
    t.setTag('BitsPerSample',       32);
    t.setTag('SampleFormat',        Tiff.SampleFormat.IEEEFP);
    t.setTag('SamplesPerPixel',     1);
    t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    t.setTag('Compression',         Tiff.Compression.None);
    
    % Rotation über ModelTransformationTag (statt Scale + Tiepoint)
    t.setTag('ModelTransformationTag', transform);
    
    % GeoKeys 
    geokeys = [
        1,    1, 0, 1;
        1024, 0, 1, 1;
        1025, 0, 1, 1;
        3072, 0, 1, epsg;
    ];
    header = [1, 1, 0, size(geokeys, 1)];
    t.setTag('GeoKeyDirectoryTag', [header; geokeys]);
    
    t.write(single(data));
    t.close();

    [fpath, fname, ~] = fileparts(filename);
    
    % QML mit NaN-Transparenz
    writeQML(fullfile(fpath, [fname '.qml']), cmap, clim);
    
    % AUX.XML für NoData = NaN
    writeAuxXML(fullfile(fpath, [fname '.tif.aux.xml']));

end


function writeQML(filename, cmapName, clim)
    if isstring(cmapName) || ischar(cmapName)
        cmap   = feval(cmapName, 16);
    else
        cmap=cmapName;
    end
    nColors = size(cmap, 1);
    
    fid = fopen(filename, 'w');
    fprintf(fid, '<!DOCTYPE qgis PUBLIC "http://mrcc.com/qgis.dtd" "SYSTEM">\n');
    fprintf(fid, '<qgis version="3.0">\n');
    fprintf(fid, '  <pipe>\n');
    fprintf(fid, '    <rasterrenderer type="singlebandpseudocolor" band="1" ');
    fprintf(fid, 'classificationMin="%f" classificationMax="%f">\n', clim(1), clim(2));
    fprintf(fid, '      <rasterTransparency>\n');
    fprintf(fid, '        <singleValuePixelList>\n');
    fprintf(fid, '          <pixelListEntry min="nan" max="nan" percentTransparent="100"/>\n');
    fprintf(fid, '        </singleValuePixelList>\n');
    fprintf(fid, '      </rasterTransparency>\n');
    fprintf(fid, '      <rastershader>\n');
    fprintf(fid, '        <colorrampshader colorRampType="INTERPOLATED" clip="0">\n');
    
    for i = 1:nColors
        val = clim(1) + (i-1) / (nColors-1) * (clim(2) - clim(1));
        r   = round(cmap(i,1) * 255);
        g   = round(cmap(i,2) * 255);
        b   = round(cmap(i,3) * 255);
        fprintf(fid, '          <item value="%f" color="#%02x%02x%02x" alpha="255" label="%.2f"/>\n', ...
                val, r, g, b, val);
    end
    
    fprintf(fid, '        </colorrampshader>\n');
    fprintf(fid, '      </rastershader>\n');
    fprintf(fid, '    </rasterrenderer>\n');
    fprintf(fid, '  </pipe>\n');
    fprintf(fid, '</qgis>\n');
    fclose(fid);
end

function writeAuxXML(filename)
    fid = fopen(filename, 'w');
    fprintf(fid, '<PAMDataset>\n');
    fprintf(fid, '  <PAMRasterBand band="1">\n');
    fprintf(fid, '    <NoDataValue>nan</NoDataValue>\n');
    fprintf(fid, '  </PAMRasterBand>\n');
    fprintf(fid, '</PAMDataset>\n');
    fclose(fid);
end