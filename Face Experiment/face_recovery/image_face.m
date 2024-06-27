function f = image_face(xj, f_id, h, w, img_title)

    Imagej = reshape(xj, [h, w]);
    f = figure;
    image(Imagej); 
    colormap(gray);
    f.Units = 'centimeters';
    f.Position = [1+3.6*fix(f_id/3)+10*mod(f_id,3) 15-4*fix(f_id/3) 3.6 4.8];
    title(img_title);
    set(gca, 'FontSize', 12, 'Fontname','times new Roman'); 
    xticks([]);
    yticks([]); 
    yticklabels([]);
    xticklabels([]);

end