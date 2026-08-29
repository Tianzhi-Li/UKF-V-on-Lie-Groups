function plot_cubic( vertices )    
fac = [1 2 3 4; 
       5 6 7 8; 
       1 2 5 8;
       3 4 7 6; 
       2 3 6 5; 
       1 4 7 8 ];
patch('Faces',fac,'Vertices',vertices,'FaceColor','white');  
grid off
hold off
axis equal
end