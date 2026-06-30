%MATLAB's graphic system -- Incompatible renders or drivers in Ubuntu

%If OpenGL is buggy, try switching to Painters or ZBuffer (older renderers):

opengl('save', 'software')  % Forces MATLAB to use software OpenGL
set(0, 'DefaultFigureRenderer', 'painters')  % Alternative to OpenGL
