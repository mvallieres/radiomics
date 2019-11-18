function [GLN] = getGLN_AertsMod(ROIonly,levels)
% USING CORRECTED METHOD (1 GLRLM MATRIX + NORMALIZED GLN)

[GLRLM] = getGLRLM(ROIonly,levels);
[textures] = getGLRLMtextures(GLRLM);
GLN = textures.GLN;

end