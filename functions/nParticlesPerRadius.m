function nParticles = nParticlesPerRadius( particleMatrix )
%NPLASMARADIUS Calculate total number of plasma particles per radial bin

% Remove dimensions with size 1 and sum particles in all velocity bins within
%   current radial bin
nParticles = sum( squeeze(particleMatrix), 1 );
    
end

