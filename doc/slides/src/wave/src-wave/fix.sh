#!/bin/sh
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'range\(0,\s+Nx\+1\)' 'Ix' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'range\(1,\s+Nx\)' 'Ix[1:-1]' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'range\(0,\s+Ny\+1\)' 'Iy' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'range\(1,\s+Ny\)' 'Iy[1:-1]' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'range\(1,\s+Nt\)' 'IIt' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'i = 0' 'i = Ix[0]' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'i = Nx' 'i = Ix[-1]' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'j = 0' 'j = Iy[0]' {} \;
find . -name '*.py' -exec grep -c Nx {} \; -exec doconce subst 'j = Ny' 'j = Iy[-1]' {} \;


