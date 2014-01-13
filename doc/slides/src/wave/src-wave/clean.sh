#!/bin/sh
for dir in `/bin/ls`; do
  if [ -d $dir ]; then
    cd $dir
    rm -rf *~ frame* *.html *.so *.o build *.pyc *.pyf wave2D_u0_loop_cy.c wave2D_u0_loop_c_cy.c tmp_build*
    cd ..
  fi
done
