#!/bin/bash
set -x
export DISPLAY=:99.0
export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
which Xvfb
Xvfb :99 -screen 0 1024x768x24 > /dev/null 2>&1 &
sleep 3
set +x
exec "$@"
