#!/usr/bin/env bash
# ==============================================================================
#  Entrypoint for the mechanogenomic-virtual-cell container.
#  Provides convenient subcommands; anything else is executed as-is.
# ==============================================================================
set -e

case "$1" in
  test)
    echo ">> Running validation suite (17 checks)"
    exec python test/test_virtual_cell.py
    ;;
  demo)
    echo ">> Running the VirtualCell demo"
    exec python src/mvirtual_cell.py
    ;;
  benchmark)
    exec python src/benchmark.py
    ;;
  sensitivity)
    exec python src/sensitivity.py
    ;;
  figures)
    echo ">> Regenerating figures"
    exec python results/make_figures.py
    ;;
  jupyter)
    echo ">> Starting Jupyter (open the printed URL; notebook in demo/)"
    exec jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser \
         --allow-root --NotebookApp.token=''
    ;;
  bash|sh)
    exec bash
    ;;
  *)
    # run any arbitrary command passed to the container
    exec "$@"
    ;;
esac
