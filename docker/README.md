# Reproducible environment (Docker)

This directory holds everything needed to run the project in a fully reproducible
container: pinned Python, pinned dependencies, headless rendering. The model,
tests, figures, renderers, and the interactive notebook run identically on any
machine — no local Python setup or version conflicts.

```
docker/
├── Dockerfile             # the image recipe
├── docker-compose.yml     # convenience wrapper (notebook, test)
├── docker-entrypoint.sh   # subcommands: test / demo / jupyter / bash / ...
└── README.md              # this guide
```
The `.dockerignore` lives at the repository root (Docker requires it there).

Published image: **`pipelinesinmegen/mvcell:v0.1`** (Docker Hub).

## Quick start (pull the prebuilt image)

No build needed — pull it straight from Docker Hub:

```bash
docker pull pipelinesinmegen/mvcell:v0.1
docker run --rm pipelinesinmegen/mvcell:v0.1            # run the validation suite
```

## Build it yourself

Build **from the repository root**, pointing at the Dockerfile with `-f`
(the `.` is the build context = the whole project):

```bash
git clone https://github.com/Danpc11/mechanogenomic-virtual-cell.git
cd mechanogenomic-virtual-cell
docker build -t pipelinesinmegen/mvcell:v0.1 -f docker/Dockerfile .
```

## Run

The default command runs the validation suite, which proves the build is correct
end-to-end:

```bash
docker run --rm pipelinesinmegen/mvcell:v0.1
```

Built-in subcommands:

```bash
docker run --rm pipelinesinmegen/mvcell:v0.1 demo          # VirtualCell demo
docker run --rm pipelinesinmegen/mvcell:v0.1 benchmark     # full model vs baselines
docker run --rm pipelinesinmegen/mvcell:v0.1 sensitivity   # local + global sensitivity
docker run --rm pipelinesinmegen/mvcell:v0.1 figures       # regenerate figures
docker run --rm -it pipelinesinmegen/mvcell:v0.1 bash      # interactive shell
```

Run any arbitrary command:

```bash
docker run --rm pipelinesinmegen/mvcell:v0.1 python src/mvirtual_cell.py
```

## Interactive notebook in the container

```bash
docker run --rm -p 8888:8888 pipelinesinmegen/mvcell:v0.1 jupyter
```

Then open the printed `http://127.0.0.1:8888/...` URL. The notebook is in
`demo/virtual_cell_demo.ipynb`.

With docker-compose (live-mounts the repo so edits persist), from the repo root:

```bash
docker compose -f docker/docker-compose.yml up jupyter
```

## Publishing to Docker Hub (maintainers)

```bash
# build (optionally tag :latest as well)
docker build -t pipelinesinmegen/mvcell:v0.1 -t pipelinesinmegen/mvcell:latest \
    -f docker/Dockerfile .

# authenticate, then push both tags
docker login
docker push pipelinesinmegen/mvcell:v0.1
docker push pipelinesinmegen/mvcell:latest
```

`pipelinesinmegen` must be your Docker Hub user or organization (not the GitHub
name), and you need push access to it.

## Notes

- Rendering runs **headless** (`PYVISTA_OFF_SCREEN=true`, `MPLBACKEND=Agg`), so
  the fluorescence, cross-section, and static PyVista renderers all work inside
  the container without a display.
- The interactive Trame 3D *app* (`--serve`) needs a browser and is best run on a
  host with a display; static PyVista export works in the container.
- The image is based on `python:3.11-slim-bookworm`. For archival reproducibility
  you can pin an exact digest: `python:3.11-slim-bookworm@sha256:<digest>`.
