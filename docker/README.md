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

## Build

Build **from the repository root**, pointing at the Dockerfile with `-f`
(the `.` is the build context = the whole project):

```bash
git clone https://github.com/Danpc11/mechanogenomic-virtual-cell.git
cd mechanogenomic-virtual-cell
docker build -t mvcell -f docker/Dockerfile .
```

## Run

The default command runs the validation suite, which proves the build is correct
end-to-end:

```bash
docker run --rm mvcell
```

Built-in subcommands:

```bash
docker run --rm mvcell demo          # VirtualCell demo
docker run --rm mvcell benchmark     # full model vs simple baselines
docker run --rm mvcell sensitivity   # local + global sensitivity
docker run --rm mvcell figures       # regenerate figures
docker run --rm -it mvcell bash      # interactive shell
```

Run any arbitrary command:

```bash
docker run --rm mvcell python src/mvirtual_cell.py
```

## Interactive notebook in the container

```bash
docker run --rm -p 8888:8888 mvcell jupyter
```

Then open the printed `http://127.0.0.1:8888/...` URL. The notebook is in
`demo/virtual_cell_demo.ipynb`.

With docker-compose (live-mounts the repo so edits persist), from the repo root:

```bash
docker compose -f docker/docker-compose.yml up jupyter
```

## Notes

- Rendering runs **headless** (`PYVISTA_OFF_SCREEN=true`, `MPLBACKEND=Agg`), so
  the fluorescence, cross-section, and static PyVista renderers all work inside
  the container without a display.
- The interactive Trame 3D *app* (`--serve`) needs a browser and is best run on a
  host with a display; static PyVista export works in the container.
- The image is based on `python:3.11-slim-bookworm`. For archival reproducibility
  you can pin an exact digest: `python:3.11-slim-bookworm@sha256:<digest>`.
