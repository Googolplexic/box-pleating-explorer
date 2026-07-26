# Box-Pleating Explorer

Tools and datasets for characterizing the space of locally-valid box-pleated
crease patterns. Crease patterns are sampled with Wave Function Collapse over a
catalogue of valid box-pleated vertex templates, with local validity (Maekawa,
Kawasaki, big-little-big) enforced by construction.

This project builds on the [`box-pleating`](https://github.com/Googolplexic/box-pleating) library.

## Status

Early development. The `mvp/` directory contains the initial degree-4, 4x4-grid
Wave Function Collapse sampler and rendering scripts used to validate the
approach.

## Layout

- `mvp/` — Wave Function Collapse sampler, vertex templates, and rendering
  scripts for a 4x4 grid
  - `templates.py` — degree-4 box-pleated vertex templates and D8 expansion
  - `wfc.py` — Wave Function Collapse sampler over the template vocabulary
  - `render_catalog.py`, `render_mosaic.py`, `verify_expansion.py`,
    `sanity_check_viz.py` — rendering and verification scripts

## Development

Install the sibling library in editable mode, then install this project:

```bash
python -m pip install -e ../box-pleating
python -m pip install -e ".[dev]"
```

Run the sampler self-checks and generate a sample mosaic:

```bash
python mvp/wfc.py
python mvp/render_mosaic.py
```

## License

MIT License - see the LICENSE file for details.
