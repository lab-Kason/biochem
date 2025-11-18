# jupyter (workspace)

This folder contains Jupyter notebooks and helper scripts for enzyme kinetics analysis.

Virtual environment
- A venv has been created at: `.venv/` (Python interpreter: `.venv/bin/python`).
- To activate the venv (macOS / zsh):

```bash
source .venv/bin/activate
```

Install dependencies

To install the pinned dependencies (created from the venv used here):

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
```

Notes
- I renamed `import numpy as np.py` to `example_numpy_usage.py` to avoid shadowing the real `numpy` package when running code from this directory.
- If you prefer a different filename, rename `example_numpy_usage.py` accordingly.
