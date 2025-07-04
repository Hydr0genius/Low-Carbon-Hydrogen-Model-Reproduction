# Low-Carbon Hydrogen Optimization Model

This repository contains the code and inputs for simulating and optimizing low-carbon hydrogen production across multiple possibly Accounting Rules and Input Parameters. See Methods

## 📦 Requirements

- Python 3.10+ (recommended to use a virtual environment)
- Gurobi (if running the optimization; requires a license - free for academics. Alternative: use open source optimizers such as Pulp, which will increase optimization times)
- See `requirements.txt` for all packages and versions

## 🚀 Quick start

```bash
git clone https://github.com/Hydr0genius/Low-Carbon-Hydrogen-Model-Reproduction.git
cd <Low-Carbon-Hydrogen-Model-Reproduction>
python -m venv .venv
.venv\Scripts\activate          # macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
jupyter lab                     # or: jupyter notebook
