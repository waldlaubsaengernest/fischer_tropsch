# Installation

## Clone the repository

```bash
git clone https://link
cd fischer_tropsch_1
```

## Create a Python virtual environment

Creating a virtual environment is recommended.

```bash
python3 -m venv venv
```

Activate it:

**Linux/macOS**

```bash
source venv/bin/activate
```

**Windows**

```bat
venv\Scripts\activate.bat
```

## Install the repository

```bash
pip install -e .
```

The folder data includes all the data and the scripts to create the figures in the paper. The structures are given in the folder data/structures.
