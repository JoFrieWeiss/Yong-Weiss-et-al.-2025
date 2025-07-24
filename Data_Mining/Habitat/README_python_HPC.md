
# Python Module and Package Setup on HPC

To use Python on the HPC system, follow these steps:

1. **Check available Python modules**

```bash
module avail python
```

2. **Load the desired Python version**

```bash
module load python/3.10.4
```

3. **Start a bash shell (if not already in bash)**

```bash
bash
```

4. **Check which Python packages are already installed**

```bash
pip list
```

5. **Install any missing packages**

```bash
pip install re pandas beautifulsoup4 requests time
```

> **Note:** If you do not have root or system-wide installation rights on the HPC, install packages only for your user:

```bash
pip install --user pandas beautifulsoup4 requests
```

This will ensure you have the necessary Python environment and packages for your work.

6. **Now use a shell script to start your scraping script on the HPC**


**Author:** Josefine Friederike Weiß
**Date:** 2025
