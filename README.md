# FTREX
Simple FTIR spectrum analyzer, and viewer with a nice GUI

**FTIR viewer and figure generator designed for R&D, data analysis, and publication‑ready graphics.**

---

## Table of Contents

1. [Description](#description)
2. [Requirements](#requirements)
3. [Initial Installation (first time only)](#initial-installation-first-time-only)
4. [Running the Application](#running-the-application)
5. [Interface and Operation Modes](#interface-and-operation-modes)

   * Basic Mode (Quick View)
   * Advanced Mode (R&D / Publication)
6. [Step‑by‑Step Workflow](#step-by-step-workflow)
7. [Peak Analysis & Annotations](#peak-analysis--annotations)
8. [Exporting Results](#exporting-results)
9. [Common Issues & Troubleshooting](#common-issues--troubleshooting)
10. [Contributing](#contributing)
11. [License & Contact](#license--contact)

---

## Description

**FTREX** is a **Streamlit‑based** application for fast, interactive visualization, processing, and high‑quality export of FTIR spectra. It includes both simple quick‑view features and advanced scientific processing tools.

---

## Requirements

* Python 3.8+ (recommended 3.9–3.11)
* pip
* Works on Windows / macOS / Linux

Supported file formats: `.dpt`, `.csv`, `.txt` (two numerical columns: wavenumber – intensity).

---

## Initial Installation (first time only) 

> **Important:** the `kaleido` package is required for exporting high‑resolution images.

1. Open your terminal (PowerShell, CMD, macOS Terminal, etc.).
2. Navigate to the project folder.
3. Run the following commands:

```bash
pip install -r requirements.txt
```

---

## Running the Application 

Inside the project folder, run:

```bash
streamlit run src/FTIR-SPECTRUM-EXPLORER.py
```

A browser tab will open automatically with the application.

---

## Interface and Operation Modes 

The app includes **two main modes**, selectable at the top of the left sidebar.

### A. Basic Mode (Quick View)

* Designed for fast visualization and new users.
* Simplified controls: automatic baseline correction & smoothing.
* **Automatic Stacking:** multiple loaded files are automatically offset to avoid overlap.
* Background grid enabled for easy coordinate reading.

### B. Advanced Mode (R&D / Publication)

* Designed for research use and generating publication‑ready figures.
* **Full preprocessing suite:**

  * Baseline correction: `Polynomial` and `Rubber-band`.
  * Signal transformation: Absorbance ↔ Transmittance (%).
  * Normalization: Area (integral) or Vector (L2).
* **Publication aesthetics:**

  * `Show Background Grid`: off by default for clean white figures.
  * `Hide Y-Ticks`: removes Y‑axis numbering when normalization makes it irrelevant.
  * `V-Lines`: enter wavenumbers (e.g., `1700, 2900`) to automatically draw labeled dashed vertical lines.

---

## Step‑by‑Step Workflow

### 1. Upload Spectrum Files 

Drag and drop `.dpt`, `.csv`, or `.txt` files into **Upload Spectrum Files** on the sidebar.
The program attempts to auto‑detect headers. Ensure your data columns are numeric.

### 2. Adjust Visualization

* In *Advanced Mode*, use **Manual Offset** to separate spectra vertically.
* Add peak markers by typing wavenumbers in **V-Lines**.
* Select a compound in **Annotate Compound** (e.g., *Urea*) to highlight characteristic regions.

### 3. Peak Analysis 

* Peaks are detected automatically.
* Use **Min Prominence** to control sensitivity.
* Filter peaks by region (e.g., `3000-2800`).
* A peak table is displayed at the bottom of the page.

### 4. Export Results 

* **Peak table:** click **Download CSV**.
* **Figures:** adjust **Export Scale** (recommended: `3` or `4` ≈ 300 DPI).
* Export formats available: **PNG**, **SVG**, **PDF**, **TIFF**.

---

## Peak Analysis & Annotations

* Automatic peak detection with adjustable width and prominence.
* Compound annotations highlight chemically relevant regions.
* All annotations are included in exports when using **Advanced Mode**.

---

## Common Issues & Troubleshooting

### ❗ "Export failed... kaleido"

You must install the image engine. Run:

```bash
pip install -U kaleido
```

### ❗ Mode switching IndexError

If it still appears, check your input files—some `.dpt` files contain corrupted or non-numeric lines.

### ❗ Blank plot

Likely due to unusual text headers in `.dpt` or `.txt` files. Clean or reformat the file.

### ❗ Slow or empty export

Lower the **Export Scale** and increase gradually. Large multi-spectrum figures require significant RAM.

---

## Contributing

1. Fork the repository.
2. Create a branch: `feature/your-feature`.
3. Submit a pull request describing your changes.

Ideas for contributions:

* Support new file formats (e.g., Bruker OPUS).
* Improve `.dpt` header handling.
* Expand the compound annotation database.

---

## License

This project is based on Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License.
You can find the license file on this repository. 

