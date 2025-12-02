## -*- coding: utf-8 -*-
"""
FTIR Spectrum Explorer (Final Fixes)
Requires: streamlit, pandas, numpy, scipy, plotly, pillow, kaleido
"""

from __future__ import annotations
import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
from scipy.signal import find_peaks, savgol_filter
# NOTA: trapz está en numpy, no en scipy.integrate en versiones nuevas
import plotly.graph_objects as go
import plotly.io as pio
from functools import lru_cache
import io, json
from PIL import Image

# ---------------------------
# 0. CONSTANTS & PALETTES
# ---------------------------

CURATED_HIGH_CONTRAST_PALETTE = [
    '#000000', '#1E88E5', '#D81B60', '#43A047', '#FB8C00', '#8E24AA',
    '#FDD835', '#00ACC1', '#6D4C41', '#C62828', '#2E7D32', '#1565C0'
]

# ---------------------------
# 1. COMPOUNDS JSON LOADER
# ---------------------------

@st.cache_data
def load_compounds_json(path: str = "compounds.json") -> Dict[str, List[List[float]]]:
    """Load compound signature ranges from the external JSON file."""
    try:
        with open(path, 'r') as f:
            data = json.load(f)
            # Convert list of lists to list of tuples for consistency
            return {k: [tuple(r) for r in v] for k, v in data.items()}
    except (FileNotFoundError, json.JSONDecodeError):
        return {}

def parse_ranges(input_str: str) -> List[Tuple[float, float]]:
    """Parse comma/semicolon-separated ranges (e.g., '3400-3300, 1700')."""
    if not input_str:
        return []
    
    ranges = []
    input_str = input_str.replace(';', ',').replace(' ', '')
    
    for item in input_str.split(','):
        if not item:
            continue
        try:
            if '-' in item:
                start, end = map(float, item.split('-'))
                ranges.append((max(start, end), min(start, end)))
            else:
                wavenumber = float(item)
                ranges.append((wavenumber + 5, wavenumber - 5)) 
        except ValueError:
            continue
    return ranges

# ---------------------------
# 2. I/O & PARSING
# ---------------------------

@st.cache_data
def read_dpt_bytes(file_bytes: bytes, filename: str) -> Optional[pd.DataFrame]:
    """Robustly read DPT/TXT/CSV files."""
    try:
        s = str(file_bytes, 'utf-8', errors='replace')
        data = io.StringIO(s)
        sep = ',' if ',' in s else (';' if ';' in s else r'\s+')

        df = pd.read_csv(data, sep=sep, header=None, comment='#', skipinitialspace=True)
        
        if df.shape[1] < 2:
            return None
        
        # Auto-detect Wavenumber
        if (df.iloc[:, 0].max() - df.iloc[:, 0].min()) > (df.iloc[:, 1].max() - df.iloc[:, 1].min()):
            wavenumber = df.iloc[:, 0]
            intensity = df.iloc[:, 1]
        else:
            wavenumber = df.iloc[:, 1]
            intensity = df.iloc[:, 0]

        result_df = pd.DataFrame({
            'Wavenumber': wavenumber,
            'Intensity': intensity
        }).dropna().sort_values(by='Wavenumber', ascending=False).reset_index(drop=True)
        
        return result_df
    except Exception:
        return None

# ---------------------------
# 3. PREPROCESSING
# ---------------------------

def rubber_band_baseline(df: pd.DataFrame) -> pd.DataFrame:
    """Simple simplified baseline correction."""
    y = df['Intensity'].values
    corrected_intensity = y - y.min()
    return df.assign(Intensity=corrected_intensity)

def baseline_correct(df: pd.DataFrame, method: str) -> pd.DataFrame:
    if method == "rubber":
        return rubber_band_baseline(df)
    return df

def smooth_signal(df: pd.DataFrame, window: int, poly: int) -> pd.DataFrame:
    if window % 2 == 0:
        window += 1
    smoothed = savgol_filter(df['Intensity'].values, window, poly)
    return df.assign(Intensity=smoothed)

def transform_signal(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    intensity = df['Intensity'].values
    if mode == "Transmittance (%) -> Absorbance":
        # Avoid log(0)
        intensity = np.where(intensity <= 0, 1e-6, intensity)
        T = intensity / 100
        A = -np.log10(T)
        return df.assign(Intensity=A)
    elif mode == "Absorbance -> Transmittance (%)":
        T = 10**(-intensity)
        return df.assign(Intensity=T * 100)
    return df

def normalize_signal(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    y = df['Intensity'].values
    if mode == "Min-Max (0-1)":
        mn, mx = y.min(), y.max()
        if mx - mn != 0:
            y_norm = (y - mn) / (mx - mn)
        else:
            y_norm = y
    elif mode == "Divide by max":
        mx = y.max()
        y_norm = y / mx if mx != 0 else y
    elif mode == "Vector (L2)":
        norm = np.linalg.norm(y)
        y_norm = y / norm if norm != 0 else y
    elif mode == "Area (Integral)":
        area = np.trapz(y, x=df['Wavenumber'].values) # Use np.trapz
        y_norm = y / abs(area) if area != 0 else y
    else:
        return df
    return df.assign(Intensity=y_norm)

# ---------------------------
# 4. PEAK ANALYSIS (CORREGIDO)
# ---------------------------

def detect_peaks(df: pd.DataFrame, type: str, prominence: float, distance_cm: float) -> pd.DataFrame:
    """
    Detect peaks using find_peaks. Fixed IndexError by ensuring strict integer arrays.
    """
    y = df['Intensity'].values
    x = df['Wavenumber'].values
    
    total_range = x.max() - x.min()
    if total_range == 0:
        return pd.DataFrame()
        
    points_per_cm = len(x) / total_range
    distance_points = int(distance_cm * points_per_cm) if distance_cm > 0 else 1
    
    # Initialize as empty INTEGER arrays
    indices_up = np.array([], dtype=int)
    indices_down = np.array([], dtype=int)
    
    if type == "up" or type == "both":
        idx, _ = find_peaks(y, prominence=abs(y.max() - y.min()) * prominence, distance=distance_points)
        indices_up = idx.astype(int)

    if type == "down" or type == "both":
        idx, _ = find_peaks(-y, prominence=abs(y.max() - y.min()) * prominence, distance=distance_points)
        indices_down = idx.astype(int)

    # Concatenate safely
    peak_indices = np.concatenate([indices_up, indices_down])
    
    if len(peak_indices) == 0:
        return pd.DataFrame()
        
    # Use indices to get data
    peak_wavenumbers = x[peak_indices]
    peak_intensities = y[peak_indices]
    peak_types = ['up'] * len(indices_up) + ['down'] * len(indices_down)
    
    return pd.DataFrame({
        'Wavenumber': peak_wavenumbers,
        'Intensity': peak_intensities,
        'Type': peak_types
    }).sort_values(by='Wavenumber', ascending=False)

def generate_publication_table(peak_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    if not peak_data:
        return pd.DataFrame()
    
    # Build list of series for concatenation
    cols = {}
    max_len = 0
    
    for sample, df in peak_data.items():
        if not df.empty:
            # Sort by position
            peaks = df['Wavenumber'].round(0).astype(int).tolist()
            cols[sample] = peaks
            if len(peaks) > max_len:
                max_len = len(peaks)
        else:
            cols[sample] = []

    # Pad with empty strings/NaN
    for k in cols:
        cols[k] = cols[k] + [''] * (max_len - len(cols[k]))
        
    return pd.DataFrame(cols)
 # ---------------------------
# 5. UI/CONTROLS
# ---------------------------

# ---------------------------
# 5. UI/CONTROLS
# ---------------------------

def sidebar_controls(compounds: Dict[str, List[Tuple[float, float]]]):
    st.sidebar.header("🔬 FTIR Spectrum Explorer")
    
    analysis_mode = st.sidebar.radio(
        "Select Interface Mode", 
        ["Basic (Quick View)", "Advanced (R&D / Publication)"], 
        index=0
    )
    
    st.sidebar.markdown("---")
    uploaded_files = st.sidebar.file_uploader("Upload Spectrum Files", type=["dpt", "txt", "csv"], accept_multiple_files=True)
    
    controls = {"uploaded_files": uploaded_files, "mode": analysis_mode}
    
    st.sidebar.markdown("---")

    if analysis_mode == "Basic (Quick View)":
        st.sidebar.subheader("Core Preprocessing")
        controls["baseline_method"] = st.sidebar.selectbox("Baseline Correction", options=["rubber", "None"], index=0)
        controls["smoothing_enabled"] = st.sidebar.checkbox("Enable Smoothing", value=True)
        controls["normalization_mode"] = st.sidebar.selectbox("Normalization", options=["None", "Divide by max", "Min-Max (0-1)"], index=0)
        
        st.sidebar.subheader("Visualization")
        controls["auto_offset"] = st.sidebar.checkbox("Automatic Stacking", value=True)
        controls["manual_offset"] = 0.0 
        
        controls["title"] = st.sidebar.text_input("Plot Title", "FTIR Spectra Analysis")
        controls["y_axis_label"] = st.sidebar.selectbox("Y-axis Label", options=["Absorbance (a.u.)", "Transmittance (%)"], index=0)
        controls["export_scale"] = st.sidebar.slider("Export Scale", 1, 5, 2)
        
        # Defaults for Basic Mode
        controls["theme"] = "Light"
        controls["show_grid"] = True # Basic users usually like grids for reading
        controls["legend_position"] = "Right"
        controls["sg_window"], controls["sg_poly"] = 11, 2
        controls["transform_mode"] = "None"
        controls["hide_yaxis_ticks"] = False
        controls["show_peaks"] = True 
        controls["min_prominence"], controls["min_distance_cm"] = 0.05, 10.0
        controls["peak_filter_input"] = ""
        controls["compound_choice"] = "None"
        controls["vline_input"] = ""
        controls["peak_decimals"] = 0
        controls["peak_type"] = "up" 

    elif analysis_mode == "Advanced (R&D / Publication)":
        st.sidebar.subheader("Advanced Preprocessing")
        controls["transform_mode"] = st.sidebar.selectbox("Signal Transformation", options=["None", "Transmittance (%) -> Absorbance", "Absorbance -> Transmittance (%)"], index=0)
        controls["baseline_method"] = st.sidebar.selectbox("Baseline Correction", options=["rubber", "polynomial", "iterative", "None"], index=0)
        controls["smoothing_enabled"] = st.sidebar.checkbox("Enable Smoothing", value=True)
        
        if controls["smoothing_enabled"]:
            controls["sg_window"] = st.sidebar.slider("Smoothing Window", 3, 51, 11, step=2) 
            controls["sg_poly"] = st.sidebar.slider("Polynomial Order", 1, 5, 2)
        else:
            controls["sg_window"], controls["sg_poly"] = 11, 2

        controls["normalization_mode"] = st.sidebar.selectbox("Normalization Mode", options=["None", "Min-Max (0-1)", "Divide by max", "Vector (L2)", "Area (Integral)"], index=0)
        
        st.sidebar.subheader("Visualization")
        controls["auto_offset"] = st.sidebar.checkbox("Automatic Stacking", value=False)
        controls["manual_offset"] = st.sidebar.slider("Manual Offset", 0.0, 2.0, 0.0, 0.05) if not controls["auto_offset"] else 0.0
        
        controls["title"] = st.sidebar.text_input("Plot Title", "FTIR Spectra — Publication Ready")
        c1, c2 = st.sidebar.columns(2)
        controls["y_axis_label"] = c1.selectbox("Y-axis Label", options=["Absorbance (a.u.)", "Transmittance (%)", "Custom"], index=0)
        controls["hide_yaxis_ticks"] = c2.checkbox("Hide Y-Ticks", value=False)
        if controls["y_axis_label"] == "Custom":
            controls["y_axis_label"] = st.sidebar.text_input("Custom Label", "Intensity")

        # NEW: Grid Control
        controls["show_grid"] = st.sidebar.checkbox("Show Background Grid", value=False) # Default False for clean look
        controls["theme"] = st.sidebar.selectbox("Theme", options=["Light", "Scientific", "Dark"], index=1) # Default Scientific (Cleaner)
        controls["legend_position"] = st.sidebar.selectbox("Legend", options=["Right", "Inside Top-Right", "Inside Bottom-Left"], index=0)
        
        st.sidebar.subheader("Peaks & Annotations")
        controls["show_peaks"] = st.sidebar.checkbox("Show Peak Markers", value=True)
        controls["compound_choice"] = st.sidebar.selectbox("Annotate Compound", options=["None"] + list(compounds.keys()), index=0)
        
        # V-Lines input
        controls["vline_input"] = st.sidebar.text_input("V-Lines (e.g. 1700, 2900)", "")
        
        controls["peak_type"] = st.sidebar.selectbox("Peak Type", options=["up", "down", "both"], index=0)
        controls["min_prominence"] = st.sidebar.slider("Min Prominence", 0.0, 1.0, 0.05, 0.01)
        controls["min_distance_cm"] = st.sidebar.slider("Min Distance (cm⁻¹)", 5.0, 50.0, 10.0)
        controls["peak_filter_input"] = st.sidebar.text_input("Filter Ranges", "")
        controls["peak_decimals"] = st.sidebar.number_input("Peak Decimals", 0, 3, 1)
        controls["export_scale"] = st.sidebar.slider("Export Scale", 1, 5, 4)

    return controls

# ---------------------------
# 6. PLOTTING
# ---------------------------

def build_figure(data: Dict[str, pd.DataFrame], peak_data: Dict[str, pd.DataFrame], controls: Dict[str, Any], compounds: Dict[str, List[Tuple[float, float]]]) -> go.Figure:
    fig = go.Figure()
    
    # Theme Setup
    if controls["theme"] == "Dark":
        template = "plotly_dark"; line_color = '#ffffff'; grid_color='#444444'
    elif controls["theme"] == "Scientific":
        template = "plotly_white"; line_color = '#000000'; grid_color='#E5E5E5'
    else: # Light
        template = "plotly_white"; line_color = '#000000'; grid_color='#E5E5E5'

    legend_map = {
        "Right": dict(yanchor="top", y=1, xanchor="left", x=1.02),
        "Inside Top-Right": dict(yanchor="top", y=0.99, xanchor="right", x=0.99),
        "Inside Bottom-Left": dict(yanchor="bottom", y=0.01, xanchor="left", x=0.01)
    }
    
    current_offset = 0.0
    color_idx = 0
    prev_max = 0.0
    filenames = list(data.keys())

    for i, fname in enumerate(filenames):
        df = data[fname]
        if df.empty: continue
        
        # Offset Logic
        if controls["auto_offset"] and i > 0:
            gap = prev_max * 0.15
            current_offset += prev_max + gap
        elif not controls["auto_offset"]:
            current_offset = controls["manual_offset"] * i
        
        y_shifted = df['Intensity'].values + current_offset
        prev_max = df['Intensity'].max() # Track max for next iteration

        # Trace
        fig.add_trace(go.Scatter(
            x=df['Wavenumber'], y=y_shifted, mode='lines', name=fname,
            line=dict(color=CURATED_HIGH_CONTRAST_PALETTE[color_idx % len(CURATED_HIGH_CONTRAST_PALETTE)], width=2.0)
        ))
        
        # Peaks
        if controls["show_peaks"] and fname in peak_data and not peak_data[fname].empty:
            pk = peak_data[fname]
            fig.add_trace(go.Scatter(
                x=pk['Wavenumber'], y=pk['Intensity'] + current_offset,
                mode='markers', showlegend=False,
                marker=dict(color=CURATED_HIGH_CONTRAST_PALETTE[color_idx % len(CURATED_HIGH_CONTRAST_PALETTE)], size=7, symbol='circle-open')
            ))
        color_idx += 1

    # Annotations (V-LINES WITH LABELS)
    # Parse ranges but take the first value as the line position
    vlines = parse_ranges(controls["vline_input"])
    for x_start, x_end in vlines:
        # Draw the line
        fig.add_vline(x=x_start, line_width=1, line_dash="dash", line_color="gray", opacity=0.8)
        
        # Draw the Label (Number) at the top
        fig.add_annotation(
            x=x_start, y=1, yref="paper", # Position at top of plot area
            text=f"<b>{int(x_start)}</b>",
            showarrow=False,
            font=dict(size=10, color="#555555"),
            yshift=10 # Shift slightly up
        )
    
    # Compound shading
    if controls["compound_choice"] != "None":
        for rng in compounds.get(controls["compound_choice"], []):
            fig.add_vrect(x0=rng[1], x1=rng[0], fillcolor="orange", opacity=0.1, line_width=0)

    # Layout
    show_grid = controls.get("show_grid", False) # Default to false if missing
    
    fig.update_layout(
        template=template,
        title=dict(text=controls["title"], font=dict(size=18)),
        xaxis=dict(
            title='Wavenumber (cm<sup>-1</sup>)', 
            linecolor=line_color, 
            mirror=True, 
            showgrid=show_grid, # CONTROLLED BY CHECKBOX
            gridcolor=grid_color,
            autorange='reversed'
        ),
        yaxis=dict(
            title=controls["y_axis_label"],
            linecolor=line_color, 
            mirror=True, 
            showgrid=False # Usually disabled in Y for FTIR
        ),
        legend=legend_map.get(controls["legend_position"], {}),
        hovermode="x unified", height=600, margin=dict(l=50, r=50, t=60, b=50)
    )

    if controls["hide_yaxis_ticks"]:
        fig.update_layout(yaxis=dict(showticklabels=False))
        
    return fig

# ---------------------------
# 7. EXPORT HELPERS
# ---------------------------
def fig_to_bytes(fig: go.Figure, fmt: str = "png", scale: int = 2) -> bytes:
    try:
        # This requires 'kaleido' to be installed in the environment
        return pio.to_image(fig, format=fmt, scale=scale)
    except Exception as e:
        # Generic fallback message if import fails
        raise RuntimeError(f"Export failed. Please ensure 'kaleido' is installed (pip install -U kaleido). Error: {str(e)}")
 # ---------------------------
# 8. MAIN
# ---------------------------

def main():
    st.set_page_config(page_title="FTIR Explorer", layout="wide", page_icon="🧪")
    st.title("🔬 FTIR Spectrum Explorer — Enterprise")

    compounds = load_compounds_json("compounds.json")
    controls = sidebar_controls(compounds)
    uploaded = controls["uploaded_files"]
    
    if not uploaded:
        st.info("Please upload spectrum files (.dpt/.txt/.csv).")
        return 

    processed_data = {}
    peak_data = {}
    
    for file in uploaded:
        df = read_dpt_bytes(file.getvalue(), file.name)
        if df is None: continue
            
        # Pipeline
        df = transform_signal(df, controls.get("transform_mode", "None"))
        df = baseline_correct(df, controls["baseline_method"])
        if controls["smoothing_enabled"]:
            df = smooth_signal(df, controls["sg_window"], controls["sg_poly"])
        df = normalize_signal(df, controls["normalization_mode"])
        
        processed_data[file.name] = df
        
        # Peak Detection (Always runs if show_peaks is True, which is default for Basic now)
        if controls.get("show_peaks", False) or controls["mode"] == "Advanced (R&D / Publication)":
            peaks = detect_peaks(df, controls["peak_type"], controls["min_prominence"], controls["min_distance_cm"])
            
            # Filter
            fr = parse_ranges(controls["peak_filter_input"])
            if fr and not peaks.empty:
                mask = pd.Series(False, index=peaks.index)
                for s, e in fr:
                    mask |= (peaks['Wavenumber'] <= s) & (peaks['Wavenumber'] >= e)
                peaks = peaks[mask]
            
            peak_data[file.name] = peaks

    if not processed_data:
        st.error("No data processed.")
        return

    # Plot
    fig = build_figure(processed_data, peak_data, controls, compounds)
    st.plotly_chart(fig, use_container_width=True)

    # Export Table
    st.markdown("---")
    pub_df = generate_publication_table(peak_data)
    
    c1, c2 = st.columns([1, 1])
    
    with c1:
        st.subheader("📊 Peak Table")
        if not pub_df.empty:
            # Decimal formatting
            d = controls.get("peak_decimals", 0)
            st.dataframe(pub_df, use_container_width=True)
            st.download_button("📥 Download CSV", pub_df.to_csv(index=False).encode("utf-8"), "ftir_peaks.csv", "text/csv")
        else:
            st.info("No peaks detected.")

    with c2:
        st.subheader("🖼 Export Figure")
        xc1, xc2, xc3, xc4 = st.columns(4)
        try:
            if xc1.button("PNG"):
                st.download_button("Download PNG", fig_to_bytes(fig, "png", controls["export_scale"]), "plot.png", "image/png")
            if xc2.button("SVG"):
                st.download_button("Download SVG", fig_to_bytes(fig, "svg", controls["export_scale"]), "plot.svg", "image/svg+xml")
            if xc3.button("PDF"):
                st.download_button("Download PDF", fig_to_bytes(fig, "pdf", controls["export_scale"]), "plot.pdf", "application/pdf")
            if xc4.button("TIFF"):
                # TIFF requires PIL conversion usually, direct output might depend on backend
                # Simple fallback to PNG->TIFF if direct fails
                png_dat = fig_to_bytes(fig, "png", controls["export_scale"])
                img = Image.open(io.BytesIO(png_dat))
                buf = io.BytesIO()
                img.save(buf, format="TIFF")
                st.download_button("Download TIFF", buf.getvalue(), "plot.tiff", "image/tiff")
        except Exception as e:
            st.error(f"Export failed: {e}")

if __name__ == '__main__':
    main()  