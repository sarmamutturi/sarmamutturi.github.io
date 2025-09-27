import streamlit as st
import requests
import json
import xml.etree.ElementTree as ET
import pandas as pd
import time
import random
import io
from typing import List, Dict, Any

# --- EXACT MONOISOTOPIC MASSES (M) ---
M_C = 12.000000000  # Exact mass of C-12
M_H = 1.007825032   # Exact mass of H-1
M_O = 15.994914619  # Exact mass of O-16

# Define Adduct Masses
ADDUCT_MASSES = {
    '[M+H]+ (Proton)': 1.007276,
    '[M+Na]+ (Sodium)': 22.989769,
    '[M+K]+ (Potassium)': 38.963704,
}

# --- Common Name Lookup Helper ---
def get_common_name(cid, formula):
    """Provides common names for key CIDs or formulas, focusing on saccharides."""
    common_names = {
        '5460037': 'Isomaltotriose',
        '439242': 'Maltotriose',
    }
    
    if str(cid) in common_names:
        return common_names[str(cid)]
        
    c_count = 0
    import re
    match = re.search(r'C(\d+)H(\d+)O(\d+)', formula)
    if match:
        c_count = int(match.group(1))

    if c_count == 6 and formula == 'C6H12O6':
        return 'Hexose (e.g., Glucose, Fructose)'
    elif c_count == 12 and formula == 'C12H22O11':
        return 'Disaccharide (e.g., Sucrose, Maltose)'
    elif c_count == 18:
        return 'Trisaccharide (e.g., Maltotriose, Isomaltotriose)'
    elif c_count == 24:
        return 'Tetrasaccharide'
    elif c_count >= 30:
        return 'Oligosaccharide / Larger Compound'
        
    return 'N/A'

# --- Core Search Functions ---
def generate_formulas_for_mass(neutral_mass, tolerance):
    """Generates all C, H, O elemental formulas matching the neutral mass."""
    
    C_min, C_max = 5, 50  
    O_min, O_max = 5, 30
    H_max = int(2 * C_max + 2 - O_min)  
    
    valid_formulas = {}
    
    for c in range(C_min, C_max + 1):
        for o in range(O_min, O_max + 1):
            
            M_c_o = c * M_C + o * M_O
            M_rem = neutral_mass - M_c_o
            
            h_min_calc = int((M_rem - tolerance) / M_H)
            h_max_calc = int((M_rem + tolerance) / M_H)
            
            h_min = max(0, h_min_calc)
            h_max = min(H_max, h_max_calc)
            
            for h in range(h_min, h_max + 1):
                if h % 2 != 0: 
                     continue
                if h > 2 * c + 2:
                    continue

                exact_mass = c * M_C + h * M_H + o * M_O
                if abs(exact_mass - neutral_mass) <= tolerance:
                    formula = f"C{c}H{h}O{o}"
                    if formula not in valid_formulas or abs(exact_mass - neutral_mass) < abs(valid_formulas[formula] - neutral_mass):
                         valid_formulas[formula] = exact_mass
                         
    return list(valid_formulas.keys())


def search_pubchem_by_formulas(formulas, mz_adduct, neutral_mass, adduct_label, progress_bar, formula_count):
    """Searches PubChem for all generated molecular formulas."""
    all_results = []
    
    for i, formula in enumerate(formulas):
        progress = (i + 1) / formula_count
        progress_bar.progress(progress)
        
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={formula}[MolecularFormula] AND C AND H AND O&retmax=50"
        
        try:
            response = requests.get(search_url, timeout=30)
            response.raise_for_status()
            root = ET.fromstring(response.text)
            cids = [cid.text for cid in root.findall(".//Id")]
            
            if not cids:
                continue

            cid_list_str = ",".join(cids)
            properties_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_list_str}/property/MolecularFormula,CanonicalSMILES,MolecularWeight,IUPACName/JSON"
            
            time.sleep(1 + random.uniform(0, 0.5)) # Be polite to API
            prop_response = requests.get(properties_url, timeout=30)
            prop_response.raise_for_status()
            prop_data = prop_response.json()
            
            for props in prop_data.get('PropertyTable', {}).get('Properties', []):
                 cid_str = str(props.get('CID'))
                 formula_str = props.get('MolecularFormula', 'N/A')
                 
                 all_results.append({
                    "Observed m/z": f"{mz_adduct:.6f}",
                    "Calculated Neutral Mass": f"{neutral_mass:.4f}",
                    "Molecular Formula": formula_str,
                    "Exact Mass": props.get('MolecularWeight', 'N/A'),
                    "PubChem CID": cid_str,
                    "IUPAC Name": props.get('IUPACName', 'N/A'),
                    "Common Name": get_common_name(cid_str, formula_str),
                    "Adduct": adduct_label,
                    "Source Database": "Formula Prediction Layer"
                })
        except Exception:
            continue
            
    return all_results

# --- STREAMLIT DRIVER FUNCTION ---
def run_mz_formula_prediction_batch_st(mz_values, adduct_mass, adduct_label, tolerance, status_placeholder):
    
    final_combined_results = []
    total_mz = len(mz_values)
    
    for idx, mz_value in enumerate(mz_values):
        neutral_mass = mz_value - adduct_mass
        
        status_placeholder.info(f"🔍 **Processing {adduct_label}** for m/z: **{mz_value:.6f}** ({idx+1}/{total_mz})")

        # STAGE 1: FORMULA PREDICTION
        with st.spinner(f"Stage 1: Generating formulas for neutral mass {neutral_mass:.4f}..."):
            potential_formulas = generate_formulas_for_mass(neutral_mass, tolerance)

        if not potential_formulas:
            status_placeholder.warning(f"❌ No plausible C, H, O formulas found for m/z {mz_value:.6f}.")
            continue

        # STAGE 2: SEARCH PUBCHEM
        formula_count = len(potential_formulas)
        status_placeholder.info(f"Stage 2: Searching PubChem with {formula_count} potential formulas...")
        
        progress_bar = st.progress(0.0)
        
        all_results = search_pubchem_by_formulas(
            potential_formulas, mz_value, neutral_mass, adduct_label, progress_bar, formula_count
        )
        progress_bar.empty() # Clear the progress bar after completion
        
        if all_results:
            df_current = pd.DataFrame(all_results).drop_duplicates(subset=["PubChem CID", "Molecular Formula"])
            final_combined_results.append(df_current)
            status_placeholder.success(f"✅ Found {len(df_current)} unique compounds for m/z {mz_value:.6f}.")
        else:
            status_placeholder.info(f"--- No compounds found for m/z {mz_value:.6f}. ---")
            
    # --- FINAL OUTPUT CONSOLIDATION ---
    if final_combined_results:
        df_final = pd.concat(final_combined_results, ignore_index=True)
        df_final = df_final.sort_values(by="Observed m/z")
        
        cols = [
            "Observed m/z", "Calculated Neutral Mass", "Molecular Formula", "Exact Mass", 
            "Common Name", "IUPAC Name", "PubChem CID", "Adduct", "Source Database"
        ]
        
        df_final = df_final[cols] 
        
        return df_final
    else:
        return pd.DataFrame()


# ====================================================================
# STREAMLIT UI LAYOUT
# ====================================================================

st.set_page_config(
    page_title="LC-MS Formula Prediction Tool",
    layout="wide"
)

st.title("🧪 LC-MS Formula & Compound Prediction Tool")
st.markdown("Enter m/z values to generate possible $\\text{C}_x\\text{H}_y\\text{O}_z$ formulas and search PubChem.")

# --- SIDEBAR FOR INPUTS ---
with st.sidebar:
    st.header("Parameters")
    
    # 1. Adduct Selection
    adduct_label = st.selectbox(
        "Select Adduct Type",
        options=list(ADDUCT_MASSES.keys()),
        index=1, # Default to [M+Na]+
        help="Select the ion used for ionization (e.g., [M+Na]+ is common for saccharides)."
    )
    adduct_mass = ADDUCT_MASSES[adduct_label]
    st.info(f"Adduct Mass: **{adduct_mass:.6f} u**")

    # 2. Tolerance Input
    tolerance = st.number_input(
        "Mass Tolerance (Da)",
        min_value=0.001,
        max_value=0.100,
        value=0.050,  # Default 50 mDa
        step=0.005,
        format="%.4f",
        help="Maximum mass deviation allowed for formula calculation (e.g., 0.05 Da)."
    )

    # 3. m/z Input
    default_mz_values = "527.16, 689.20, 365.12, 851.24"
    mz_input = st.text_area(
        "List of m/z values (comma-separated)",
        value=default_mz_values,
        height=150,
        help="Enter one or more m/z values, separated by commas or new lines."
    )

    # 4. RUN Button
    run_button = st.button("🚀 Run Analysis", type="primary")
    
    # --- ATTRIBUTION / CONTACT INFORMATION ---
    st.markdown("---")
    st.markdown(
        """
        **Developed by Sarma Mutturi Lab**
        * **Contact:** <smutturi.cftri@csir.res.in>
        * **Site:** [https://sarmamutturi.github.io/](https://sarmamutturi.github.io/)
        """
    )
    st.markdown("---")

# --- MAIN AREA FOR RESULTS ---
status_placeholder = st.empty()
status_placeholder.info("Enter m/z values and click 'Run Analysis' to start.")

if run_button:
    # ... (Rest of the execution logic)
    mz_values_clean = []
    
    # Clean and parse m/z values
    try:
        # Replace commas/newlines with spaces, split, and convert to float
        mz_str_list = mz_input.replace(',', ' ').replace('\n', ' ').split()
        for mz_str in mz_str_list:
            mz_values_clean.append(float(mz_str.strip()))
            
    except ValueError:
        status_placeholder.error("Invalid input detected in m/z list. Please enter valid numbers only.")
        st.stop()
        
    if not mz_values_clean:
        status_placeholder.error("Please enter at least one m/z value to analyze.")
        st.stop()

    status_placeholder.empty() # Clear initial status
    st.subheader(f"Results for {len(mz_values_clean)} m/z values ({adduct_label})")
    
    # Run the core logic
    df_results = run_mz_formula_prediction_batch_st(
        mz_values=mz_values_clean,
        adduct_mass=adduct_mass,
        adduct_label=adduct_label,
        tolerance=tolerance,
        status_placeholder=status_placeholder
    )

    if not df_results.empty:
        # 1. Display the DataFrame
        st.dataframe(df_results, use_container_width=True)
        
        # 2. Download Button
        @st.cache_data
        def convert_df_to_excel(df):
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='Compound Results')
            processed_data = output.getvalue()
            return processed_data

        excel_data = convert_df_to_excel(df_results)
        
        st.download_button(
            label="💾 Download Results as Excel (.xlsx)",
            data=excel_data,
            file_name="lcms_batch_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
        status_placeholder.success(f"✅ Analysis Complete! Total unique compounds found: **{len(df_results)}**.")
    else:
        status_placeholder.error("No compounds were found matching the criteria for the given m/z values.")
