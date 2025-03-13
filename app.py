import streamlit as st
import math
import numpy as np
import csv
import matplotlib.pyplot as plt
import io
import pandas as pd

# --------------------------
# Function Definitions
# --------------------------
def calculate_ft_mic(dose, fu, vd, mic, clt, di):
    try:
        term1 = math.log((dose * fu) / (vd * mic))
        term2 = vd / clt
        term3 = 100 / di
        ft_mic = term1 * term2 * term3
        return max(ft_mic, 0)
    except (ValueError, ZeroDivisionError):
        return "Invalid parameters. Please check your inputs."

def calculate_lognormal_params(mean, std):
    variance = std ** 2
    mu = math.log(mean ** 2 / math.sqrt(variance + mean ** 2))
    sigma = math.sqrt(math.log(1 + variance / (mean ** 2)))
    return mu, sigma

def monte_carlo_simulation(dose_intervals, fu_mean, fu_std, vd_mean, vd_std, clt_mean, clt_std, mics, num_patients, target):
    fu_mu, fu_sigma = calculate_lognormal_params(fu_mean, fu_std)
    vd_mu, vd_sigma = calculate_lognormal_params(vd_mean, vd_std)
    clt_mu, clt_sigma = calculate_lognormal_params(clt_mean, clt_std)

    results = {}
    for mic in mics:
        results[mic] = []
        for dose, di in dose_intervals:
            success_count = 0
            for _ in range(num_patients):
                fu = max(np.random.lognormal(mean=fu_mu, sigma=fu_sigma), 0)
                vd = max(np.random.lognormal(mean=vd_mu, sigma=vd_sigma), 0)
                clt = max(np.random.lognormal(mean=clt_mu, sigma=clt_sigma), 0)
                ft_mic = calculate_ft_mic(dose, fu, vd, mic, clt, di)
                if ft_mic >= target:
                    success_count += 1
            pta = (success_count / num_patients) * 100
            results[mic].append((dose, di, pta))
    return results

def calculate_cfr(results, mic_distribution):
    cfr_results = {}
    for mic, combinations in results.items():
        for dose, di, pta in combinations:
            if (dose, di) not in cfr_results:
                cfr_results[(dose, di)] = 0
            cfr_results[(dose, di)] += pta * mic_distribution.get(mic, 0)
    return cfr_results

# --------------------------
# Streamlit App Code
# --------------------------
def main():
    # Header and description are always at the top
    st.title("Pharmacokinetic PTA and CFR Simulator")
    st.write("This app calculates the Probability of Target Attainment (PTA) and Cumulative Fraction Response (CFR) for various dosing regimens.")

    # --------------------------
    # Sidebar: Common Input Parameters
    # --------------------------
    st.sidebar.header("Input Parameters")
    
    # Pharmacokinetic parameters
    fu_mean = st.sidebar.number_input("Fraction Unbound (fu) Mean", value=0.073, format="%.3f")
    fu_std = st.sidebar.number_input("Fraction Unbound (fu) Std", value=0.032, format="%.3f")
    vd_mean = st.sidebar.number_input("Volume of Distribution (L) Mean", value=7.8, format="%.1f")
    vd_std = st.sidebar.number_input("Volume of Distribution (L) Std", value=5.4, format="%.1f")
    clt_mean = st.sidebar.number_input("Total Clearance (L/h) Mean", value=0.83, format="%.2f")
    clt_std = st.sidebar.number_input("Total Clearance (L/h) Std", value=0.83, format="%.2f")
    
    # MIC inputs
    mics_input = st.sidebar.text_input("MIC Values (comma separated)", "0.125,0.25,0.5,1.0,2.0,4.0,8.0,16.0,32.0")
    target = st.sidebar.number_input("PK/PD Target (%fT > MIC)", value=55, step=1)
    num_patients = st.sidebar.number_input("Number of Simulated Patients", value=5000, step=100)
    
    # MIC distribution input
    mic_distribution_input = st.sidebar.text_area(
        "MIC Distribution (format: mic:fraction, separated by commas)",
        value="0.125:0.0,0.25:0.0,0.5:0.0,1.0:0.011,2.0:0.338,4.0:0.599,8.0:0.023,16.0:0.026,32.0:0.003"
    )
    
    # Parse MIC values
    try:
        mic_list = [float(x.strip()) for x in mics_input.split(",")]
    except Exception:
        st.error("Please enter valid MIC values separated by commas.")
        return
    
    # Parse MIC distribution
    mic_distribution = {}
    try:
        pairs = mic_distribution_input.split(",")
        for pair in pairs:
            mic_key, fraction = pair.split(":")
            mic_distribution[float(mic_key.strip())] = float(fraction.strip())
    except Exception:
        st.error("Please enter valid MIC distribution values in the format mic:fraction separated by commas.")
        return

    # --------------------------
    # Sidebar: Custom Regimen Input
    # --------------------------
    # Initialize session state for regimen combinations if not already present
    if "dose_intervals" not in st.session_state:
        st.session_state.dose_intervals = []
    
    st.sidebar.subheader("Add Regimen Combination")
    dose_input = st.sidebar.number_input("Dose (mg)", value=2000, step=100)
    # Interval must be between 1 and 24 hours
    interval_input = st.sidebar.number_input("Interval (h)", value=24, min_value=1, max_value=24, step=1)
    
    if st.sidebar.button("Add Regimen"):
        regimen = (dose_input, interval_input)
        if regimen not in st.session_state.dose_intervals:
            st.session_state.dose_intervals.append(regimen)
    
    # Display current regimen combinations as a neat bullet list
    if st.session_state.dose_intervals:
        regimen_list = "\n".join(
            ["- Dose: {} mg, Interval: {} h".format(dose, di) for dose, di in st.session_state.dose_intervals]
        )
        st.sidebar.markdown("**Current Regimen Combinations:**\n" + regimen_list)
    else:
        st.sidebar.info("No regimen combinations added yet.")
    
    # Ensure at least one regimen exists before simulation
    if not st.session_state.dose_intervals:
        st.info("Please add at least one regimen combination from the sidebar.")
        return

    # --------------------------
    # Run Simulation
    # --------------------------
    if st.sidebar.button("Run Simulation"):
        with st.spinner("Simulating..."):
            results = monte_carlo_simulation(
                st.session_state.dose_intervals, fu_mean, fu_std, vd_mean, vd_std,
                clt_mean, clt_std, mic_list, int(num_patients), target
            )
            cfr_results = calculate_cfr(results, mic_distribution)
        st.success("Simulation complete!")
        
        # --------------------------
        # Display PTA Results Table (MICs as rows, regimens as columns)
        # --------------------------
        mic_list_sorted = sorted(mic_list)
        pta_matrix = {}
        for mic in mic_list_sorted:
            row = {}
            for dose, di in st.session_state.dose_intervals:
                regimen = f"Dose: {dose} mg, Interval: {di} h"
                pta_val = next((pta for d, i, pta in results[mic] if d == dose and i == di), None)
                row[regimen] = round(pta_val, 2) if pta_val is not None else None
            pta_matrix[mic] = row
        df_pta_matrix = pd.DataFrame.from_dict(pta_matrix, orient="index")
        df_pta_matrix.index.name = "MIC (mg/L)"
        
        st.header("PTA Results")
        st.dataframe(df_pta_matrix)
        
        # --------------------------
        # Display CFR Results Table
        # --------------------------
        cfr_data = [
            {"Regimen": f"Dose: {dose} mg, Interval: {di} h", "CFR (%)": round(cfr, 2)}
            for (dose, di), cfr in cfr_results.items()
        ]
        df_cfr = pd.DataFrame(cfr_data)
        
        st.header("CFR Results")
        html_table = df_cfr.to_html(index=False)
        st.markdown(html_table, unsafe_allow_html=True)
        
        # --------------------------
        # Display PTA vs MIC Plot
        # --------------------------
        st.header("PTA vs MIC Plot")
        fig1 = plt.figure(figsize=(10, 6))
        for dose, di in set((dose, di) for combos in results.values() for dose, di, _ in combos):
            y_values = []
            for mic in mic_list_sorted:
                for d, interval, pta in results.get(mic, []):
                    if d == dose and interval == di:
                        y_values.append(pta)
            if len(y_values) == len(mic_list_sorted):
                plt.plot(mic_list_sorted, y_values, label=f"Dose: {dose} mg, Interval: {di} h")
        plt.xscale("log")
        plt.xticks(mic_list_sorted, [str(mic) for mic in mic_list_sorted])
        plt.gca().set_xticks(mic_list_sorted)
        plt.gca().set_xticklabels([str(mic) for mic in mic_list_sorted])
        plt.minorticks_off()
        plt.xlabel("MIC Values (mg/L)")
        plt.ylabel("Percentage Achieved Target (%)")
        plt.title("Percentage of Subjects Achieving Target %fT>MIC")
        plt.axhline(y=95, color='red', linestyle='--', linewidth=1.5, label="95% Target")
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.legend()
        plt.tight_layout()
        st.pyplot(fig1)
        
        # --------------------------
        # Display CFR Bar Chart
        # --------------------------
        st.header("CFR Bar Chart")
        fig2 = plt.figure(figsize=(12, 6))
        regimens = [f"Dose: {dose} mg, Interval: {di} h" for dose, di in cfr_results.keys()]
        cfr_values = list(cfr_results.values())
        plt.bar(regimens, cfr_values)
        plt.xlabel("Regimens")
        plt.ylabel("CFR (%)")
        plt.title("Cumulative Fraction Response (CFR) by Regimen")
        plt.xticks(rotation=45, ha="right")
        plt.ylim(0, 100)
        plt.axhline(y=95, color='red', linestyle='--', linewidth=1.5, label="95% Target")
        plt.tight_layout()
        st.pyplot(fig2)
        
        # --------------------------
        # Download Results as CSV
        # --------------------------
        st.header("Download Results")
        pta_csv = io.StringIO()
        writer = csv.writer(pta_csv)
        writer.writerow(["MIC (mg/L)"] + [f"Dose: {dose} mg, Interval: {di} h" for dose, di in st.session_state.dose_intervals])
        for mic in mic_list_sorted:
            row = [mic]
            for dose, di in st.session_state.dose_intervals:
                pta_val = next((pta for d, i, pta in results[mic] if d == dose and i == di), 0)
                row.append(pta_val)
            writer.writerow(row)
        st.download_button("Download PTA CSV", data=pta_csv.getvalue(), file_name="pta_results_matrix.csv", mime="text/csv")
        
        cfr_csv = io.StringIO()
        writer = csv.writer(cfr_csv)
        writer.writerow(["Dose (mg)", "Interval (h)", "CFR (%)"])
        for (dose, di), cfr in cfr_results.items():
            writer.writerow([dose, di, cfr])
        st.download_button("Download CFR CSV", data=cfr_csv.getvalue(), file_name="cfr_results.csv", mime="text/csv")

if __name__ == "__main__":
    main()
