import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from streamlit_option_menu import option_menu
import base64

# Page configuration
st.set_page_config(
    page_title="Enzyme Inhibitors in Drug Development",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS for beautiful styling
def local_css():
    st.markdown("""
    <style>
    .main-header {
        font-size: 3.5rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
        font-weight: 700;
        background: linear-gradient(45deg, #2E86AB, #A23B72);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    .section-header {
        font-size: 2rem;
        color: #2E86AB;
        border-left: 5px solid #A23B72;
        padding-left: 1rem;
        margin: 2rem 0 1rem 0;
    }
    .info-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 15px;
        color: white;
        margin: 1rem 0;
    }
    .mechanism-card {
        background: white;
        padding: 1.5rem;
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        margin: 1rem 0;
        border-left: 4px solid #A23B72;
    }
    .stButton>button {
        background: linear-gradient(45deg, #2E86AB, #A23B72);
        color: white;
        border: none;
        padding: 0.5rem 2rem;
        border-radius: 25px;
        font-weight: 600;
    }
    .drug-card {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 10px;
        border: 1px solid #e9ecef;
        margin: 0.5rem 0;
    }
    </style>
    """, unsafe_allow_html=True)

local_css()

# Header with navigation
def create_header():
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown('<div class="main-header">🧬 Enzyme Inhibitors in Drug Development</div>', unsafe_allow_html=True)
        
        # Navigation
        selected = option_menu(
            menu_title=None,
            options=["Overview", "Mechanisms", "Case Studies", "Kinetics", "Calculator", "References"],
            icons=["house", "gear", "book", "graph-up", "calculator", "journal-text"],
            menu_icon="cast",
            default_index=0,
            orientation="horizontal",
            styles={
                "container": {"padding": "0!important", "background-color": "#f8f9fa"},
                "icon": {"color": "#2E86AB", "font-size": "18px"}, 
                "nav-link": {"font-size": "16px", "text-align": "center", "margin":"0px", "--hover-color": "#eee"},
                "nav-link-selected": {"background-color": "#2E86AB"},
            }
        )
    return selected

# Overview Section
def show_overview():
    st.markdown('<div class="section-header">📊 Executive Summary</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        <div class="info-card">
        <h3>🚀 Why Enzyme Inhibitors Matter</h3>
        Enzyme inhibitors represent one of the most successful strategies in modern drug development, 
        targeting specific enzymes involved in disease pathways with remarkable precision and efficacy.
        </div>
        """, unsafe_allow_html=True)
        
        st.write("""
        ### Key Advantages:
        - **High Specificity**: Target only disease-related enzymes
        - **Predictable Kinetics**: Well-understood mechanisms
        - **Proven Success**: Multiple blockbuster drugs
        - **Versatile Applications**: Cancer, infections, metabolic diseases
        """)
    
    with col2:
        # Quick stats
        stats_data = {
            'Category': ['Approved Drugs', 'Clinical Trials', 'Market Value', 'Success Rate'],
            'Value': ['250+', '800+', '$150B+', '15%']
        }
        df_stats = pd.DataFrame(stats_data)
        fig = px.bar(df_stats, x='Value', y='Category', orientation='h',
                    color='Category', color_discrete_sequence=px.colors.qualitative.Set3)
        fig.update_layout(showlegend=False, height=300, margin=dict(l=0, r=0, t=0, b=0))
        st.plotly_chart(fig, width='stretch')

# Interactive Mechanisms Section
def show_mechanisms():
    st.markdown('<div class="section-header">🔬 Inhibition Mechanisms</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        mechanism = st.selectbox(
            "Select Inhibition Type:",
            ["Competitive Inhibition", "Non-competitive Inhibition", 
             "Uncompetitive Inhibition", "Mixed Inhibition"],
            key="mechanism_select"
        )
        
        st.markdown("""
        <div class="mechanism-card">
        <h4>🎯 How it works:</h4>
        """, unsafe_allow_html=True)
        
        if mechanism == "Competitive Inhibition":
            st.write("""
            - Inhibitor competes with substrate for active site
            - Resembles substrate structure
            - Effect can be overcome by increasing substrate concentration
            - **Example**: Statins (HMG-CoA reductase inhibitors)
            """)
        elif mechanism == "Non-competitive Inhibition":
            st.write("""
            - Binds to enzyme at site other than active site
            - Reduces enzyme efficiency without blocking substrate binding
            - Cannot be overcome by substrate concentration
            - **Example**: Heavy metal ions
            """)
        elif mechanism == "Uncompetitive Inhibition":
            st.write("""
            - Binds only to enzyme-substrate complex
            - Common in multi-substrate reactions
            - **Example**: Lithium for certain enzymes
            """)
        else:
            st.write("""
            - Combination of competitive and non-competitive features
            - Binds to both enzyme and enzyme-substrate complex
            - **Example**: Many kinase inhibitors
            """)
        st.markdown("</div>", unsafe_allow_html=True)
    
    with col2:
        # Interactive kinetics plot
        st.subheader("📈 Kinetic Effects")
        
        km = st.slider("Km value (affinity)", 0.1, 10.0, 1.0, 0.1, key="km_slider")
        vmax = st.slider("Vmax (maximum velocity)", 1, 100, 50, 1, key="vmax_slider")
        
        # Generate kinetic curves
        substrate = np.linspace(0.1, 20, 100)
        velocity_no_inhibitor = vmax * substrate / (km + substrate)
        
        if mechanism == "Competitive Inhibition":
            alpha = 2  # Inhibition factor
            velocity_inhibitor = vmax * substrate / (km * alpha + substrate)
        elif mechanism == "Non-competitive Inhibition":
            alpha = 2
            velocity_inhibitor = (vmax / alpha) * substrate / (km + substrate)
        else:
            velocity_inhibitor = vmax * substrate / (km * 1.5 + substrate * 1.2)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=substrate, y=velocity_no_inhibitor, 
                               name='No Inhibitor', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=substrate, y=velocity_inhibitor, 
                               name='With Inhibitor', line=dict(color='red', dash='dash')))
        
        fig.update_layout(
            title=f"Michaelis-Menten Kinetics - {mechanism}",
            xaxis_title="Substrate Concentration [S]",
            yaxis_title="Reaction Velocity v",
            height=400
        )
        st.plotly_chart(fig, width='stretch')

# IC50/Ki Calculator Section
def show_calculator():
    st.markdown('<div class="section-header">🧮 IC50 & Ki Calculator</div>', unsafe_allow_html=True)
    
    st.write("""Calculate inhibition constants and understand drug potency metrics.""")
    
    # Create tabs for different calculators
    tab1, tab2, tab3 = st.tabs(["IC50 Calculator", "Ki Calculator", "Dose-Response Curve"])
    
    with tab1:
        st.subheader("IC50 Calculator")
        st.write("**IC50** (Half maximal inhibitory concentration): The concentration of inhibitor required to reduce enzyme activity by 50%.")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("#### Input Parameters")
            
            # Activity measurements
            v0 = st.number_input("Initial Velocity (v₀)", min_value=0.0, value=100.0, step=1.0, 
                               help="Enzyme velocity without inhibitor")
            
            st.markdown("#### Inhibitor Concentrations & Activities")
            num_points = st.slider("Number of data points", 3, 10, 5)
            
            concentrations = []
            activities = []
            
            for i in range(num_points):
                col_a, col_b = st.columns(2)
                with col_a:
                    conc = st.number_input(f"[I]_{i+1} (µM)", min_value=0.0, value=float((i+1)*2), 
                                          step=0.1, key=f"conc_{i}")
                    concentrations.append(conc)
                with col_b:
                    act = st.number_input(f"Activity_{i+1} (%)", min_value=0.0, max_value=100.0, 
                                        value=float(max(10, 100 - i*18)), step=1.0, key=f"act_{i}")
                    activities.append(act)
        
        with col2:
            st.markdown("#### Results")
            
            if len(concentrations) >= 3:
                # Fit dose-response curve (Hill equation)
                # y = Bottom + (Top - Bottom) / (1 + (IC50/x)^HillSlope)
                # Simplified: assume Hill slope = 1
                
                # Calculate IC50 using interpolation
                conc_array = np.array(concentrations)
                act_array = np.array(activities)
                
                # Sort by concentration
                sorted_indices = np.argsort(conc_array)
                conc_sorted = conc_array[sorted_indices]
                act_sorted = act_array[sorted_indices]
                
                # Find IC50 (50% activity)
                if len(act_sorted) > 1 and act_sorted.max() > 50 and act_sorted.min() < 50:
                    ic50 = np.interp(50, act_sorted[::-1], conc_sorted[::-1])
                    
                    st.success(f"### IC50 = {ic50:.2f} µM")
                    
                    # Generate smooth curve
                    conc_smooth = np.linspace(0.01, max(concentrations)*1.2, 100)
                    # Hill equation with calculated IC50
                    hill_slope = 1.0
                    top = max(activities)
                    bottom = min(activities)
                    act_smooth = bottom + (top - bottom) / (1 + (conc_smooth/ic50)**hill_slope)
                    
                    # Plot
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=conc_sorted, y=act_sorted, mode='markers',
                                           name='Data', marker=dict(size=10, color='red')))
                    fig.add_trace(go.Scatter(x=conc_smooth, y=act_smooth, mode='lines',
                                           name='Fit', line=dict(color='blue')))
                    fig.add_hline(y=50, line_dash="dash", line_color="green", 
                                annotation_text="IC50")
                    
                    fig.update_layout(
                        title="Dose-Response Curve",
                        xaxis_title="Inhibitor Concentration (µM)",
                        yaxis_title="Activity (%)",
                        height=400,
                        xaxis_type="log"
                    )
                    st.plotly_chart(fig, width='stretch')
                    
                    # Interpretation
                    st.info(f"""**Interpretation:**
- At {ic50:.2f} µM, the enzyme activity is reduced to 50%
- Lower IC50 = More potent inhibitor
- Typical potency ranges:
  - Very potent: < 0.1 µM
  - Potent: 0.1-1 µM  
  - Moderate: 1-10 µM
  - Weak: > 10 µM
                    """)
                else:
                    st.warning("Not enough data points crossing 50% activity to calculate IC50")
    
    with tab2:
        st.subheader("Ki Calculator (Inhibition Constant)")
        st.write("**Ki**: Dissociation constant of the enzyme-inhibitor complex. Lower Ki = Stronger binding.")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("#### Input Parameters")
            
            inhibition_type = st.selectbox(
                "Inhibition Type",
                ["Competitive", "Non-competitive", "Uncompetitive"],
                key="ki_type"
            )
            
            ic50_input = st.number_input("IC50 (µM)", min_value=0.01, value=5.0, step=0.1,
                                        help="From IC50 assay")
            substrate_conc = st.number_input("[S] Substrate Concentration (µM)", 
                                           min_value=0.0, value=10.0, step=1.0)
            km_input = st.number_input("Km (µM)", min_value=0.01, value=5.0, step=0.1,
                                      help="Michaelis constant")
        
        with col2:
            st.markdown("#### Results")
            
            # Cheng-Prusoff equation for competitive inhibition:
            # Ki = IC50 / (1 + [S]/Km)
            
            if inhibition_type == "Competitive":
                ki = ic50_input / (1 + substrate_conc / km_input)
                st.success(f"### Ki = {ki:.3f} µM")
                
                st.markdown("**Calculation:**")
                st.latex(r"K_i = \frac{IC_{50}}{1 + \frac{[S]}{K_m}}")
                
                st.info(f"""**Cheng-Prusoff Equation (Competitive)**
- IC50 = {ic50_input} µM
- [S] = {substrate_conc} µM
- Km = {km_input} µM
- **Ki = {ki:.3f} µM**

Ki represents the true binding affinity of the inhibitor to the enzyme.
                """)
                
            elif inhibition_type == "Non-competitive":
                ki = ic50_input
                st.success(f"### Ki = {ki:.3f} µM")
                
                st.info(f"""**Non-competitive Inhibition**
- For non-competitive inhibitors: Ki ≈ IC50
- **Ki = {ki:.3f} µM**

The inhibitor binds to a site different from the active site.
                """)
            
            else:  # Uncompetitive
                ki = ic50_input * (1 + km_input / substrate_conc)
                st.success(f"### Ki = {ki:.3f} µM")
                
                st.latex(r"K_i = IC_{50} \times \left(1 + \frac{K_m}{[S]}\right)")
                
                st.info(f"""**Uncompetitive Inhibition**
- IC50 = {ic50_input} µM
- Km = {km_input} µM  
- [S] = {substrate_conc} µM
- **Ki = {ki:.3f} µM**

The inhibitor only binds to the enzyme-substrate complex.
                """)
    
    with tab3:
        st.subheader("Dose-Response Curve Generator")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("#### Curve Parameters")
            
            top_activity = st.slider("Top Activity (%)", 0, 100, 100, 1)
            bottom_activity = st.slider("Bottom Activity (%)", 0, 100, 0, 1)
            ic50_curve = st.number_input("IC50 (µM)", min_value=0.01, value=1.0, step=0.1, key="ic50_curve")
            hill_slope = st.slider("Hill Slope", 0.5, 4.0, 1.0, 0.1,
                                 help="Steepness of the curve")
            
            conc_range_max = st.number_input("Max Concentration (µM)", min_value=0.1, value=100.0, step=1.0)
        
        with col2:
            st.markdown("#### Generated Curve")
            
            # Generate dose-response curve
            concentrations_curve = np.logspace(-3, np.log10(conc_range_max), 100)
            
            # Hill equation: y = Bottom + (Top - Bottom) / (1 + (x/IC50)^HillSlope)
            response = bottom_activity + (top_activity - bottom_activity) / \
                      (1 + (concentrations_curve / ic50_curve)**hill_slope)
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=concentrations_curve, y=response, mode='lines',
                                   line=dict(color='purple', width=3)))
            fig.add_hline(y=50, line_dash="dash", line_color="green",
                        annotation_text=f"IC50 = {ic50_curve} µM")
            fig.add_vline(x=ic50_curve, line_dash="dash", line_color="green")
            
            fig.update_layout(
                title="Dose-Response Curve",
                xaxis_title="Inhibitor Concentration (µM)",
                yaxis_title="Activity (%)",
                height=400,
                xaxis_type="log"
            )
            st.plotly_chart(fig, width='stretch')
            
            st.markdown("**Hill Equation:**")
            st.latex(r"y = Bottom + \frac{Top - Bottom}{1 + \left(\frac{[I]}{IC_{50}}\right)^{h}}")
            st.write(f"where h = {hill_slope} (Hill slope)")

# References Section
def show_references():
    st.markdown('<div class="section-header">📚 References & Resources</div>', unsafe_allow_html=True)
    
    st.write("""This educational tool is based on established principles in biochemistry and pharmacology.""")
    
    # Create tabs for different reference categories
    tab1, tab2, tab3, tab4 = st.tabs(["Key Papers", "Textbooks", "Online Resources", "Citation"])
    
    with tab1:
        st.subheader("📄 Key Research Papers")
        
        st.markdown("""#### Enzyme Inhibition Mechanisms
1. **Copeland, R.A.** (2013). "Evaluation of Enzyme Inhibitors in Drug Discovery." 
   *Methods of Biochemical Analysis*, 46, 1-265. 
   [DOI: 10.1002/9781118540398](https://doi.org/10.1002/9781118540398)

2. **Cheng, Y., & Prusoff, W.H.** (1973). "Relationship between the inhibition constant (Ki) and the concentration of inhibitor which causes 50 per cent inhibition (IC50) of an enzymatic reaction." 
   *Biochemical Pharmacology*, 22(23), 3099-3108.
   [DOI: 10.1016/0006-2952(73)90196-2](https://doi.org/10.1016/0006-2952(73)90196-2)

3. **Segel, I.H.** (1975). "Enzyme Kinetics: Behavior and Analysis of Rapid Equilibrium and Steady-State Enzyme Systems." 
   *Wiley-Interscience*, New York.

#### Drug Development Case Studies

4. **Istvan, E.S., & Deisenhofer, J.** (2001). "Structural mechanism for statin inhibition of HMG-CoA reductase." 
   *Science*, 292(5519), 1160-1164.
   [DOI: 10.1126/science.1059344](https://doi.org/10.1126/science.1059344)

5. **Kohl, N.E., et al.** (1988). "Active human immunodeficiency virus protease is required for viral infectivity." 
   *PNAS*, 85(13), 4686-4690.
   [DOI: 10.1073/pnas.85.13.4686](https://doi.org/10.1073/pnas.85.13.4686)

6. **Druker, B.J., et al.** (2001). "Efficacy and safety of a specific inhibitor of the BCR-ABL tyrosine kinase in chronic myeloid leukemia." 
   *New England Journal of Medicine*, 344(14), 1031-1037.
   [DOI: 10.1056/NEJM200104053441401](https://doi.org/10.1056/NEJM200104053441401)
        """)
    
    with tab2:
        st.subheader("📖 Recommended Textbooks")
        
        st.markdown("""1. **Berg, J.M., Tymoczko, J.L., & Stryer, L.** (2019). 
   *Biochemistry* (9th ed.). W.H. Freeman and Company.
   - Chapter 8: Enzymes: Basic Concepts and Kinetics
   - Chapter 9: Catalytic Strategies

2. **Nelson, D.L., & Cox, M.M.** (2017). 
   *Lehninger Principles of Biochemistry* (7th ed.). W.H. Freeman.
   - Chapter 6: Enzymes

3. **Silverman, R.B., & Holladay, M.W.** (2014). 
   *The Organic Chemistry of Drug Design and Drug Action* (3rd ed.). Academic Press.
   - Chapter 4: Enzyme Inhibition and Inactivation

4. **Copeland, R.A.** (2005). 
   *Evaluation of Enzyme Inhibitors in Drug Discovery: A Guide for Medicinal Chemists and Pharmacologists*. 
   Wiley-Interscience.

5. **Wlodawer, A., & Vondrasek, J.** (1998). 
   *Annual Review of Biophysics and Biomolecular Structure*, 27, 249-284.
   - "Inhibitors of HIV-1 protease: a major success of structure-assisted drug design"
        """)
    
    with tab3:
        st.subheader("🌐 Online Resources")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""#### Databases & Tools
- **[PubChem](https://pubchem.ncbi.nlm.nih.gov/)** - Chemical information database
- **[DrugBank](https://go.drugbank.com/)** - Drug and drug target database  
- **[BRENDA](https://www.brenda-enzymes.org/)** - Enzyme information system
- **[ChEMBL](https://www.ebi.ac.uk/chembl/)** - Bioactive molecules database
- **[PDB](https://www.rcsb.org/)** - Protein Data Bank

#### Educational Resources
- **[Khan Academy - Enzyme Kinetics](https://www.khanacademy.org/science/biology/energy-and-enzymes)**
- **[NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/)** - Free biochemistry textbooks
            """)
        
        with col2:
            st.markdown("""#### Organizations
- **[FDA - Drug Development](https://www.fda.gov/drugs/development-approval-process-drugs)**
- **[WHO Essential Medicines](https://www.who.int/groups/expert-committee-on-selection-and-use-of-essential-medicines)**
- **[ACS Chemical Biology](https://pubs.acs.org/journal/acbcct)**

#### Journals
- *Nature Reviews Drug Discovery*
- *Journal of Medicinal Chemistry*
- *Biochemistry*
- *Drug Discovery Today*
            """)
    
    with tab4:
        st.subheader("📝 How to Cite This Tool")
        
        st.markdown("""If you use this interactive educational tool in your work, please cite as:
        """)
        
        st.code("""Enzyme Inhibitors in Drug Development: An Interactive Educational Tool
(2025). Created with Streamlit.
Available at: [Your Deployment URL]
Accessed: [Date]""", language="text")
        
        st.markdown("**BibTeX Format:**")
        st.code("""@misc{enzyme_inhibitors_tool_2025,
  title={Enzyme Inhibitors in Drug Development: An Interactive Educational Tool},
  author={[Your Name]},
  year={2025},
  howpublished={\\url{[Your Deployment URL]}},
  note={Interactive web application built with Streamlit}
}""", language="bibtex")
        
        st.markdown("---")
        st.info("""**Disclaimer:** This tool is designed for educational purposes only. 
All data and calculations should be verified with primary literature sources 
and are not intended for clinical or commercial drug development decisions.""")

# Case Studies Section
def show_case_studies():
    st.markdown('<div class="section-header">💊 Successful Drug Case Studies</div>', unsafe_allow_html=True)
    
    case_study = st.radio(
        "Select Drug Case Study:",
        ["Statins (Cholesterol)", "HIV Protease Inhibitors", "ACE Inhibitors (Blood Pressure)", 
         "Kinase Inhibitors (Cancer)", "COX-2 Inhibitors (Pain)"],
        horizontal=True
    )
    
    if case_study == "Statins (Cholesterol)":
        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("""
            <div class="drug-card">
            <h4>💊 Atorvastatin (Lipitor)</h4>
            <b>Target:</b> HMG-CoA reductase<br>
            <b>Mechanism:</b> Competitive inhibition<br>
            <b>Disease:</b> Hypercholesterolemia<br>
            <b>Impact:</b> One of best-selling drugs ever
            </div>
            """, unsafe_allow_html=True)
        with col2:
            # Simple efficacy chart
            data = {'Year': [1990, 2000, 2010, 2020],
                   'Cardiac Deaths per 1000': [8.5, 6.2, 4.1, 2.8]}
            df = pd.DataFrame(data)
            fig = px.line(df, x='Year', y='Cardiac Deaths per 1000', 
                         title="Impact on Cardiac Mortality")
            st.plotly_chart(fig, width='stretch')

# Main application flow
def main():
    selected_section = create_header()
    
    if selected_section == "Overview":
        show_overview()
    elif selected_section == "Mechanisms":
        show_mechanisms()
    elif selected_section == "Case Studies":
        show_case_studies()
    elif selected_section == "Kinetics":
        st.markdown('<div class="section-header">📊 Interactive Kinetics Explorer</div>', unsafe_allow_html=True)
        # Add kinetic explorer content
    elif selected_section == "Calculator":
        show_calculator()
    elif selected_section == "References":
        show_references()
    else:
        st.markdown('<div class="section-header">🚀 Future Directions</div>', unsafe_allow_html=True)
        # Add future trends content
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #666;'>
    <i>Interactive Educational Poster - Enzyme Inhibitors in Drug Development</i><br>
    Created with Streamlit | For educational purposes
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
