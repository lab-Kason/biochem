import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from streamlit_option_menu import option_menu

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
        st.markdown('<p style="text-align: center; font-size: 1.1rem; color: #666; margin-bottom: 1rem;">An Interactive Educational Platform for Understanding Drug Development</p>', unsafe_allow_html=True)
        
        # Complete Application Guide
        with st.expander("🎯 **START HERE: Complete User Guide**", expanded=False):
            st.markdown("""
            ### Welcome! 👋
            
            This interactive platform helps you understand enzyme inhibitors in drug development.
            Perfect for students, researchers, and anyone interested in biochemistry and pharmacology.
            
            ---
            
            ### 📋 Six Interactive Sections:
            
            **1. 🏠 Overview**
            - Statistics on enzyme inhibitor drugs
            - Market impact and therapeutic applications
            - *Best for: Understanding the big picture*
            
            **2. ⚙️ Mechanisms**
            - Interactive simulator for 4 inhibition types
            - Real-time kinetic curve visualization
            - Adjustable parameters (Km, Vmax, inhibitor strength)
            - *Best for: Learning how inhibitors work*
            
            **3. 📚 Case Studies**
            - 5 blockbuster drugs (Statins, HIV protease, ACE, Kinase, COX-2)
            - Clinical trial data with interactive charts
            - *Best for: Real-world applications*
            
            **4. 📊 Kinetics**
            - Lineweaver-Burk plots
            - Eadie-Hofstee transformation
            - Michaelis-Menten simulator
            - *Best for: Advanced kinetics analysis*
            
            **5. 🧮 Calculator**
            - IC50 calculator (with CSV export)
            - Ki calculator (Cheng-Prusoff equations)
            - Dose-response curve generator
            - *Best for: Analyzing your experimental data*
            
            **6. 📖 References**
            - 26 peer-reviewed papers
            - 12 textbooks
            - Online databases and resources
            - *Best for: Citations and further reading*
            
            ---
            
            ### 💡 Quick Start Tips:
            
            - **Beginners:** Start with Overview → Mechanisms → Case Studies
            - **Students with data:** Go to Calculator to analyze results
            - **In-depth learners:** Work through all sections sequentially
            - **Each section has its own guide** - Look for "📖 How to Use" expandable boxes
            
            **All tools are interactive - adjust sliders to see immediate results!**
            
            *Select a section from the menu below to begin ⬇️*
            """)
        
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
    
    # User Guide
    with st.expander("📖 How to Use This Section", expanded=False):
        st.markdown("""
        **Overview Section Guide:**
        
        This section provides key statistics about enzyme inhibitors in modern medicine.
        
        **What you'll see:**
        - 📊 **Success statistics** - How many FDA drugs target enzymes
        - 🎯 **Market data** - Financial impact of enzyme inhibitor drugs
        - 📈 **Interactive bar chart** - Visual breakdown of therapeutic applications
        
        **Simply scroll down to explore the statistics and visualizations.**
        
        All data is referenced from peer-reviewed publications (see References section).
        """)
    
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
        
        st.caption("""*Data sources: FDA Drug Approvals Database (2024); ClinicalTrials.gov; 
        Evaluate Pharma Market Reports (2023). Success rate from DiMasi et al. (2016). 
        See References section for full citations.*""")

# Interactive Mechanisms Section
def show_mechanisms():
    st.markdown('<div class="section-header">🔬 Inhibition Mechanisms</div>', unsafe_allow_html=True)
    
    # User Guide for Mechanisms
    with st.expander("📖 How to Use This Interactive Simulator", expanded=False):
        st.markdown("""
        **Interactive Enzyme Inhibition Simulator**
        
        **Step-by-Step Guide:**
        1. **Select inhibition type** from dropdown (left panel)
        2. **Read the mechanism description** to understand how it works
        3. **Adjust Km slider** (0.1-10.0) - substrate binding affinity
        4. **Adjust Vmax slider** (1-100) - maximum reaction velocity
        5. **Toggle "Show Inhibitor Effect"** checkbox to compare curves
        6. **Use Inhibitor Strength (α) slider** (1.0-5.0) to see dose-dependent effects
        
        **Understanding the Parameters:**
        - **Km:** Lower values = enzyme has higher affinity for substrate
        - **Vmax:** Higher values = enzyme can work faster
        - **α (alpha):** Inhibitor strength factor (higher = stronger inhibition)
        
        **Interpreting the Plot:**
        - **Blue curve:** Normal enzyme activity (no inhibitor)
        - **Red dashed curve:** Activity with inhibitor present
        - Watch how different mechanisms affect curve shape!
        
        **What to Observe:**
        - **Competitive:** Can't reach same Vmax at low [S], but can at high [S]
        - **Non-competitive:** Vmax reduced, but Km appears unchanged
        - **Uncompetitive:** Both Vmax and apparent Km reduced
        - **Mixed:** Combination of competitive and non-competitive effects
        
        *Try different combinations to see how each parameter affects enzyme kinetics!*
        """)
    
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
        
        # Add inhibitor strength control
        show_inhibitor = st.checkbox("Show Inhibitor Effect", value=True, key="show_inh_mech")
        if show_inhibitor:
            inhibitor_strength = st.slider("Inhibitor Strength (alpha)", 1.0, 5.0, 2.0, 0.1, 
                                          key="alpha_slider",
                                          help="Higher values = stronger inhibition")
        
        # Generate kinetic curves
        substrate = np.linspace(0.1, 20, 100)
        velocity_no_inhibitor = vmax * substrate / (km + substrate)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=substrate, y=velocity_no_inhibitor, 
                               name='No Inhibitor', line=dict(color='blue')))
        
        if show_inhibitor:
            alpha = inhibitor_strength  # Use slider value
            
            if mechanism == "Competitive Inhibition":
                velocity_inhibitor = vmax * substrate / (km * alpha + substrate)
            elif mechanism == "Non-competitive Inhibition":
                velocity_inhibitor = (vmax / alpha) * substrate / (km + substrate)
            elif mechanism == "Uncompetitive Inhibition":
                velocity_inhibitor = (vmax / alpha) * substrate / (km / alpha + substrate)
            else:  # Mixed Inhibition
                alpha_prime = alpha * 0.75  # Non-competitive component (75% of alpha)
                velocity_inhibitor = (vmax / alpha_prime) * substrate / ((km * alpha / alpha_prime) + substrate)
            
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
    
    # User Guide for Calculator
    with st.expander("📖 How to Use the Calculator Tools", expanded=False):
        st.markdown("""
        **Three Powerful Calculation Tools - Choose a Tab:**
        
        ---
        
        ### **Tab 1: IC50 Calculator** 🧪
        
        **Purpose:** Calculate IC50 from your experimental data
        
        **Step-by-Step:**
        1. Select number of data points (3-10 recommended, 5 is good start)
        2. Enter your **inhibitor concentrations** (μM) in left column
        3. Enter corresponding **activity percentages** (0-100%) in right column
        4. IC50 is automatically calculated and displayed with a dose-response curve
        5. Click **"Download Results as CSV"** to export data
        
        **Tips:**
        - Your data must cross 50% activity for calculation to work
        - Use log-spaced concentrations for better curve fitting (e.g., 0.1, 1, 10, 100)
        - Duplicate concentrations will trigger a warning
        - Lower IC50 = more potent inhibitor
        
        **Interpreting Results:**
        - Very potent: < 0.1 μM
        - Potent: 0.1-1 μM
        - Moderate: 1-10 μM
        - Weak: > 10 μM
        
        ---
        
        ### **Tab 2: Ki Calculator** ⚖️
        
        **Purpose:** Convert IC50 to Ki (inhibition constant) using Cheng-Prusoff equations
        
        **Step-by-Step:**
        1. Select your **inhibition type** (Competitive, Non-competitive, or Uncompetitive)
        2. Enter your **IC50 value** (μM) from experiments
        3. Enter **substrate concentration [S]** (μM) used in assay
        4. Enter **Km value** (μM) for your enzyme
        5. Ki is automatically calculated with the appropriate formula
        
        **Formulas Used:**
        - Competitive: Ki = IC50 / (1 + [S]/Km)
        - Non-competitive: Ki = IC50
        - Uncompetitive: Ki = IC50 / (1 + Km/[S])
        
        **Important:** Units must be consistent (all μM or all nM)
        
        ---
        
        ### **Tab 3: Dose-Response Curve Generator** 📈
        
        **Purpose:** Generate theoretical dose-response curves for presentations or teaching
        
        **Step-by-Step:**
        1. Set **Top Activity** (usually 100% for no inhibitor)
        2. Set **Bottom Activity** (usually 0% for complete inhibition)
        3. Enter desired **IC50 value** (μM)
        4. Adjust **Hill Slope** (1.0 is standard, higher = steeper curve)
        5. Set **Max Concentration** range for X-axis
        6. View generated curve instantly
        
        **Uses:**
        - Creating example curves for presentations
        - Understanding Hill equation behavior
        - Comparing different IC50 values visually
        
        ---
        
        **All calculators provide instant results as you adjust parameters!**
        """)
    
    st.write("""Calculate inhibition constants and understand drug potency metrics.""")
    
    # Create tabs for different calculators
    tab1, tab2, tab3 = st.tabs(["IC50 Calculator", "Ki Calculator", "Dose-Response Curve"])
    
    with tab1:
        st.subheader("IC50 Calculator")
        st.write("**IC50** (Half maximal inhibitory concentration): The concentration of inhibitor required to reduce enzyme activity by 50%.")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("#### Input Parameters")
            
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
                # Validate input data
                conc_array = np.array(concentrations)
                act_array = np.array(activities)
                
                # Check for duplicate concentrations
                unique_conc = np.unique(conc_array)
                if len(unique_conc) < len(conc_array):
                    st.warning("⚠️ Warning: Duplicate concentration values detected. This may affect curve fitting accuracy.")
                
                # Fit dose-response curve (Hill equation)
                # y = Bottom + (Top - Bottom) / (1 + (IC50/x)^HillSlope)
                # Simplified: assume Hill slope = 1
                
                # Calculate IC50 using interpolation
                # Sort by concentration
                sorted_indices = np.argsort(conc_array)
                conc_sorted = conc_array[sorted_indices]
                act_sorted = act_array[sorted_indices]
                
                # Find IC50 (50% activity)
                if len(act_sorted) > 1 and act_sorted.max() > 50 and act_sorted.min() < 50:
                    # Check if activities decrease with concentration (typical inhibition)
                    if act_sorted[0] > act_sorted[-1]:
                        # Activities decrease: reverse for interpolation
                        ic50 = np.interp(50, act_sorted[::-1], conc_sorted[::-1])
                    else:
                        # Activities increase: interpolate directly
                        ic50 = np.interp(50, act_sorted, conc_sorted)
                    
                    st.success(f"### IC50 = {ic50:.2f} µM")
                    
                    # Generate smooth curve
                    max_conc = max(max(concentrations), 1.0)  # Ensure minimum range
                    conc_smooth = np.linspace(0.01, max_conc*1.2, 100)
                    # Hill equation with calculated IC50
                    hill_slope = 1.0
                    
                    # Determine top and bottom based on actual data pattern
                    if act_sorted[0] > act_sorted[-1]:
                        # Activities decrease with concentration (typical inhibition)
                        top = act_sorted[0]
                        bottom = act_sorted[-1]
                        # Hill equation for inhibition: activity decreases
                        act_smooth = bottom + (top - bottom) / (1 + (conc_smooth/ic50)**hill_slope)
                    else:
                        # Activities increase with concentration (atypical)
                        top = act_sorted[-1]
                        bottom = act_sorted[0]
                        # Inverse Hill equation
                        act_smooth = bottom + (top - bottom) / (1 + (ic50/conc_smooth)**hill_slope)
                    
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
                    
                    # Export data option
                    results_df = pd.DataFrame({
                        'Concentration_uM': conc_sorted,
                        'Activity_percent': act_sorted,
                        'IC50_uM': [ic50] * len(conc_sorted)
                    })
                    csv = results_df.to_csv(index=False)
                    st.download_button(
                        label="📅 Download Results as CSV",
                        data=csv,
                        file_name="ic50_results.csv",
                        mime="text/csv"
                    )
                else:
                    st.warning("⚠️ **Cannot calculate IC50:** Data must cross the 50% activity threshold. "
                             f"Current range: {act_sorted.min():.1f}% to {act_sorted.max():.1f}%. "
                             "Please adjust your data points to include values both above and below 50%.")
    
    with tab2:
        st.subheader("Ki Calculator (Inhibition Constant)")
        st.write("**Ki**: Dissociation constant of the enzyme-inhibitor complex. Lower Ki = Stronger binding.")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("#### Input Parameters")
            
            st.info("💡 **Important:** Ensure IC50, [S], and Km are all in the same units (µM).")
            
            inhibition_type = st.selectbox(
                "Inhibition Type",
                ["Competitive", "Non-competitive", "Uncompetitive"],
                key="ki_type"
            )
            
            ic50_input = st.number_input("IC50 (µM)", min_value=0.01, value=5.0, step=0.1,
                                        help="From IC50 assay")
            substrate_conc = st.number_input("[S] Substrate Concentration (µM)", 
                                           min_value=0.01, value=10.0, step=1.0,
                                           help="Actual substrate concentration used in the assay")
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
                # For uncompetitive inhibition: Ki' = IC50 / (1 + [S]/Km)
                ki = ic50_input / (1 + substrate_conc / km_input)
                st.success(f"### Ki = {ki:.3f} µM")
                
                st.latex(r"K_i = \frac{IC_{50}}{1 + \frac{[S]}{K_m}}")
                
                st.info(f"""**Uncompetitive Inhibition**
- IC50 = {ic50_input} µM
- Km = {km_input} µM  
- [S] = {substrate_conc} µM
- **Ki = {ki:.3f} µM**

The inhibitor only binds to the enzyme-substrate complex (ES).
Ki represents the dissociation constant for the ESI complex.
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
            
            # Validation
            if bottom_activity > top_activity:
                st.warning("⚠️ Bottom activity is greater than top activity. Curve will be inverted.")
            
            # Generate dose-response curve
            concentrations_curve = np.logspace(-3, np.log10(conc_range_max), 100)
            
            # Hill equation: y = Bottom + (Top - Bottom) / (1 + (x/IC50)^HillSlope)
            response = bottom_activity + (top_activity - bottom_activity) / \
                      (1 + (concentrations_curve / ic50_curve)**hill_slope)
            
            # Calculate actual IC50 activity level (midpoint)
            ic50_activity_level = (top_activity + bottom_activity) / 2
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=concentrations_curve, y=response, mode='lines',
                                   line=dict(color='purple', width=3)))
            fig.add_hline(y=ic50_activity_level, line_dash="dash", line_color="green",
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
    
    # User Guide for References
    with st.expander("📖 How to Use This Reference Library", expanded=False):
        st.markdown("""
        **Comprehensive Scientific References**
        
        All data, equations, and claims in this application are supported by peer-reviewed scientific literature.
        
        **What's included in each tab:**
        
        ---
        
        **Tab 1: Key Papers** 📄
        - **26+ peer-reviewed publications** from top journals
        - All citations in **APA 7th edition format**
        - **DOI links** for direct access to papers
        - Organized by topic:
          - Enzyme inhibition theory
          - Statins and cholesterol drugs
          - HIV protease inhibitors
          - ACE inhibitors
          - Kinase inhibitors
          - COX-2 inhibitors
        
        **How to use:**
        - Click DOI links to access full papers (may require institutional access)
        - Copy citations for your own reports or publications
        - All references are authoritative sources from journals like NEJM, Nature, Blood, JAMA
        
        ---
        
        **Tab 2: Textbooks** 📚
        - **12 recommended textbooks** for deeper learning
        - Classic biochemistry references (Berg, Lehninger, Voet)
        - Specialized enzyme kinetics texts (Cornish-Bowden)
        - Medicinal chemistry and drug design books
        - Pharmacology references (Goodman & Gilman)
        
        **Includes specific chapter recommendations** for relevant topics
        
        ---
        
        **Tab 3: Online Resources** 🌐
        - **Free databases:**
          - PubChem (chemical structures)
          - DrugBank (drug information)
          - BRENDA (enzyme data)
          - ChEMBL (bioactivity data)
          - PDB (protein structures)
        - **Educational resources:**
          - Khan Academy tutorials
          - NCBI Bookshelf (free textbooks)
        - **Professional organizations**
        
        **All resources are free to access!**
        
        ---
        
        **Tab 4: Citation** 📝
        - How to cite this educational tool
        - GitHub repository link
        - License information
        
        ---
        
        **Note:** All references have been verified for accuracy and accessibility. DOI links are functional as of December 2025.
        """)
    
    st.write("""This educational tool is based on established principles in biochemistry and pharmacology.""")
    
    # Create tabs for different reference categories
    tab1, tab2, tab3, tab4 = st.tabs(["Key Papers", "Textbooks", "Online Resources", "Citation"])
    
    with tab1:
        st.subheader("📄 Key Research Papers")
        
        st.info("**Note:** All references are formatted in APA 7th edition style. DOI links are provided for journal articles and have been verified as functional. Textbooks typically do not have DOI identifiers.")
        
        st.markdown("""
#### Enzyme Inhibition Theory & Methods

1. Copeland, R. A. (2013). Evaluation of enzyme inhibitors in drug discovery: A guide for medicinal chemists and pharmacologists (2nd ed.). *Methods of Biochemical Analysis*, *46*, 1–265. https://doi.org/10.1002/9781118540398

2. Cheng, Y., & Prusoff, W. H. (1973). Relationship between the inhibition constant (K1) and the concentration of inhibitor which causes 50 per cent inhibition (I50) of an enzymatic reaction. *Biochemical Pharmacology*, *22*(23), 3099–3108. https://doi.org/10.1016/0006-2952(73)90196-2

3. Segel, I. H. (1993). *Enzyme kinetics: Behavior and analysis of rapid equilibrium and steady-state enzyme systems*. Wiley-Interscience.

4. Cornish-Bowden, A. (2012). *Fundamentals of enzyme kinetics* (4th ed.). Wiley-Blackwell.

5. Morrison, J. F., & Walsh, C. T. (1988). The behavior and significance of slow-binding enzyme inhibitors. *Advances in Enzymology and Related Areas of Molecular Biology*, *61*, 201–301.

6. DiMasi, J. A., Grabowski, H. G., & Hansen, R. W. (2016). Innovation in the pharmaceutical industry: New estimates of R&D costs. *Journal of Health Economics*, *47*, 20–33. https://doi.org/10.1016/j.jhealeco.2016.01.012

#### Statins (HMG-CoA Reductase Inhibitors)

7. Istvan, E. S., & Deisenhofer, J. (2001). Structural mechanism for statin inhibition of HMG-CoA reductase. *Science*, *292*(5519), 1160–1164. https://doi.org/10.1126/science.1059344

8. Endo, A. (2010). A historical perspective on the discovery of statins. *Proceedings of the Japan Academy, Series B*, *86*(5), 484–493. https://doi.org/10.2183/pjab.86.484

9. Heart Protection Study Collaborative Group. (2002). MRC/BHF Heart Protection Study of cholesterol lowering with simvastatin in 20,536 high-risk individuals: A randomised placebo-controlled trial. *The Lancet*, *360*(9326), 7–22. https://doi.org/10.1016/S0140-6736(02)09327-3

10. Ridker, P. M., Danielson, E., Fonseca, F. A., Genest, J., Gotto, A. M., Jr., Kastelein, J. J., Koenig, W., Libby, P., Lorenzatti, A. J., MacFadyen, J. G., Nordestgaard, B. G., Shepherd, J., Willerson, J. T., & Glynn, R. J. (2008). Rosuvastatin to prevent vascular events in men and women with elevated C-reactive protein. *New England Journal of Medicine*, *359*(21), 2195–2207. https://doi.org/10.1056/NEJMoa0807646

#### HIV Protease Inhibitors

11. Palella, F. J., Jr., Delaney, K. M., Moorman, A. C., Loveless, M. O., Fuhrer, J., Satten, G. A., Aschman, D. J., & Holmberg, S. D. (1998). Declining morbidity and mortality among patients with advanced human immunodeficiency virus infection. *New England Journal of Medicine*, *338*(13), 853–860. https://doi.org/10.1056/NEJM199803263381301

12. Kohl, N. E., Emini, E. A., Schleif, W. A., Davis, L. J., Heimbach, J. C., Dixon, R. A., Scolnick, E. M., & Sigal, I. S. (1988). Active human immunodeficiency virus protease is required for viral infectivity. *Proceedings of the National Academy of Sciences*, *85*(13), 4686–4690. https://doi.org/10.1073/pnas.85.13.4686

13. Wlodawer, A., & Vondrasek, J. (1998). Inhibitors of HIV-1 protease: A major success of structure-assisted drug design. *Annual Review of Biophysics and Biomolecular Structure*, *27*, 249–284. https://doi.org/10.1146/annurev.biophys.27.1.249

14. Gulick, R. M., Mellors, J. W., Havlir, D., Eron, J. J., Gonzalez, C., McMahon, D., Richman, D. D., Valentine, F. T., Jonas, L., Meibohm, A., Emini, E. A., & Chodakewitz, J. A. (1997). Treatment with indinavir, zidovudine, and lamivudine in adults with human immunodeficiency virus infection and prior antiretroviral therapy. *New England Journal of Medicine*, *337*(11), 734–739. https://doi.org/10.1056/NEJM199709113371102

15. Flexner, C. (1998). HIV-protease inhibitors. *New England Journal of Medicine*, *338*(18), 1281–1293. https://doi.org/10.1056/NEJM199804303381808

#### ACE Inhibitors

16. Cushman, D. W., & Ondetti, M. A. (1991). History of the design of captopril and related inhibitors of angiotensin converting enzyme. *Hypertension*, *17*(4), 589–592. https://doi.org/10.1161/01.HYP.17.4.589

17. Pfeffer, M. A., Braunwald, E., Moyé, L. A., Basta, L., Brown, E. J., Jr., Cuddy, T. E., Davis, B. R., Geltman, E. M., Goldman, S., Flaker, G. C., Klein, M., Lamas, G. A., Packer, M., Rouleau, J., Rouleau, J. L., Rutherford, J., Wertheimer, J. H., & Hawkins, C. M. (1992). Effect of captopril on mortality and morbidity in patients with left ventricular dysfunction after myocardial infarction: Results of the survival and ventricular enlargement trial. *New England Journal of Medicine*, *327*(10), 669–677. https://doi.org/10.1056/NEJM199209033271001

18. Yusuf, S., Sleight, P., Pogue, J., Bosch, J., Davies, R., & Dagenais, G. (2000). Effects of an angiotensin-converting-enzyme inhibitor, ramipril, on cardiovascular events in high-risk patients. *New England Journal of Medicine*, *342*(3), 145–153. https://doi.org/10.1056/NEJM200001203420301

#### Kinase Inhibitors

19. Druker, B. J., Talpaz, M., Resta, D. J., Peng, B., Buchdunger, E., Ford, J. M., Lydon, N. B., Kantarjian, H., Capdeville, R., Ohno-Jones, S., & Sawyers, C. L. (2001). Efficacy and safety of a specific inhibitor of the BCR-ABL tyrosine kinase in chronic myeloid leukemia. *New England Journal of Medicine*, *344*(14), 1031–1037. https://doi.org/10.1056/NEJM200104053441401

20. Cohen, P. (2002). Protein kinases—the major drug targets of the twenty-first century? *Nature Reviews Drug Discovery*, *1*(4), 309–315. https://doi.org/10.1038/nrd773

21. Deininger, M., Buchdunger, E., & Druker, B. J. (2005). The development of imatinib as a therapeutic agent for chronic myeloid leukemia. *Blood*, *105*(7), 2640–2653. https://doi.org/10.1182/blood-2004-08-3097

22. Hochhaus, A., Larson, R. A., Guilhot, F., Radich, J. P., Branford, S., Hughes, T. P., Baccarani, M., Deininger, M. W., Cervantes, F., Fujihara, S., Ortmann, C. E., Menssen, H. D., Kantarjian, H., O'Brien, S. G., & Druker, B. J. (2017). Long-term outcomes of imatinib treatment for chronic myeloid leukemia. *New England Journal of Medicine*, *376*(10), 917–927. https://doi.org/10.1056/NEJMoa1609324

#### COX-2 Inhibitors

23. Vane, J. R., & Botting, R. M. (1998). Mechanism of action of nonsteroidal anti-inflammatory drugs. *American Journal of Medicine*, *104*(3A), 2S–8S. https://doi.org/10.1016/S0002-9343(97)00203-9

24. Bombardier, C., Laine, L., Reicin, A., Shapiro, D., Burgos-Vargas, R., Davis, B., Day, R., Ferraz, M. B., Hawkey, C. J., Hochberg, M. C., Kvien, T. K., & Schnitzer, T. J. (2000). Comparison of upper gastrointestinal toxicity of rofecoxib and naproxen in patients with rheumatoid arthritis. *New England Journal of Medicine*, *343*(21), 1520–1528. https://doi.org/10.1056/NEJM200011233432103

25. Silverstein, F. E., Faich, G., Goldstein, J. L., Simon, L. S., Pincus, T., Whelton, A., Makuch, R., Eisen, G., Agrawal, N. M., Stenson, W. F., Burr, A. M., Zhao, W. W., Kent, J. D., Lefkowith, J. B., Verburg, K. M., & Geis, G. S. (2000). Gastrointestinal toxicity with celecoxib vs nonsteroidal anti-inflammatory drugs for osteoarthritis and rheumatoid arthritis: The CLASS study. *JAMA*, *284*(10), 1247–1255. https://doi.org/10.1001/jama.284.10.1247

26. FitzGerald, G. A. (2004). Coxibs and cardiovascular disease. *New England Journal of Medicine*, *351*(17), 1709–1711. https://doi.org/10.1056/NEJMp048288
        """)
    
    with tab2:
        st.subheader("📖 Recommended Textbooks")
        
        st.markdown("""**All references in APA 7th edition format.**

#### Biochemistry & Enzyme Kinetics

1. Berg, J. M., Tymoczko, J. L., Gatto, G. J., Jr., & Stryer, L. (2019). *Biochemistry* (9th ed.). W. H. Freeman and Company.
   - Chapter 8: Enzymes: Basic Concepts and Kinetics
   - Chapter 9: Catalytic Strategies
   - Chapter 10: Regulatory Strategies

2. Nelson, D. L., & Cox, M. M. (2021). *Lehninger principles of biochemistry* (8th ed.). W. H. Freeman.
   - Chapter 6: Enzymes
   - Chapter 12: Biosignaling

3. Voet, D., Voet, J. G., & Pratt, C. W. (2016). *Fundamentals of biochemistry: Life at the molecular level* (5th ed.). Wiley.
   - Chapter 12: Enzyme Kinetics, Inhibition, and Control

4. Cornish-Bowden, A. (2012). *Fundamentals of enzyme kinetics* (4th ed.). Wiley-Blackwell.

5. Price, N. C., & Stevens, L. (1999). *Fundamentals of enzymology: The cell and molecular biology of catalytic proteins* (3rd ed.). Oxford University Press.

#### Drug Design & Medicinal Chemistry

6. Silverman, R. B., & Holladay, M. W. (2014). *The organic chemistry of drug design and drug action* (3rd ed.). Academic Press.
   - Chapter 4: Enzyme Inhibition and Inactivation
   - Chapter 5: Receptor Binding and Drug Design

7. Copeland, R. A. (2013). *Evaluation of enzyme inhibitors in drug discovery: A guide for medicinal chemists and pharmacologists* (2nd ed.). Wiley-Interscience.

8. Patrick, G. L. (2017). *An introduction to medicinal chemistry* (6th ed.). Oxford University Press.
   - Chapter 7: Enzymes as Drug Targets
   - Chapter 20: Anticancer Agents

9. Wermuth, C. G., Aldous, D., Raboisson, P., & Rognan, D. (Eds.). (2015). *The practice of medicinal chemistry* (4th ed.). Academic Press.

10. Lemke, T. L., Williams, D. A., Roche, V. F., & Zito, S. W. (2020). *Foye's principles of medicinal chemistry* (8th ed.). Wolters Kluwer.

#### Pharmacology

11. Brunton, L. L., Hilal-Dandan, R., & Knollmann, B. C. (Eds.). (2018). *Goodman & Gilman's: The pharmacological basis of therapeutics* (13th ed.). McGraw-Hill Education.

12. Katzung, B. G., & Trevor, A. J. (2021). *Basic & clinical pharmacology* (15th ed.). McGraw-Hill Education.
   - Section II: Autonomic Drugs
   - Section V: Drugs That Act in the Central Nervous System
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
        
        st.markdown("""If you use this interactive educational tool in your teaching, research, or publications, 
please cite it using one of the following formats:""")
        
        st.markdown("#### APA 7th Edition (Recommended)")
        st.code("""Chiu, K. (2025). Enzyme inhibitors in drug development: An interactive educational tool [Interactive web application]. Streamlit. https://lab-kason-biochem.streamlit.app""", language="text")
        
        st.markdown("#### MLA 9th Edition")
        st.code("""Chiu, Kason. "Enzyme Inhibitors in Drug Development: An Interactive Educational Tool." Streamlit, 2025, lab-kason-biochem.streamlit.app. Accessed [Date].""", language="text")
        
        st.markdown("#### Chicago 17th Edition")
        st.code("""Chiu, Kason. 2025. "Enzyme Inhibitors in Drug Development: An Interactive Educational Tool." Interactive web application. Streamlit. https://lab-kason-biochem.streamlit.app.""", language="text")
        
        st.markdown("#### BibTeX Format")
        st.code("""@misc{chiu2025enzyme,
  title={Enzyme Inhibitors in Drug Development: An Interactive Educational Tool},
  author={Chiu, Kason},
  year={2025},
  howpublished={\\url{https://lab-kason-biochem.streamlit.app}},
  note={Interactive web application built with Streamlit and Python. 
        Includes IC50/Ki calculators, kinetics simulators, and drug case studies}
}""", language="bibtex")
        
        st.markdown("---")
        
        st.markdown("#### Additional Information")
        st.info("""**Version:** 1.0 (November 2025)
**Platform:** Streamlit 1.51.0
**Technologies:** Python, NumPy, Pandas, Plotly
**License:** Educational use permitted with attribution
**Repository:** https://github.com/lab-Kason/biochem
        """)
        
        st.markdown("---")
        
        st.warning("""**Disclaimer:** This tool is designed for educational purposes only. 
All data, calculations, and clinical information should be verified with primary 
literature sources and are not intended for clinical or commercial drug development 
decisions. The case study data presented are approximate values derived from published 
literature for illustrative purposes. Always consult current medical literature and 
guidelines for clinical applications.""")

# Case Studies Section
def show_case_studies():
    st.markdown('<div class="section-header">💊 Successful Drug Case Studies</div>', unsafe_allow_html=True)
    
    # User Guide for Case Studies
    with st.expander("📖 How to Explore Case Studies", expanded=False):
        st.markdown("""
        **Real-World Drug Success Stories**
        
        **What's included:**
        This section showcases **5 blockbuster enzyme inhibitor drugs** that revolutionized medicine.
        
        **How to navigate:**
        1. **Select a drug class** from the radio buttons above
        2. **Read the drug card** (left panel) for key facts:
           - Drug name and target enzyme
           - Inhibition mechanism
           - FDA approval date
           - Market impact
        3. **Study "How It Works"** section to understand the mechanism
        4. **Analyze interactive charts** (right panel) showing:
           - Clinical efficacy data
           - Patient outcomes over time
           - Comparative potency
        
        **The 5 Case Studies:**
        
        **1. Statins (Cholesterol)** 💊
        - Example: Lipitor (Atorvastatin)
        - Best-selling drug of all time
        - Reduces cardiac deaths significantly
        
        **2. HIV Protease Inhibitors** 🦠
        - Transformed HIV from fatal to manageable
        - Increased life expectancy from ~1 year to near-normal
        - Structure-based drug design success story
        
        **3. ACE Inhibitors (Blood Pressure)** ❤️
        - Examples: Lisinopril, Enalapril
        - One of most prescribed drug classes
        - Inspired by snake venom peptides
        
        **4. Kinase Inhibitors (Cancer)** 🎗️
        - Example: Gleevec (Imatinib)
        - Revolutionized cancer treatment
        - 95% remission rate in CML
        
        **5. COX-2 Inhibitors (Pain)** 🩹
        - Selective pain relief
        - Reduced GI side effects vs traditional NSAIDs
        - Example of targeted drug design
        
        **Each case includes:**
        - ✅ Mechanism of action
        - ✅ Clinical trial data
        - ✅ Interactive visualizations
        - ✅ Real-world impact statistics
        - ✅ References to original research
        
        *Click through each drug to see how enzyme inhibitor design led to life-saving medications!*
        """)
    
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
            <b>First Approved:</b> 1996 (FDA)<br>
            <b>Peak Sales:</b> $12.9 billion/year<br>
            <b>Impact:</b> Best-selling drug of all time
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("""**How It Works:**
Statins competitively inhibit HMG-CoA reductase, the rate-limiting enzyme in cholesterol 
biosynthesis. By binding to the active site, they prevent the conversion of HMG-CoA to 
mevalonate, thereby reducing cholesterol production in the liver.

**Clinical Impact:**
- Reduces LDL cholesterol by 39-60%
- Decreases cardiovascular events by 25-35%
- Prevents ~10,000 deaths annually in the US alone

**Key Studies:** Heart Protection Study (2002), JUPITER Trial (2008)
            """)
        
        with col2:
            # Efficacy chart
            data = {'Year': [1990, 1995, 2000, 2005, 2010, 2015, 2020],
                   'Cardiac Deaths per 100,000': [321, 298, 257, 216, 179, 155, 134]}
            df = pd.DataFrame(data)
            fig = px.line(df, x='Year', y='Cardiac Deaths per 100,000', 
                         title="Impact of Statins on Cardiac Mortality (US)",
                         markers=True)
            fig.update_layout(height=350)
            st.plotly_chart(fig, width='stretch')
            st.caption("""*Data source: CDC Wonder Database, Age-Adjusted Death Rates. 
            Trends correlate with statin introduction and widespread adoption (Istvan & Deisenhofer, 2001; 
            Heart Protection Study, 2002). See References section for full citations.*""")
            
            # Market comparison
            statin_data = pd.DataFrame({
                'Drug': ['Atorvastatin', 'Rosuvastatin', 'Simvastatin', 'Pravastatin'],
                'LDL Reduction (%)': [39, 45, 35, 28],
                'Type': ['High Potency', 'High Potency', 'Moderate', 'Moderate']
            })
            fig2 = px.bar(statin_data, x='Drug', y='LDL Reduction (%)', 
                         color='Type', title='Comparative Potency of Statins')
            fig2.update_layout(height=300)
            st.plotly_chart(fig2, width='stretch')
            st.caption("""*Representative data based on clinical trial meta-analyses. 
            LDL-C reduction at standard doses (Istvan & Deisenhofer, 2001; Heart Protection Study, 2002; 
            Ridker et al., 2008). See References section for full citations.*""")
    
    elif case_study == "HIV Protease Inhibitors":
        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("""
            <div class="drug-card">
            <h4>💊 Ritonavir, Saquinavir, Indinavir</h4>
            <b>Target:</b> HIV-1 Protease<br>
            <b>Mechanism:</b> Competitive inhibition<br>
            <b>Disease:</b> HIV/AIDS<br>
            <b>First Approved:</b> 1995 (Saquinavir)<br>
            <b>Impact:</b> Transformed HIV from fatal to manageable chronic disease
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("""**How It Works:**
HIV protease inhibitors mimic the transition state of the natural peptide substrate, 
binding tightly to the enzyme's active site. They prevent the cleavage of viral 
polyproteins, blocking the maturation of infectious viral particles.

**Clinical Impact:**
- Reduces viral load by >90% when combined with other antiretrovirals
- Increased life expectancy from ~1 year to near-normal
- Death rate decreased by 80% after introduction (1996-1998)
- Part of HAART (Highly Active Antiretroviral Therapy)

**Design Success:** Structure-based drug design using X-ray crystallography
            """)
        
        with col2:
            # HIV survival timeline
            survival_data = pd.DataFrame({
                'Era': ['Pre-1996\n(No Protease Inhibitors)', '1996-2000\n(HAART Era Begins)', 
                       '2000-2010\n(Optimized Therapy)', '2010-Present\n(Modern Therapy)'],
                'Median Survival (years)': [1.5, 8, 20, 35],
                'Order': [1, 2, 3, 4]
            })
            fig = px.bar(survival_data, x='Era', y='Median Survival (years)',
                        title='HIV Survival: Impact of Protease Inhibitors',
                        color='Median Survival (years)',
                        color_continuous_scale='Viridis')
            fig.update_layout(height=350, showlegend=False)
            st.plotly_chart(fig, width='stretch')
            st.caption("""*Data sources: Palella et al. (1998) NEJM; Antiretroviral Therapy Cohort Collaboration (2008); 
            UNAIDS Global AIDS Update (2020). See References section for full citations.*""")
            
            # Drug potency comparison
            st.markdown("**Protease Inhibitor Potency (IC50 values):**")
            pi_data = pd.DataFrame({
                'Drug': ['Ritonavir', 'Saquinavir', 'Indinavir', 'Lopinavir'],
                'IC50 (nM)': [15, 0.4, 0.56, 1.3]
            })
            fig2 = px.bar(pi_data, x='Drug', y='IC50 (nM)', 
                         title='Lower IC50 = More Potent',
                         log_y=True)
            fig2.update_layout(height=300)
            st.plotly_chart(fig2, width='stretch')
            st.caption("""*Data source: Flexner (1998) HIV-protease inhibitors. NEJM 338(18):1281-1293. 
            See References section for full citation and DOI.*""")
    
    elif case_study == "ACE Inhibitors (Blood Pressure)":
        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("""
            <div class="drug-card">
            <h4>💊 Lisinopril, Enalapril, Captopril</h4>
            <b>Target:</b> Angiotensin-Converting Enzyme (ACE)<br>
            <b>Mechanism:</b> Competitive inhibition<br>
            <b>Disease:</b> Hypertension, Heart Failure<br>
            <b>First Approved:</b> 1981 (Captopril)<br>
            <b>Impact:</b> One of most prescribed drug classes worldwide
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("""**How It Works:**
ACE inhibitors block the conversion of angiotensin I to angiotensin II, a potent 
vasoconstrictor. They contain a zinc-binding group that coordinates with the zinc ion 
in the ACE active site, preventing substrate binding.

**Clinical Impact:**
- Reduces blood pressure by 10-15 mmHg (systolic)
- Decreases heart failure mortality by 20-30%
- Protects kidney function in diabetic patients
- Used by >40 million Americans annually

**Design Inspiration:** Based on snake venom peptides (Bothrops jararaca)
            """)
        
        with col2:
            # Blood pressure reduction
            bp_data = pd.DataFrame({
                'Weeks': [0, 2, 4, 8, 12],
                'Systolic BP (mmHg)': [160, 152, 145, 138, 135],
                'Diastolic BP (mmHg)': [98, 94, 90, 86, 84]
            })
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=bp_data['Weeks'], y=bp_data['Systolic BP (mmHg)'],
                                   name='Systolic', mode='lines+markers', line=dict(color='red')))
            fig.add_trace(go.Scatter(x=bp_data['Weeks'], y=bp_data['Diastolic BP (mmHg)'],
                                   name='Diastolic', mode='lines+markers', line=dict(color='blue')))
            fig.update_layout(title='Typical Blood Pressure Response to ACE Inhibitors',
                            xaxis_title='Treatment Duration (weeks)',
                            yaxis_title='Blood Pressure (mmHg)',
                            height=350)
            st.plotly_chart(fig, width='stretch')
            st.caption("""*Representative data based on clinical trials: ALLHAT (2002), HOPE Study (Yusuf et al., 2000). 
            See References section for full citations.*""")
            
            # Cardiovascular outcomes
            outcome_data = pd.DataFrame({
                'Outcome': ['Heart Attack', 'Stroke', 'Heart Failure', 'Death'],
                'Risk Reduction (%)': [20, 25, 30, 23]
            })
            fig2 = px.bar(outcome_data, x='Outcome', y='Risk Reduction (%)',
                         title='Cardiovascular Risk Reduction with ACE Inhibitors',
                         color='Risk Reduction (%)', color_continuous_scale='Greens')
            fig2.update_layout(height=300, showlegend=False)
            st.plotly_chart(fig2, width='stretch')
            st.caption("""*Data compiled from: SAVE Trial (Pfeffer et al., 1992), HOPE Study (Yusuf et al., 2000), 
            SOLVD Trial. Meta-analysis values averaged across major trials. See References section for full citations.*""")
    
    elif case_study == "Kinase Inhibitors (Cancer)":
        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("""
            <div class="drug-card">
            <h4>💊 Imatinib (Gleevec), Gefitinib, Erlotinib</h4>
            <b>Target:</b> Tyrosine Kinases (BCR-ABL, EGFR)<br>
            <b>Mechanism:</b> Competitive inhibition (ATP-binding site)<br>
            <b>Disease:</b> Chronic Myeloid Leukemia, Lung Cancer<br>
            <b>First Approved:</b> 2001 (Imatinib)<br>
            <b>Impact:</b> Paradigm shift toward targeted cancer therapy
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("""**How It Works:**
Kinase inhibitors compete with ATP for the enzyme's binding site, preventing 
phosphorylation of target proteins. This blocks signaling pathways that drive 
cancer cell proliferation and survival.

**Clinical Impact - Imatinib for CML:**
- 10-year survival rate: 83% (vs. 20% before 2001)
- Complete cytogenetic response: 87% of patients
- Transformed CML from terminal to chronic condition
- Led to development of 50+ kinase inhibitor drugs

**Precision Medicine:** First major success of targeted molecular therapy
            """)
        
        with col2:
            # CML survival comparison
            cml_data = pd.DataFrame({
                'Treatment': ['Pre-Imatinib\n(1990s)', 'Imatinib\n(2001+)', 'Second-gen\n(2006+)'],
                '5-Year Survival (%)': [30, 89, 93],
                'Order': [1, 2, 3]
            })
            fig = px.bar(cml_data, x='Treatment', y='5-Year Survival (%)',
                        title='CML Survival: The Imatinib Revolution',
                        color='5-Year Survival (%)',
                        color_continuous_scale='Blues')
            fig.update_layout(height=350, showlegend=False)
            st.plotly_chart(fig, width='stretch')
            st.caption("""*Data sources: Druker et al. (2001) NEJM; Hochhaus et al. (2017) 10-year follow-up study. 
            See References section for full citations.*""")
            
            # Kinase inhibitor selectivity
            st.markdown("**Selectivity Profile:**")
            selectivity_data = pd.DataFrame({
                'Target': ['BCR-ABL', 'PDGFR', 'c-KIT', 'Off-targets'],
                'IC50 (nM)': [260, 380, 410, 5000],
                'Type': ['Primary', 'Secondary', 'Secondary', 'Non-target']
            })
            fig2 = px.bar(selectivity_data, x='Target', y='IC50 (nM)',
                         title='Imatinib Selectivity (Lower = More Potent)',
                         color='Type', log_y=True)
            fig2.update_layout(height=300)
            st.plotly_chart(fig2, width='stretch')
            st.caption("""*Data source: Deininger et al. (2005) The development of imatinib as a therapeutic agent. 
            Blood 105(7):2640-2653. See References section for full citation.*""")
    
    else:  # COX-2 Inhibitors (Pain)
        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("""
            <div class="drug-card">
            <h4>💊 Celecoxib (Celebrex), Rofecoxib (Vioxx)</h4>
            <b>Target:</b> Cyclooxygenase-2 (COX-2)<br>
            <b>Mechanism:</b> Selective competitive inhibition<br>
            <b>Disease:</b> Pain, Inflammation, Arthritis<br>
            <b>First Approved:</b> 1998 (Celecoxib)<br>
            <b>Impact:</b> Safer alternative to traditional NSAIDs
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("""**How It Works:**
COX-2 inhibitors selectively block cyclooxygenase-2, the enzyme induced during 
inflammation, while sparing COX-1 (important for stomach protection). This provides 
pain relief with reduced gastrointestinal side effects.

**Clinical Impact:**
- Similar pain relief to traditional NSAIDs
- 50-60% reduction in serious GI complications
- Reduced gastric ulcers: 0.4% vs. 1.4% (traditional NSAIDs)
- Used by millions for arthritis management

**Selectivity:** COX-2/COX-1 selectivity ratio >100:1

**Note:** Rofecoxib withdrawn in 2004 due to cardiovascular risks
            """)
        
        with col2:
            # Selectivity comparison
            cox_data = pd.DataFrame({
                'Drug': ['Celecoxib', 'Traditional NSAIDs', 'Aspirin'],
                'COX-2 Inhibition': [90, 75, 65],
                'COX-1 Inhibition': [10, 80, 95]
            })
            fig = go.Figure()
            fig.add_trace(go.Bar(name='COX-2 Inhibition', x=cox_data['Drug'], 
                               y=cox_data['COX-2 Inhibition'], marker_color='green'))
            fig.add_trace(go.Bar(name='COX-1 Inhibition', x=cox_data['Drug'], 
                               y=cox_data['COX-1 Inhibition'], marker_color='red'))
            fig.update_layout(title='COX-2 Selectivity: Reducing GI Side Effects',
                            yaxis_title='% Inhibition',
                            barmode='group',
                            height=350)
            st.plotly_chart(fig, width='stretch')
            st.caption("""*Data source: Vane & Botting (1998) Mechanism of action of NSAIDs. 
            American Journal of Medicine 104(3A):2S-8S. See References section for full citation.*""")
            
            # Side effect comparison
            side_effects = pd.DataFrame({
                'Side Effect': ['GI Ulcers', 'GI Bleeding', 'Dyspepsia'],
                'Traditional NSAIDs (%)': [1.4, 1.0, 15],
                'COX-2 Inhibitors (%)': [0.4, 0.3, 8]
            })
            fig2 = go.Figure()
            fig2.add_trace(go.Bar(name='Traditional NSAIDs', 
                                x=side_effects['Side Effect'], 
                                y=side_effects['Traditional NSAIDs (%)']))
            fig2.add_trace(go.Bar(name='COX-2 Inhibitors', 
                                x=side_effects['Side Effect'], 
                                y=side_effects['COX-2 Inhibitors (%)']))
            fig2.update_layout(title='Gastrointestinal Side Effects Comparison',
                             yaxis_title='Incidence (%)',
                             barmode='group',
                             height=300)
            st.plotly_chart(fig2, width='stretch')
            st.caption("""*Data sources: CLASS Study (Silverstein et al., 2000) & VIGOR Trial (Bombardier et al., 2000). 
            See References section for full citations.*""")

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
        
        # User Guide for Kinetics
        with st.expander("📖 How to Use the Kinetics Explorer", expanded=False):
            st.markdown("""
            **Advanced Enzyme Kinetics Analysis Tools**
            
            **Three Different Analysis Methods - Choose a Tab:**
            
            ---
            
            ### **Tab 1: Lineweaver-Burk Plot** 📐
            
            **Purpose:** Linear transformation of Michaelis-Menten equation to determine Km and Vmax
            
            **How to use:**
            1. Adjust **Vmax slider** (10-200 μmol/min) - maximum velocity
            2. Adjust **Km slider** (0.5-20 mM) - substrate affinity
            3. Check **"Add Inhibitor"** to compare with/without inhibition
            4. Select **inhibition type** (Competitive, Non-competitive, Uncompetitive)
            5. Set **Ki value** (inhibitor dissociation constant)
            6. Set **[I]** (inhibitor concentration)
            
            **Interpreting the plot:**
            - **Y-intercept = 1/Vmax**
            - **X-intercept = -1/Km**
            - **Slope = Km/Vmax**
            
            **Distinguishing inhibition types:**
            - **Competitive:** Lines intersect on Y-axis (same Vmax, different Km)
            - **Non-competitive:** Lines intersect on X-axis (different Vmax, same Km)
            - **Uncompetitive:** Parallel lines (both Vmax and Km change proportionally)
            
            ---
            
            ### **Tab 2: Eadie-Hofstee Plot** 📊
            
            **Purpose:** Alternative linear transformation with more uniform error distribution
            
            **How to use:**
            1. Adjust **Vmax** and **Km** parameters
            2. View linear relationship: v vs v/[S]
            
            **Advantages:**
            - Better statistical properties than Lineweaver-Burk
            - More uniform error distribution across the plot
            - Y-intercept directly gives Vmax
            - Slope = -Km
            
            **When to use:** Preferred for accurate parameter determination from experimental data
            
            ---
            
            ### **Tab 3: Kinetic Simulator** 🧪
            
            **Purpose:** Comprehensive Michaelis-Menten simulation with calculated parameters
            
            **How to use:**
            1. Set **Vmax** (μmol/min) - maximum reaction velocity
            2. Set **Km** (mM) - Michaelis constant
            3. Set **[E]** enzyme concentration (μM)
            4. Adjust **substrate range** (1-100 mM)
            5. View calculated **kcat** and **kcat/Km** values
            
            **What you'll see:**
            - **Upper plot:** Michaelis-Menten curve with Vmax and Km indicators
            - **Lower plot:** Enzyme saturation percentage vs [S]
            
            **Calculated constants:**
            - **kcat (turnover number):** Number of substrate molecules converted per enzyme per minute
            - **kcat/Km (specificity constant):** Measure of catalytic efficiency
            - Comparison to diffusion limit (~10⁸-10⁹ M⁻¹s⁻¹)
            
            **Applications:**
            - Understanding enzyme saturation
            - Comparing enzyme efficiency
            - Visualizing substrate concentration effects
            
            ---
            
            **All plots update in real-time as you adjust parameters!**
            """)
        
        st.write("""Explore enzyme kinetics through interactive plots and simulations. 
        Visualize Michaelis-Menten and Lineweaver-Burk plots under different inhibition conditions.""")
        
        # Create tabs for different kinetic analyses
        tab1, tab2, tab3 = st.tabs(["Lineweaver-Burk Plot", "Eadie-Hofstee Plot", "Kinetic Simulator"])
        
        with tab1:
            st.subheader("Lineweaver-Burk (Double Reciprocal) Plot")
            
            col1, col2 = st.columns([1, 2])
            
            with col1:
                st.markdown("#### Parameters")
                
                vmax_lb = st.slider("Vmax (µmol/min)", 10, 200, 100, 5, key="vmax_lb")
                km_lb = st.slider("Km (mM)", 0.5, 20.0, 5.0, 0.5, key="km_lb")
                
                show_inhibitor_lb = st.checkbox("Add Inhibitor", value=True, key="show_inh_lb")
                
                if show_inhibitor_lb:
                    inh_type_lb = st.selectbox(
                        "Inhibition Type",
                        ["Competitive", "Non-competitive", "Uncompetitive"],
                        key="inh_type_lb"
                    )
                    
                    ki_lb = st.slider("Ki (mM)", 0.5, 10.0, 2.0, 0.5, key="ki_lb")
                    inhibitor_conc = st.slider("[I] Inhibitor Concentration (mM)", 
                                              0.0, 10.0, 3.0, 0.5, key="inh_conc_lb")
                
                st.markdown("""**Equation:**""")
                st.latex(r"\frac{1}{v} = \frac{K_m}{V_{max}} \cdot \frac{1}{[S]} + \frac{1}{V_{max}}")
                
                st.info("""**Interpretation:**
- **Y-intercept:** 1/Vmax
- **X-intercept:** -1/Km
- **Slope:** Km/Vmax
                """)
            
            with col2:
                # Generate substrate concentrations
                substrate_conc = np.array([0.5, 1, 2, 4, 8, 16, 32])  # mM
                
                # Calculate velocities (no inhibitor)
                velocity_no_inh = vmax_lb * substrate_conc / (km_lb + substrate_conc)
                
                # Lineweaver-Burk transformation (avoid division by zero)
                # Filter out any zero values
                if np.any(substrate_conc == 0) or np.any(velocity_no_inh == 0):
                    st.error("⚠️ Cannot create Lineweaver-Burk plot: division by zero detected")
                else:
                    reciprocal_s = 1 / substrate_conc
                    reciprocal_v_no_inh = 1 / velocity_no_inh
                    
                    fig = go.Figure()
                    
                    # Plot without inhibitor
                    fig.add_trace(go.Scatter(
                        x=reciprocal_s, 
                        y=reciprocal_v_no_inh,
                        mode='lines+markers',
                        name='No Inhibitor',
                        line=dict(color='blue', width=2),
                        marker=dict(size=8)
                    ))
                    
                    # Add inhibitor conditions
                    if show_inhibitor_lb:
                        alpha = 1 + (inhibitor_conc / ki_lb)
                        
                        if inh_type_lb == "Competitive":
                            velocity_inh = vmax_lb * substrate_conc / (km_lb * alpha + substrate_conc)
                        elif inh_type_lb == "Non-competitive":
                            velocity_inh = (vmax_lb / alpha) * substrate_conc / (km_lb + substrate_conc)
                        else:  # Uncompetitive
                            velocity_inh = (vmax_lb / alpha) * substrate_conc / (km_lb / alpha + substrate_conc)
                        
                        reciprocal_v_inh = 1 / velocity_inh
                        
                        fig.add_trace(go.Scatter(
                            x=reciprocal_s, 
                            y=reciprocal_v_inh,
                            mode='lines+markers',
                            name=f'With {inh_type_lb} Inhibitor',
                            line=dict(color='red', width=2, dash='dash'),
                            marker=dict(size=8)
                        ))
                    
                    fig.update_layout(
                        title='Lineweaver-Burk Plot',
                        xaxis_title='1/[S] (1/mM)',
                        yaxis_title='1/v (min/µmol)',
                        height=500,
                        hovermode='closest'
                    )
                    
                    # Add grid
                    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
                    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
                    
                    st.plotly_chart(fig, width='stretch')
                    st.caption("*Theoretical simulation based on Michaelis-Menten enzyme kinetics equations (Segel, 1993).*")
                    
                    # Display calculated parameters
                    st.markdown("**Calculated Parameters:**")
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("Vmax", f"{vmax_lb} µmol/min")
                        st.metric("Km", f"{km_lb} mM")
                    with col_b:
                        st.metric("Vmax/Km (Efficiency)", f"{vmax_lb/km_lb:.2f}")
                        if show_inhibitor_lb:
                            st.metric("Ki", f"{ki_lb} mM")
        
        with tab2:
            st.subheader("Eadie-Hofstee Plot")
            
            col1, col2 = st.columns([1, 2])
            
            with col1:
                st.markdown("#### Parameters")
                
                vmax_eh = st.slider("Vmax (µmol/min)", 10, 200, 100, 5, key="vmax_eh")
                km_eh = st.slider("Km (mM)", 0.5, 20.0, 5.0, 0.5, key="km_eh")
                
                st.markdown("**Equation:**")
                st.latex(r"v = -K_m \cdot \frac{v}{[S]} + V_{max}")
                
                st.info("""**Interpretation:**
- **Y-intercept:** Vmax
- **Slope:** -Km
- More uniform error distribution than Lineweaver-Burk
                """)
            
            with col2:
                # Generate data
                substrate_conc_eh = np.linspace(0.5, 30, 50)
                velocity_eh = vmax_eh * substrate_conc_eh / (km_eh + substrate_conc_eh)
                v_over_s = velocity_eh / substrate_conc_eh
                
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=v_over_s,
                    y=velocity_eh,
                    mode='lines+markers',
                    name='Eadie-Hofstee',
                    line=dict(color='purple', width=2),
                    marker=dict(size=6)
                ))
                
                fig.update_layout(
                    title='Eadie-Hofstee Plot',
                    xaxis_title='v/[S] (µmol/min/mM)',
                    yaxis_title='v (µmol/min)',
                    height=500
                )
                
                st.plotly_chart(fig, width='stretch')
                st.caption("*Theoretical plot based on enzyme kinetics equations (Cornish-Bowden, 2012).*")
                
                st.markdown(f"""**Results:**
- **Vmax:** {vmax_eh} µmol/min (y-intercept)
- **Km:** {km_eh} mM (negative of slope)
- **Catalytic efficiency (kcat/Km):** {vmax_eh/km_eh:.2f}
                """)
        
        with tab3:
            st.subheader("Enzyme Kinetics Simulator")
            
            st.write("""Simulate enzyme reactions and observe how different parameters affect reaction velocity.""")
            
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown("#### Enzyme Properties")
                
                sim_vmax = st.number_input("Vmax (µmol/min)", min_value=1.0, value=100.0, step=5.0)
                sim_km = st.number_input("Km (mM)", min_value=0.1, value=5.0, step=0.5)
                enzyme_conc = st.number_input("[E] Enzyme Concentration (µM)", 
                                             min_value=0.1, value=1.0, step=0.1)
                
                st.markdown("#### Reaction Conditions")
                substrate_range = st.slider("Substrate Concentration Range (mM)", 
                                           1.0, 100.0, 50.0, 1.0)
                temperature = st.slider("Temperature (°C)", 20, 40, 37, 1,
                                      help="For reference only - not used in calculation")
                
                # Calculate kcat
                if enzyme_conc > 0:
                    kcat = sim_vmax / enzyme_conc
                    kcat_km = kcat / sim_km
                    st.markdown(f"""**Calculated Constants:**
- **kcat (turnover number):** {kcat:.1f} min⁻¹
- **kcat/Km (specificity constant):** {kcat_km:.2f} min⁻¹mM⁻¹
- **Diffusion limit:** ~10⁸ to 10⁹ M⁻¹s⁻¹
                    """)
                else:
                    st.warning("⚠️ Enzyme concentration must be > 0 to calculate kcat")
            
            with col2:
                # Generate simulation data
                substrate_sim = np.linspace(0.1, substrate_range, 200)
                velocity_sim = sim_vmax * substrate_sim / (sim_km + substrate_sim)
                
                # Create figure with multiple views
                fig = go.Figure()
                
                # Michaelis-Menten curve
                fig.add_trace(go.Scatter(
                    x=substrate_sim,
                    y=velocity_sim,
                    mode='lines',
                    name='Reaction Velocity',
                    line=dict(color='green', width=3)
                ))
                
                # Add Vmax line
                fig.add_hline(y=sim_vmax, line_dash="dash", line_color="red",
                            annotation_text=f"Vmax = {sim_vmax}")
                
                # Add Km indicator
                v_at_km = sim_vmax / 2
                fig.add_hline(y=v_at_km, line_dash="dot", line_color="orange",
                            annotation_text=f"Vmax/2 at Km = {sim_km}")
                fig.add_vline(x=sim_km, line_dash="dot", line_color="orange")
                
                fig.update_layout(
                    title='Michaelis-Menten Kinetics Simulation',
                    xaxis_title='[S] Substrate Concentration (mM)',
                    yaxis_title='v Reaction Velocity (µmol/min)',
                    height=400
                )
                
                st.plotly_chart(fig, width='stretch')
                st.caption("*Simulated using Michaelis-Menten equation: v = Vmax[S]/(Km + [S]).*")
                
                # Saturation analysis
                saturation = (substrate_sim / (sim_km + substrate_sim)) * 100
                
                fig2 = go.Figure()
                fig2.add_trace(go.Scatter(
                    x=substrate_sim,
                    y=saturation,
                    mode='lines',
                    fill='tonexty',
                    line=dict(color='blue', width=2)
                ))
                
                fig2.update_layout(
                    title='Enzyme Saturation',
                    xaxis_title='[S] Substrate Concentration (mM)',
                    yaxis_title='Enzyme Saturation (%)',
                    height=300
                )
                
                st.plotly_chart(fig2, width='stretch')
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
