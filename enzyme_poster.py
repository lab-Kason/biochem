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
    .pipeline-stage {
        background: linear-gradient(45deg, #2E86AB, #A23B72);
        color: white;
        padding: 1rem;
        border-radius: 10px;
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
            options=["Overview", "Mechanisms", "Case Studies", "Kinetics", "Pipeline", "Future Trends"],
            icons=["house", "gear", "book", "graph-up", "clock", "rocket"],
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

# ==================== PAGE 1: OVERVIEW ====================
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
        st.plotly_chart(fig, use_container_width=True)

# ==================== PAGE 2: MECHANISMS ====================
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
        st.plotly_chart(fig, use_container_width=True)

# ==================== PAGE 3: CASE STUDIES ====================
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
            data = {'Year': [1990, 2000, 2010, 2020],
                   'Cardiac Deaths per 1000': [8.5, 6.2, 4.1, 2.8]}
            df = pd.DataFrame(data)
            fig = px.line(df, x='Year', y='Cardiac Deaths per 1000', 
                         title="Impact on Cardiac Mortality")
            st.plotly_chart(fig, use_container_width=True)
    
    elif case_study == "HIV Protease Inhibitors":
        col1, col2 = st.columns([1, 1])
        with col1:
            st.markdown("""
            <div class="drug-card">
            <h4>💊 Ritonavir</h4>
            <b>Target:</b> HIV-1 protease<br>
            <b>Mechanism:</b> Competitive inhibition<br>
            <b>Disease:</b> HIV/AIDS<br>
            <b>Impact:</b> Transformed HIV into manageable chronic condition
            </div>
            """, unsafe_allow_html=True)
        with col2:
            data = {'Year': [1995, 2000, 2005, 2010, 2020],
                   'HIV Deaths (thousands)': [50, 22, 12, 8, 3]}
            df = pd.DataFrame(data)
            fig = px.line(df, x='Year', y='HIV Deaths (thousands)', 
                         title="Reduction in HIV Mortality")
            st.plotly_chart(fig, use_container_width=True)
    
    elif case_study == "ACE Inhibitors (Blood Pressure)":
        st.markdown("""
        <div class="drug-card">
        <h4>💊 Lisinopril</h4>
        <b>Target:</b> Angiotensin-converting enzyme (ACE)<br>
        <b>Mechanism:</b> Competitive inhibition<br>
        <b>Disease:</b> Hypertension, Heart Failure<br>
        <b>Impact:</b> First-line therapy for hypertension worldwide
        </div>
        """, unsafe_allow_html=True)

# ==================== PAGE 4: KINETICS ====================
def show_kinetics():
    st.markdown('<div class="section-header">📈 Interactive Kinetics Explorer</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("Kinetic Parameters")
        
        # Michaelis-Menten parameters
        km = st.slider("Michaelis Constant (Km)", 0.1, 10.0, 1.0, 0.1, key="kinetics_km")
        vmax = st.slider("Maximum Velocity (Vmax)", 1.0, 100.0, 50.0, 1.0, key="kinetics_vmax")
        
        # Inhibitor parameters
        inhibitor_present = st.checkbox("Add Inhibitor", value=True)
        if inhibitor_present:
            ki = st.slider("Inhibition Constant (Ki)", 0.1, 5.0, 1.0, 0.1)
            inhibitor_type = st.selectbox("Inhibitor Type", ["Competitive", "Non-competitive", "Uncompetitive"])
            inhibitor_concentration = st.slider("[Inhibitor]", 0.1, 10.0, 1.0, 0.1)
    
    with col2:
        # Generate kinetic data
        substrate_conc = np.linspace(0.1, 20, 100)
        
        # No inhibitor curve
        velocity_no_inhib = vmax * substrate_conc / (km + substrate_conc)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=substrate_conc, y=velocity_no_inhib,
            name='No Inhibitor',
            line=dict(color='blue', width=3)
        ))
        
        if inhibitor_present:
            # Calculate inhibited velocity based on type
            if inhibitor_type == "Competitive":
                alpha = 1 + (inhibitor_concentration / ki)
                velocity_inhib = vmax * substrate_conc / (km * alpha + substrate_conc)
            elif inhibitor_type == "Non-competitive":
                alpha = 1 + (inhibitor_concentration / ki)
                velocity_inhib = (vmax / alpha) * substrate_conc / (km + substrate_conc)
            else:  # Uncompetitive
                alpha = 1 + (inhibitor_concentration / ki)
                velocity_inhib = vmax * substrate_conc / (km + substrate_conc * alpha)
            
            fig.add_trace(go.Scatter(
                x=substrate_conc, y=velocity_inhib,
                name=f'With {inhibitor_type} Inhibitor',
                line=dict(color='red', width=3, dash='dash')
            ))
        
        fig.update_layout(
            title="Michaelis-Menten Kinetics",
            xaxis_title="Substrate Concentration [S] (mM)",
            yaxis_title="Reaction Velocity v (μM/min)",
            height=500
        )
        st.plotly_chart(fig, use_container_width=True)

# ==================== PAGE 5: PIPELINE ====================
def show_pipeline():
    st.markdown('<div class="section-header">⏳ Drug Development Pipeline</div>', unsafe_allow_html=True)
    
    # Interactive pipeline stages
    stages = [
        {"name": "Target Identification", "duration": "1-2 years", "success_rate": "High"},
        {"name": "Lead Discovery", "duration": "1-2 years", "success_rate": "Medium"},
        {"name": "Preclinical Testing", "duration": "2-3 years", "success_rate": "60%"},
        {"name": "Phase I Clinical Trials", "duration": "1-2 years", "success_rate": "50%"},
        {"name": "Phase II Clinical Trials", "duration": "2-3 years", "success_rate": "30%"},
        {"name": "Phase III Clinical Trials", "duration": "3-4 years", "success_rate": "25%"},
        {"name": "Regulatory Approval", "duration": "1-2 years", "success_rate": "85%"},
        {"name": "Post-Market Surveillance", "duration": "Ongoing", "success_rate": "N/A"}
    ]
    
    selected_stage = st.selectbox("Select Development Stage to Explore:", [s["name"] for s in stages])
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        for stage in stages:
            if stage["name"] == selected_stage:
                st.markdown(f"""
                <div class="pipeline-stage">
                <h4>🎯 {stage['name']}</h4>
                <b>Duration:</b> {stage['duration']}<br>
                <b>Success Rate:</b> {stage['success_rate']}<br>
                <b>Key Activities:</b> Target validation, assay development
                </div>
                """, unsafe_allow_html=True)
    
    with col2:
        # Pipeline visualization
        pipeline_data = {
            'Stage': [s["name"] for s in stages],
            'Duration (years)': [float(s["duration"].split('-')[0]) for s in stages],
            'Success Rate (%)': [100 if s["success_rate"] == "High" else 
                               80 if s["success_rate"] == "Medium" else
                               60 if s["success_rate"] == "60%" else
                               50 if s["success_rate"] == "50%" else
                               30 if s["success_rate"] == "30%" else
                               25 if s["success_rate"] == "25%" else
                               85 if s["success_rate"] == "85%" else 0 for s in stages]
        }
        
        fig = px.bar(pipeline_data, x='Stage', y='Success Rate (%)',
                    title="Success Rates Across Development Stages",
                    color='Success Rate (%)')
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)

# ==================== PAGE 6: FUTURE TRENDS ====================
def show_future_trends():
    st.markdown('<div class="section-header">🚀 Future Directions & Emerging Technologies</div>', unsafe_allow_html=True)
    
    trend = st.selectbox(
        "Explore Future Trends:",
        ["AI in Drug Discovery", "PROTAC Technology", "Allosteric Inhibitors", 
         "Personalized Medicine", "Multi-target Inhibitors"]
    )
    
    if trend == "AI in Drug Discovery":
        st.markdown("""
        <div class="info-card">
        <h3>🤖 Artificial Intelligence in Enzyme Inhibitor Discovery</h3>
        </div>
        """, unsafe_allow_html=True)
        
        st.write("""
        **Current Applications:**
        - Predictive modeling of inhibitor-enzyme interactions
        - High-throughput virtual screening
        - De novo drug design using generative AI
        
        **Impact:**
        - Reduced discovery time from years to months
        - Higher success rates in clinical trials
        - Identification of novel binding sites
        """)
        
        # AI adoption timeline
        data = {'Year': [2020, 2022, 2024, 2026, 2028],
               'AI-Adopted Projects (%)': [15, 35, 60, 80, 95]}
        df = pd.DataFrame(data)
        fig = px.line(df, x='Year', y='AI-Adopted Projects (%)', 
                     title="Projected AI Adoption in Drug Discovery")
        st.plotly_chart(fig, use_container_width=True)
    
    elif trend == "PROTAC Technology":
        st.markdown("""
        <div class="info-card">
        <h3>🔗 PROTACs: Proteolysis-Targeting Chimeras</h3>
        </div>
        """, unsafe_allow_html=True)
        
        st.write("""
        **Revolutionary Approach:**
        - Instead of inhibiting, they tag enzymes for destruction
        - Can target "undruggable" enzymes
        - Higher specificity and lower dosing
        
        **Advantages:**
        - Overcome drug resistance
        - Target previously inaccessible enzymes
        - Longer duration of action
        """)

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
        show_kinetics()
    elif selected_section == "Pipeline":
        show_pipeline()
    else:
        show_future_trends()
    
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