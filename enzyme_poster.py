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
        st.plotly_chart(fig, use_container_width=True)

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
        st.plotly_chart(fig, use_container_width=True)

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
            st.plotly_chart(fig, use_container_width=True)

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
    elif selected_section == "Pipeline":
        st.markdown('<div class="section-header">⏳ Drug Development Pipeline</div>', unsafe_allow_html=True)
        # Add pipeline content
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
