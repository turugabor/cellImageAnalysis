import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.title('Image Analysis')
if 'count' not in st.session_state:
    st.session_state.count = 0

increment = st.button('Increment')
if increment:
    st.session_state.count += 1

st.write('Count = ', st.session_state.count)

st.write("app is not fully functional yet")