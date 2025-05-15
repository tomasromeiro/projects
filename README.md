# Projects

**Overview**: This repository contains projects I have been working on, particularly at University College London (UCL), in the context of the Computational Finance MSc. It contains scripts (Python and C++) and written reports.

Reports and code are named and numbered as per the list below.

## Project List
1. Option Pricing via Finite Difference Methods (C++) - implements various finite difference methods to price options. The program was developed as part of the assessment for the Applied Computational Finance module.
   This project was completed as part of the **Applied Computational Finance** course at UCL.
   The application implements multiple numerical schemes for solving the option pricing problem via backward-marching finite difference methods:
   - **Explicit Scheme**
   - **Fully Implicit Method**
   - **Crank-Nicolson Technique**
   - **Prices Vanilla European / American Options**

2. Numerical simulation of the two-dimensional heat equation (Python) - the primary aim of this assignment is to determine the time at which the temperature at the center of a 2D reaches a given number. A specific initial/boundary condition setup.
   This project was completed as part of the **Techniques of High-Performance Computing** course at UCL and focuses on:
   - **Numerical methods** for time-stepping (explicit vs. implicit)
   - **Code optimization**, such as vectorization or parallelization in Python
   - Use of **GPU acceleration**
  
3. Credit Risk Prediction Reassessed (Python) - a group research project for the Introduction to Machine Learning in Finance at UCL. It critically re-assesses a published credit risk prediction paper using advanced machine learning techniques. In this project we:
   - Applied XGBoost, Random Forests, SVM, and Multi-Layer Perceptron models with thorough hyperparameter tuning.
   - Developed a robust data preprocessing pipeline that avoids data leakage that occurred in the reference article.
   - Critically re-assessed a published credit risk prediction methodology, exposing unrealistic performance metrics and oversampling practices in the original work.
   - Implemented and finely tuned multiple machine learning models, establishing a realistic and credible performance benchmark for credit risk prediction.
   - Employed a comprehensive evaluation framework using multiple performance metrics and visualizations to illustrate model strengths and limitations.

5. Trading with Neural Networks: A Long–Short Approach (Python) - this project investigates the use of machine learning models to forecast daily S&P 500 returns and construct profitable long–short trading strategies. It was done in the context of the Advanced Machine Learning in Finance module at UCL. The methodology is based on a published paper and extended in time horizon and models applied. It's structured as follows:
   - **Models Used**: Logistic Regression, LSTM, CNN, and Wavelet-CNN.
   - **Data Preparation:** Rolling-window approach using 240-day return sequences for feature generation.
   - **Modeling:** Comparison of logistic regression (baseline) with deep learning models: LSTM, CNN, and Wavelet-CNN.
   - **Portfolio Construction:** Daily ranking of stocks by predicted probability to select the top 10 longs and bottom 10 shorts.
   - **Evaluation:** Analysis of classification performance and key financial metrics (annualized returns, Sharpe ratio, drawdowns), with additional breakdowns by industry.
  
6. A Network Perspective on the Global Meat Trade (Python) – an individual project for the Data Science module at UCL. It explores global trade dynamics of meat products using network analysis on international trade data. In this project I:
   - Constructed and pruned a directed trade network using bilateral meat trade flows from the Global Trade Atlas dataset (2005–2022), provided by the university.
   - Applied key network science tools including centrality measures (e.g., PageRank, betweenness), community detection, and maximum spanning trees to uncover influential countries and critical trade relationships.
   - Identified importer-exporter dependencies through a custom threshold-based analysis highlighting systemic vulnerabilities in global meat supply chains.
   - Assessed network robustness via node removal simulations and spectral analysis of the Laplacian matrix to evaluate structural resilience.
   - Built a full Python pipeline using pandas, networkx, matplotlib, and seaborn to clean, analyze, and visualize complex trade interactions at scale.
   
