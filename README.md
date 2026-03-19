# ECON-DataTranslationChallenge
DTC- Do harsher punishments for drunk driving reduce the probability of future drunk-driving offenses?

Project Overview

This project analyzes whether stricter penalties for drunk driving reduce the likelihood of future offenses. Using an econometric approach, we translate data-driven insights into clear, non-technical findings.

Business Question

Do harsher punishments for drunk driving reduce the probability of future drunk-driving offenses?

Methodology:

Applied a Regression Discontinuity Design (RDD) framework

Exploited a cutoff in Blood Alcohol Content (BAC) to identify causal effects

Conducted robustness checks including bandwidth sensitivity and density tests

Key Findings:

There is evidence of a behavioral response around the legal BAC threshold

Drivers just above the cutoff face stricter penalties, which impacts future behavior

The results suggest that policy enforcement can influence decision-making

Repository Structure
ECON-DataTranslationChallenge/
│
├── .Qmd file/        # Quarto analysis and supporting images
├── Figures/          # Visualizations used in analysis
├── R_Script/         # Data processing and modeling scripts
├── Rendered_doc/     # Final rendered report (docx)
├── Ppt/              # Presentation slides
├── .gitignore        # Ignored files configuration
├── README.md         # Project documentation
Files Description

Quarto File: Contains full analysis and narrative

R Scripts: Data cleaning, modeling, and visualization code

Figures: Graphs supporting results

Rendered Document: Final report output

Presentation: Summary slides for communication

Tools & Technologies

1.R

2.Quarto

3.ggplot2

4.Regression Discontinuity methods

How to Run

1.Open the .qmd file in RStudio

2.Install required R packages

3.Render the document using Quarto

4.Review outputs in the Rendered_doc folder

