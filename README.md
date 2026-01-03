# Regression Analysis (GDT500)

## Overview
This repository contains a fully reproducible statistical analysis using simulated data (n = 200) to demonstrate linear and logistic regression methods. The project was completed as part of an academic assignment for the **GDT500 Multivariable Analysis** course in the **Doctor of Public Health (DrPH)** programme at **Universiti Sains Malaysia (USM)**.

---

## Course & Assignment Context
- **Course:** GDT500 – Multivariable Analysis  
- **Programme:** Doctor of Public Health (DrPH)  
- **Institution:** Universiti Sains Malaysia (USM)  
- **Instructor:** Professor Kamarul Imran Musa  
- **Assessment Type:** Group assignment (simulated data)

---
## Group Information

- **Group:** Group 2
- **Members:**
  - MOHD FITTRI FAHMI BIN FAUZI
  - MOHD FAUZIE BIN ISMAIL
  - AHMAD SHAHRIL HAFIFI BIN SAAD

## Research Objectives
The objectives of this project are to:
1. Simulate realistic public health datasets suitable for regression analysis.
2. Conduct exploratory data analysis (EDA).
3. Fit and interpret multivariable regression models.
4. Examine interaction effects.
5. Perform model diagnostics and assumption checking.
6. Ensure full analytical reproducibility.

---

## Data Description
All datasets used in this project are **fully simulated** and contain **no real participant data**.

### Linear Regression Analysis
- **Outcome:** Quality of life score (`qol`, continuous)
- **Predictors:**
  - Years working (`years_working`)
  - Physical activity level (`phys_activity`)
  - Obesity status (`obesity`)
- Includes an interaction between physical activity and obesity.

### Logistic Regression Analysis
- **Outcome:** Depression status (`depression`, binary)
- **Predictors:**
  - Years working (`years_working`)
  - Physical activity level (`phys_activity`)
  - Obesity status (`obesity`)
- Includes an interaction between physical activity and obesity.

---
## Analysis Workflow
Each analysis follows the structure specified in the assignment task:

- **Part A:** Data simulation and description  
- **Part B:** Exploratory data analysis (EDA)  
- **Part C:** Regression modelling (main effects and interaction models)  
- **Part D:** Model diagnostics and assumption checking  
- **Part E:** Interpretation and conclusions  

All steps are fully scripted and documented within the Quarto reports.

---
### Data Generation

All datasets are simulated within the Quarto documents themselves:

- **Linear analysis:** Data generation code is in the "Part A" section of `linear.qmd`
- **Logistic analysis:** Data generation code is in the "Part A" section of `logistic.qmd`

The simulation uses fixed seeds (`set.seed()`) to ensure reproducibility.

---

## Published Reports (Posit Connect)
The rendered reports have been published to USM Posit Connect for easy access:

- **Linear Regression Analysis:**  
  https://posit-connect.kk.usm.my/content/dd4f266c-ca4b-45b9-91b3-a7064444874a

- **Logistic Regression Analysis:**  
  https://posit-connect.kk.usm.my/content/6f3c7710-f751-495d-baf4-80ec79ab0cc1

---

## Repository Structure

The repository is organised into **two self-contained analysis folders**, allowing each analysis to be run independently.

```
├── multiple_linear/
│   ├── linear.qmd                      # Quarto analysis document
│   ├── linear.html                     # Rendered HTML report
│   ├── references.bib                   # Bibliography file
│   ├── styles.css                       # HTML styling file
│   ├── generate_data/                   # Data generation subfolder
│   │   ├── linear_generatedata.R        # R script for data simulation
│   │   └── qol_data.csv                 # Generated quality of life dataset
│
├── multiple_logistic/
│   ├── logistic.qmd                     # Quarto analysis document
│   ├── logistic.html                    # Rendered HTML report
│   ├── references.bib                   # Bibliography file
│   ├── styles.css                       # HTML styling file
│   ├── generate_data/                   # Data generation subfolder
│   │   ├── logistic_generatedata.R      # R script for data simulation
│   │   └── data_logistic_regression.csv # Generated depression dataset
│
├── README.md                            # This file
└── LICENSE                              # MIT License
```

---

## Prerequisites

To reproduce this analysis, you will need:

- **R** (version ≥ 4.0.0)
- **RStudio** (recommended for `.Rproj` workflow)
- **Quarto** (version ≥ 1.4)

### Required R Packages

Both analyses use the following R packages:

- **Data manipulation & workflow:**
  - `tidyverse`, `dplyr`, `modelr`, `tibble`
  
- **Visualization:**
  - `ggplot2`, `ggridges`, `GGally`, `patchwork`, `corrplot`, `interactions`
  
- **Statistical modelling & diagnostics:**
  - `broom`, `broom.helpers`, `performance`, `car`, `lmtest`, `mfp`
  
- **Table & report generation:**
  - `gtsummary`, `gt`, `knitr`
  
- **Data exploration & formatting:**
  - `summarytools`, `labelled`, `DT`
  
- **Logistic regression specific:**
  - `caret`, `ResourceSelection`, `pROC`

## How to Reproduce

### Full Reproducibility

Each analysis folder is self-contained and can be reproduced independently:

1. **Navigate to the analysis folder:**
   ```bash
   cd "Assignment_Linear Regression"  # or cd Logistic
   ```

2. **Install required R packages:**
   Open R or RStudio and run:
   ```R
   install.packages(c(
     "tidyverse", "dplyr", "modelr", "tibble",
     "ggplot2", "ggridges", "GGally", "patchwork", "corrplot", "interactions",
     "broom", "broom.helpers", "performance", "car", "lmtest", "mfp",
     "gtsummary", "gt", "knitr",
     "summarytools", "labelled", "DT",
     "caret", "ResourceSelection", "pROC"
   ))
   ```


3. **Render the Quarto document:**
   ```bash
   quarto render linear.qmd  # or logistic.qmd
   ```

4. **View the output:**
   - The rendered HTML report will be generated in the same folder
   - Open `linear.html` or `logistic.html` in your browser


---

## Contributing

This is an academic assignment repository. While the code is open for educational purposes, contributions are limited to the group members listed above.

If you find errors or have suggestions, please:
1. Open an issue describing the problem
2. Reference the specific file and line number
3. Provide context for the suggested improvement

---

## Contact

For questions or clarifications regarding this analysis:

- **Primary Contact:** MOHD FITTRI FAHMI BIN FAUZI
- **Email:** fittri.fahmi@gmail.com
- **Course Instructor:** Professor Kamarul Imran Musa
- **Institution:** Universiti Sains Malaysia (USM)

---

## Acknowledgments

We would like to acknowledge:

- **Professor Kamarul Imran Musa** for supervision and guidance throughout this assignment
- **Department of Community Medicine, USM** for providing the learning environment and resources
- The **R and Quarto communities** for developing the open-source tools that made this reproducible analysis possible
- **Posit, PBC** for maintaining the Posit Connect infrastructure used for report hosting

---
## Disclaimer

This repository contains academic work submitted for the GDT500 Multivariable Analysis course. All datasets are simulated and do not contain real participant data. The analyses are for educational demonstration purposes and should not be used for clinical decision-making.  AI tools were used to assist with code structuring and language refinement where appropriate, with all analytical decisions made by the authors.

---

## Version History

- **v1.0.0 (2024-06-01):** Initial release with complete linear and logistic regression analyses.
---

## License
This repository is licensed under the MIT License.
