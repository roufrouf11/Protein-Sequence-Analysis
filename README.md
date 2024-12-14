## Protein Sequence Analysis

## Overview
This project implements a program in Perl to analyze protein sequences from the UniProt database. Using the Kyte-Doolittle hydrophobicity scale, the program predicts transmembrane regions by calculating average hydrophobicity across sliding windows. The predictions are evaluated using metrics such as accuracy and the Matthews Correlation Coefficient (MCC).

## Objectives
1. Extract protein sequences from UniProt database entries.
2. Predict transmembrane regions using hydrophobicity values.
3. Evaluate prediction performance using:
   - **Accuracy**
   - **Matthews Correlation Coefficient (MCC)**

## Features
- **Data Input**:
  - Reads a protein sequence from a `protein.swiss` file.
  - Accepts user inputs for window size (`k`) and sensitivity level.

- **Hydrophobicity Analysis**:
  - Calculates hydrophobicity scores using a sliding window approach.
  - Labels transmembrane regions based on hydrophobicity thresholds.

- **Performance Evaluation**:
  - Calculates true positives (TP), true negatives (TN), false positives (FP), and false negatives (FN).
  - Outputs metrics including accuracy and MCC.

- **Output**:
  - Saves predictions and metrics to `output.csv`.

## Implementation Details
1. **Key Components**:
   - Hash `%hyd`: Stores hydrophobicity values for amino acids.
   - Subroutines:
     - `calculate_average`: Computes the average value of a list.
     - `hydrophobicity`: Calculates window-based hydrophobicity scores and writes them to a CSV file.

2. **Workflow**:
   - Reads the protein sequence from the input file.
   - Accepts user-defined parameters for the sliding window size and sensitivity.
   - Processes the sequence to predict transmembrane regions.
   - Evaluates predictions using TP, TN, FP, and FN counts.
   - Computes accuracy and MCC.
   - Writes results to an output file.

3. **Performance Results**:
   - For window size `k = 11`: Accuracy = 97.4%, MCC = 78.8%
   - For window size `k = 15`: Accuracy = 98.0%, MCC = 79.8%
   - For window size `k = 21`: Accuracy = 98.0%, MCC = 70.3%

   **Note**: Larger window sizes do not always yield better predictions. This can be visualized through generated Excel charts.

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/your-repository.git
   cd your-repository
   ```
2. Ensure Perl is installed on your system.

3. Install required Perl modules:
   ```bash
   cpan install Math::Complex
   ```

## Usage
1. Place the protein sequence file (`protein.swiss`) in the project directory.
2. Run the script:
   ```bash
   perl script_name.pl
   ```
3. Follow on-screen prompts to input:
   - Window size (`k`).
   - Sensitivity.
4. Check the output in `output.csv` for results and metrics.

## File Structure
- `script_name.pl`: Perl script implementing the analysis.
- `protein.swiss`: Input file containing the protein sequence.
- `output.csv`: Output file containing predictions and metrics.

## Results Interpretation
- Predictions include identified transmembrane regions with associated metrics.
- Performance metrics (accuracy and MCC) are calculated to assess the quality of predictions.
