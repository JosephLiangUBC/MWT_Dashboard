# Project Description

## Project Objective
This 2-part project aims to accomplish two main goals: 
1. streamline the processing and storage of MWT data (currently only experiments on 10s ISI paradigm).
2. create an easy-to-use and accessible data dashboard for MWT data.

# Part 1: Data processing and uploading to server
This part of the project is accomplished by going through the 4 jupyter notebook files names steps 1 through 4. For the most part, required user input is minimal (but easy to perform where it is required), and the notebooks are designed to be intuitive for the non-initiated to follow. For the most part, going through the steps by running each cell one after the other should achieve success. This step is designed to run on local machines that have been set up with Python and Java SDK 8u331. 

## Step 1: Choreography Extraction through 

### Objective

The objective of Step 1 is to extract and organize worm movement sequences from Multi Worm Tracker (MWT) data in preparation for tap response modeling. This code runs a Java-based processing tool (`Chore.jar`, which is included in this Github repository) on raw MWT experiment files to extract behavioral data and generate standard output files in .dat, .trv and .txt formats for downstream analysis. In order for this to work, you'll need to have Java SDK 8u331 (@Gurmehak add link here).

### Data Sources/Inputs

MWT data organized by screens. Folders must follow the naming convention below for the code to process them accurately:

- screen_name
    - gene_allele-1
        - [plate-1].zip
        - [plate-2].zip
        - ...
    - gene_allele-2
        - [plate-1].zip
        - [plate-2].zip
        - ...
    - ...

@Gurmehak add screenshot as an example to the above format here.

### Outputs

This process generates the following files in the same folder (strain > gene_allele > date):
- `.dat` files
- `.trv` files
- `.txt` files

### Processing Workflow
@Gurmehak elaborate on this (can be short but shouldn't be missing since you have the same section for the other steps)

## Step 2: Tap Response Data Extraction 

### Objective

The objective of Step 2 is to extract tap-induced behavioral responses from worm trajectories and normalize them across experiments. This step enables consistent comparison of movement and posture metrics across strains, screens, and alleles.

The objective of Step 2 is to extract tap-induced behavioral responses from worm trajectories and normalize them across experiments. This step enables consistent comparison of movement and posture metrics across strains, screens, and alleles.

### Data Sources / Inputs

- After running the relevant cells, a widget will be displayed requiring user input
- Input is the **screen folder processed in Step 1**, which contains `.trv` files.
- Users must select a screen (e.g., `PD_Screen`, `Neuron_Genes_Screen`, etc.). 

### Outputs

- `tap_output`: Processed CSV files of DataFrames containing extracted tap data.

### Processing Workflow

#### 1. Read folder

- Use the widget to select the folder and corresponding screen name. If an appropriate screen label does not exist, then generate a new label by adding a new string to the list.
@Gurmehak insert screenshot or recording .gif of the widget under operation 

#### 2. Determine Tap Numbers and Tolerances (user input required)

- User-defined inputs:
  - `number_of_taps`: Usually 30.
  - `ISI`: Interstimulus interval (typically 10 seconds).
  - `first_tap`: Frame at which the first tap occurs (commonly 600). This must be verified in `.trv` files.

- The script dynamically computes frame indices for each tap and assigns a frame tolerance around each tap to allow for timing jitter.  
  For example, if the first tap is at frame 600:
  - Tap 1, tolerance: (598, 602)  
  - Tap 2, tolerance: (608, 612)  
  - Tap 3, tolerance: (618, 622)  
  ...
  - Tap 31 is manually adjusted and added where necessary.


#### 3. Process `.trv` Files by Strain

For each strain:
- Extract metadata such as `Date` and `Screen` from file paths.
- Rename fixed-index columns to meaningful names like `time`, `stim_rev`, `no_rev`, `dist`, `dura`, etc.
- Create new columns:
  - `prob` = `stim_rev` / (`stim_rev` + `no_rev`)
  - `speed` = `dist` / `dura`
- Assigns taps based on tolerance windows.

#### 4. Merge All DataFrames

- Concatenate all strain DataFrames.
- Split genotype into `Gene` and `Allele`. If allele is missing (e.g., N2), default to 'N2'.
- Drop missing/NA values.

#### 5. Save as CSV

- Export the merged tap response DataFrame as a CSV file.


## Step 3: Baseline and PSA Data Extraction 

### Objective 

The objective of Step 3 is to extract and summarize both **baseline behavior** (before any stimulus) and **post-stimulus tap responses** for each worm across strains and experiments. The extracted features are used to detect deviations in behavior across genotypes and support downstream modeling and visualization.

### Data Sources / Inputs

- This step takes as input the **screen folder processed in Step 1**, which contains `.dat` files for each experiment.
- Screen selection (e.g., PD_Screen, Neuron_Genes_Screen, etc.).

- Input is the **screen folder processed in Step 1**, which contains `.dat` files.
- Users must select a screen (e.g., `PD_Screen`, `Neuron_Genes_Screen`, etc.).

@Gurmehak insert screenshot or recording .gif of the widget under operation 


### Output

Two `.csv` files are generated:

1. `{Screen}_baseline_output.csv`  
   Contains worm behavior between 490s and 590s (pre-stimulus baseline).

2. `{Screen}_post_stimulus.csv`  
   Aggregates metrics aligned to each tap, summarizing behavior in post-stimulus windows.


### Processing Workflow

#### 1. Read folder 

- Select the screen folder and corresponding screen name.
@Gurmehak insert screenshot or recording .gif of the widget under operation 

#### 2. Define `bin`, Tap numbers and Tolerances (user input required)

- Define time bins (typically 1s intervals from 0–1200s). Important for x-axis and y-axis view range for graphing.
- User-defined inputs:
  - `number_of_taps`: Usually 30.
  - `ISI`: Interstimulus interval (typically 10 seconds).
  - `first_tap`: Frame at which the first tap occurs (commonly 600). This must be verified in `.trv` files.
- The script dynamically computes PSA window which is 7–9.5s after each tap. E.g., if first tap is at 600s:
    - Tap 1, tolerance: (607.0, 609.5)
    - Tap 2, tolerance: (617.0, 619.5)
    - Tap 3, tolerance: (627.0, 629.5)
    ...

#### 3. Process Data by strain (.dat files)
 
For each Strain:
- Extract metadata such as `Plate_id`, `Date`, and `Screen` from file paths.
- Rename columns to standard names.

#### 4. Merge All Strains

- Concatenate all cleaned DataFrames.
- Extract `Gene` and `Allele` from `dataset`.

#### 5. Create Baseline Dataset

- Filter rows where `490.0 <= Time <= 590.0`.
- Drop irrelevant columns.

#### 6. Create Post-Stimulus Dataset

- Filter rows corresponding to post-stimulus time windows i.e., time frames 7s - 9.5s after a tap.
- Group by metadata and tap number (['Experiment', 'Screen', 'Date', 'Plate_id', 'Gene', 'Allele', 'dataset', 'taps'])
- Aggregate metrics using mean:
  - Speed, Bias, Angular Speed, Aspect Ratio, Kink, Curve, Crab
  - Time column is aggregated using "min" instead of "mean"

#### 7. Save as CSV

- Export both baseline and post-stimulus CSVs.


## Step 4: Tstat analysis and Database Export

### Objective

The objective of Step 4 is to prepare and integrate behavioral features into a centralized PostgreSQL database. This includes merging baseline, tap, and PSA data; performing feature engineering; running statistical analyses (t-tests, MSD); and exporting final datasets.

### Data Sources / Inputs

- Baseline CSV
- Tap CSV
- PSA CSV
- Screen selection (e.g., `PD_Screen`, `Neuron_Genes_Screen`, etc.).

### Outputs

- Merged tap + PSA data
- Baseline data (unchanged from Step3)
- T-stat by gene 
- T-stat by allele
- MSD by gene
- MSD by allele
- Summarized PSA metrics

### Processing Workflow

#### 1. Read folder 

- Select the screen folder and corresponding screen name.
- Read baseline, tap, and PSA files.
@Gurmehak insert screenshot or recording .gif of the widget under operation 

### 2. Feature Engineering on Tap Data

- Add new features for `dura`, `prob` and `speed` like:
    - `init`: 1st tap
    - `recov`: 31st tap
    - `final`: mean of taps 28–30
    - `habit`: `init` - `final`
    - `recovery`: 100 * (`recov`- `init`) / `init`
    - `memory_retention_dura`: `recov`- `final`
    
### 3. Summarize PSA Data
- Reduce 31-tap PSA data per experiment to one row using aggregation.
- Group by metadata (['Experiment', 'Plate_id', 'Date', 'Screen', 'dataset', 'Gene', 'Allele'])
- Add new features for [`PSA Instantaneous Speed`, `PSA Bias`, `PSA Angular Speed`, `PSA Aspect Ratio`, `PSA Kink`, `PSA Curve`, `PSA Crab`]:
    - `initial`: 1st tap
    - `recovery`: 31st tap
    - `final`: mean of taps 28–30
    - `peak`: maximum value of the specified metric
    - `peak_tap`: tap number where the max value is observed
    - `mean`: average of the values across all 31 taps 
    - `sensitization`: `peak` - `initial`
    - `habituation`: `peak` - `final` 
    - `spontaneous_recovery`:
        - 100 * (`initial` - `recovery`)/`initial` if metric is `PSA Instantaneous Speed`, `PSA Bias`, `PSA Angular Speed`
        - 100 * (`recovery` - `initial`)/`initial` if metric is `PSA Aspect Ratio`, `PSA Kink`, `PSA Curve`, `PSA Crab`
    - `memory_retention`:
        - `final` - `recovery` if metric is `PSA Instantaneous Speed`, `PSA Bias`, `PSA Angular Speed`
        - `recovery` - `final` if metric is `PSA Aspect Ratio`, `PSA Kink`, `PSA Curve`, `PSA Crab`
   
#### 4. Aggregate by Metadata

- Create grouped DataFrames by `gene` and `allele`, conditioned on "baseline", "tap", or "psa".
- Group by metadata ['Plate_id','Date','Screen','dataset','Gene','Allele']

### 5. Mean Sample Distance (MSD) Calculation

- Each phenotype is grouped by Gene or allele (`dataset` column).
- Compute mean, SEM, and 95% CI for each metric by gene/allele.
- Normalize against N2 control.
- Merge across baseline, tap, and PSA.

### 6. T-Test Analysis

- Two-sample t-tests between each strain and N2.
- Run for all metrics in baseline, tap, and PSA datasets.
- Store results separately for genes and alleles.

#### 7. Merge & Final Export

- Combine t-test results across baseline + tap + PSA.
- Clean, rename, reorder and format final outputs

#### 8. Database Export

- Write structured datasets to PostgreSQL.
- Requires local `database.ini` for credentials.

# Part 2: Data Dashboard
