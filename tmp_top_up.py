import os, json, traceback
from tmp_software import SummaryTMP
import pandas as pd

# Getting accessions
accessions = pd.read_csv(
    'uvsx/data/curated_database/top_hits.csv')['accession']

all_df = pd.read_csv('TMP_results/summary/summary_result.csv')

# Create output directory safely
os.makedirs('TMP_results/summary', exist_ok=True)

# Get the completed accessions
done = set(all_df['accession'].astype(str))

# Classifying each accession
for acc in accessions:
    if acc not in done:

        print(f'========== {acc} ==========\n')

        try:

            # Run model
            result = SummaryTMP(acc, "gemma4:e4b")

            print(f'{result.get("thermal_range", None)}\n\n')

            # Convert result to dataframe
            df = pd.DataFrame([result])

            # Serialize complex fields safely
            for col in ["metadata", "timings", "nodes"]:
                if col in df.columns:
                    df[col] = df[col].apply(json.dumps)

        except Exception as e:

            print(f'ERROR processing accession {acc}\n\n')
            print(traceback.format_exc())

            # Dynamically create error row
            if len(all_df.columns) > 0:

                error_row = {
                    col: "ERROR"
                    for col in all_df.columns
                }

            else:
                # Fallback if very first accession fails
                error_row = {}

            # Preserve accession + error message
            error_row["accession"] = acc
            error_row["error_message"] = str(e)

            df = pd.DataFrame([error_row])

        # Append result
        all_df = pd.concat(
            [all_df, df],
            ignore_index=True,
            sort=False
        )

        # Save progress after every accession
        all_df.to_csv(
            'TMP_results/summary/summary_result.csv',
            index=False
        )


############################################################################################################################################################
# Total duration
############################################################################################################################################################

df = pd.read_csv(
    'TMP_results/summary/summary_result.csv'
)

total_duration = pd.to_numeric(
    df['duration'],
    errors='coerce'
).sum()

print(
    f'Total Duration: '
    f'{round((total_duration / 60), 2)} minutes'
)
