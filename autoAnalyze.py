import os
import math 
import pandas as pd 
import glob
import numpy as np
import scipy.stats as stats
from sklearn.linear_model import LinearRegression
import sys
from pathlib import Path

extensions = ['*.csv','*.xl','*.xlsx', '*.xlsm']

def analyze(csv_path, output_path):
    name = Path(csv_path).stem
    df = pd.read_csv(csv_path)
    os.makedirs(output_path, exist_ok=True)
    name = Path(csv_path).stem
    df = pd.read_excel(csv_path)
    #os.makedirs(output_path, exist_ok=True)
    if df.empty:
        return
    else:
        std = np.std(df['DIB Radius'], ddof=1)
        print('DIB Radius Standard Deviation:', std)
        mean = df['DIB Radius'].mean()
        print('DIB Radius Mean:', mean)

        uL = mean + (std*3) 
        print(f'Upper Limit: {uL}') 
        
        lL = mean - (std*3)
        print(f'Lower Limit: {lL}')

        df.dropna(subset=['DIB Radius'], inplace=True)

        outliers = np.where((df['DIB Radius'] > uL) | (df['DIB Radius'] < lL))
        print(f'Rows before Cleaning: {len(df)}')
        print(f'Number of Outliers: {len(outliers[0])}')

        if len(outliers[0]) < len(df) * 0.5:
            df.drop(outliers[0], axis=0, inplace=True)
            df.reset_index(drop=True, inplace=True)
            print(f'Rows remaining: {len(df)}')
        else:
            print("Too many outliers detected - skipping removal")
        
        # ===== EXTRACT osmP WITH ERROR HANDLING =====
        try:
            split_name = name.split(" ")
            sp_name = split_name[-1]

            osmP = ""

            for i in range(3):
                if len(osmP) < 3:
                    osmP += sp_name[i]

            osmP = float(osmP) / 1000
            print(f"Osmolarity extracted: {osmP} mol/L")

        except (ValueError, IndexError) as e:
            print(f"✗ ERROR: Could not extract valid osmolarity from filename: {name}")
            print(f"  Expected format with 3-digit number in filename (e.g., '...170mOsm...')")
            print(f"  Skipping this video.\n")
            cap.release()
            vid.release()
            cv2.destroyAllWindows()
            return
        
        # ===== REST OF PROCESSING =====
        df['Org. Concentr'] = None
        df.at[0, 'Org. Concentr'] = osmP

        df['Adjusted Time'] = (df['Time Stamp'] - df.loc[0, 'Time Stamp'])

        df['DIB Area'] = math.pi * (df['DIB Radius']**2)
        
        av = df.loc[0, 'Droplet 1 Volume']
        
        df['(V/V0)^2'] = (df['Droplet 1 Volume'] / av)**2

        df['Linear Permebility Section:'] = None 

        df['Init Radius'] = None
        df.at[0, 'Init Radius'] = (df.loc[0, 'Droplet 1 Radius'] / 10000)

        df['Init Vol'] = None
        df.at[0, 'Init Vol'] = (4/3) * 3.1415 * (df.loc[0, 'Init Radius']**3)

        df['Linearized DIB Radius'] = None
        slope, intercept, r, p, std_error = stats.linregress(df['Adjusted Time'], df['DIB Radius'])

        df.at[0, 'Linearized DIB Radius'] = intercept

        df['DIB Radius(cm)'] = None
        df.at[0, 'DIB Radius(cm)'] = intercept / 10000

        df['DIB Area(cm^2)'] = None
        df.at[0, 'DIB Area(cm^2)'] = math.pi * (df.loc[0, 'DIB Radius(cm)']**2)

        slope, intercept, r, p, std_error = stats.linregress(df['Adjusted Time'], df['(V/V0)^2'])
        
        summary_cols = {
            'slope (V/V0)^2': slope,
            'r^2': intercept,
            'Permeability (avg DIB Rad)': ((slope/2) * df.loc[0, 'DIB Radius(cm)']) / (df.loc[0, 'DIB Area(cm^2)'] * 0.018 * df.loc[0, 'Org. Concentr']) * 2,
            '3rd Degree Polynomial Section:': None
        }

        for col, val in summary_cols.items():
            if col not in df.columns:
                df[col] = np.nan
            df.at[0, col] = val

        x = df['Adjusted Time'].dropna().values
        y = df['DIB Area'].dropna().values

        min_len = min(len(x), len(y))
        x = x[:min_len]
        y = y[:min_len]

        coeffs = P.polyfit(x, y, 3)

        df['A'] = None
        df.at[0, 'A'] = coeffs[3]
        df['B'] = None
        df.at[0, 'B'] = coeffs[2]
        df['C'] = None
        df.at[0, 'C'] = coeffs[1]
        df['D'] = None
        df.at[0, 'D'] = coeffs[0] 
        
        df['Eval'] = ((0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'A']*df['Adjusted Time']**4)/(2*df.loc[0, 'Droplet 1 Volume'])+(2*0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'B']*df['Adjusted Time']**3)/(3*df.loc[0,'Droplet 1 Volume'])+(0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'C']*df['Adjusted Time']**2)/df.loc[0,'Droplet 1 Volume']+(2*0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'D']*df['Adjusted Time'])/df.loc[0,'Droplet 1 Volume'])

        x = df['Eval'].values
        y = df['(V/V0)^2'].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        x = x[mask]
        y = y[mask]
        
        slope, intercept = np.polyfit(x, y, 1)
        df['Permeability (slope)'] = None
        df.at[0, 'Permeability (slope)'] = slope

        df['Permeability (intercept)'] = None
        df.at[0, 'Permeability (intercept)'] = intercept

        df.to_excel(csv_path, index=False)

        print(f'[DONE] {name} processed. \n')
def main(csv_path, output_path):
    if os.path.isdir(csv_path):
        csvFiles = []
        for ext in extensions:
            csvFiles.extend(glob.glob(os.path.join(csv_path, ext)))
        for csv in csvFiles:
            analyze(csv, output_path)
    else:
        analyze(csv_path, output_path)

if __name__ == "__main__":
    
    csv_path = sys.argv[1]

    output_path = sys.argv[2]


    main(csv_path, output_path)


