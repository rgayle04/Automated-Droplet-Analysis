import os
import math 
import pandas as pd 
import glob
import numpy as np
from numpy.polynomial import polynomial as P 
import scipy.stats as stats
from collections import defaultdict
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import sys
from pathlib import Path



extensions = ['*.csv','*.xl','*.xlsx', '*.xlsm']



def swiffer(csv_path):
    #os.makedirs(opath, exist_ok=True)
    print(f'{csv_path}')

    df = pd.read_csv(csv_path)
    name = Path(csv_path).stem 
    print(f'Original Rows: {len(df)}')

    #print(f'{int(df.isna().sum())} instances of NaN data in {name}')
    print(f'Cleaning {name}...')
    for col in df.columns: 
        before = len(df)
        df = df[~df[col].astype(str).str.contains('-nan\(ind\)|#NAME\?', regex=True, case=False, na=False)]
        after = len(df)
        if before != after:
            print(f'Removed {(before - after)} rows from column {col}')
    #df = df.dropna(inplace=True)

    print(f'Cleaned Rows: {len(df)}')

    file_path= csv_path.replace('.csv', '.xlsx')
    df.to_csv(csv_path, index=False)

    df.to_excel(file_path, index=False)
    

def swapFileExt(csv_path):
    file_path= csv_path.replace('.csv', '.xlsx')

    return file_path

def analyze(csv_path):
    name = Path(csv_path).stem
    df = pd.read_excel(csv_path)
    #os.makedirs(output_path, exist_ok=True)
    if df.empty:
        return
    else:

        std = np.std(df['DIB Radius'], ddof= 1)
        print('DIB Rad Standard Deviation:', std)
        mean = df['DIB Radius'].mean()
        print('DIB Rad Mean:', mean)

        uL = mean + (std*3) 
        print(f'Upper Limit: {uL}') 
        lL = mean - (std*3)
        print(f'Lower Limit: {lL}')
        #df = df[df['DIB Radius'].notna()]

        outliers = np.where((df['DIB Radius'] > uL) | (df['DIB Radius'] < lL) | (df['DIB Radius'].isna()))

        print(f'Number of Outliers: {len(outliers[0])}')
        df.drop(outliers[0],axis=0, inplace=True)
        print(f'Rows remaining after cleaning: {len(df)}')
        split_name = name.split(" ")
        sp_name = split_name[-1]

        osmP = ""

        for i in range(3):
            if len(osmP)<3:
                osmP+=sp_name[i]
        #print(int(osmP))
        
        osmP = float(osmP)/1000

        df['Org. Concentr'] = None

        df.at[0, 'Org. Concentr'] = osmP

        df['Adjusted Time'] = (df['Time Stamp']-df.loc[0, 'Time Stamp'])

        df['DIB Area'] =  math.pi*(df['DIB Radius']**2)
        
        av = df.loc[0, 'Droplet 1 Volume']
        
        df['(V/V0)^2']= (df['Droplet 1 Volume']/av)**2

        df['Linear Permebility Section:'] = None 

        df['Init Radius'] = None
        df.at[0, 'Init Radius'] = (df.loc[0,'Droplet 1 Radius']/10000)

        #f = (df.loc[1,'Droplet 1 Radius']/10000) 

        df['Init Vol'] = None
        df.at[0, 'Init Vol'] = (4/3)*3.1415*(df.loc[0,'Init Radius']**3)

        df['Linearized DIB Radius'] = None
        slope, intercept, r, p, std_error = stats.linregress(df['Adjusted Time'], df['DIB Radius'])

        df.at[0, 'Linearized DIB Radius'] = intercept

        df['DIB Radius(cm)'] = None
        df.at[0, 'DIB Radius(cm)'] = intercept/10000

        df['DIB Area(cm^2)'] = None
        df.at[0, 'DIB Area(cm^2)'] = math.pi*(df.loc[0,'DIB Radius(cm)']**2)

        slope, intercept, r, p, std_error = stats.linregress(df['Adjusted Time'], df['(V/V0)^2'])
        
        df['slope (V/V0)^2'] = None
        df.at[0, 'slope (V/V0)^2'] = slope

        df['r^2'] = None
        df.at[0, 'r^2'] = intercept

        df['Permeability (avg DIB Rad)'] = None 

        df.at[0, 'Permeability (avg DIB Rad)'] = ((slope/2)* df.loc[0, 'DIB Radius(cm)'])/(df.loc[0, 'DIB Area(cm^2)']*0.018*df.loc[0,'Org. Concentr'])*2

        df['3rd Degree Polynomial Section:'] = None

        x = df['Adjusted Time'].dropna().values

        y = df['DIB Area'].dropna().values

        min_len = min(len(x), len(y))

        x = x[:min_len]
        y = y[:min_len]

        coeffs = P.polyfit(x,y,3)


        #X = np.column_stack([x**1, x**2, x**3])

        #X_with_intercept = np.column_stack([np.ones(len(x)), X])

        #coefficients = np.linalg.lstsq(X_with_intercept, y, rcond=None)[0]

        
#        coefficients = np.polyfit(x,y,3)
#        a3, a2, a1, a0 = coefficients

        #a = coefficients.tolist() 
        df['A'] = None
        df.at[0, 'A'] = coeffs[3]
        df['B'] = None
        df.at[0, 'B'] = coeffs[2]
        df['C'] = None
        df.at[0, 'C'] = coeffs[1]
        df['D'] = None
        df.at[0, 'D'] = coeffs[0] 
        df['Eval']=((0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'A']*df['Adjusted Time']**4)/(2*df.loc[0, 'Droplet 1 Volume'])+(2*0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'B']*df['Adjusted Time']**3)/(3*df.loc[0,'Droplet 1 Volume'])+(0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'C']*df['Adjusted Time']**2)/df.loc[0,'Droplet 1 Volume']+(2*0.018*df.loc[0, 'Org. Concentr']*df.loc[0, 'D']*df['Adjusted Time'])/df.loc[0,'Droplet 1 Volume'])

        
        #df.at[0, 'Permeability'] = ((slope/2)* df.loc[0, 'DIB Radius(cm)'])/(df.loc[0, 'DIB Area(cm^2)']*0.018*(df.loc[0, 'Org. Concentr']))

        
        #slope, intercept, r, p ,std_error = stats.linregress(df['Eval'], df['(V/V0)^2'])
        
        x = df['Eval'].values
        y = df['(V/V0)^2'].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        x = x[mask]
        y = y[mask]
        
        slope, intercept = np.polyfit(x,y,1)
        df['Permeability (slope)']= None
        df.at[0, 'Permeability (slope)'] = slope

        df['Permeability (intercept)'] = None
        df.at[0, 'Permeability (intercept)'] = intercept
        #excel_path = csv_path.replace('.csv', '.xlsx')

        df.to_excel(csv_path, index=False)

        print(f'[DONE] {name} processed. \n')
'''
def analyzeFiles(path):
    print(f'{path}')
    os.makedirs(path, exist_ok=True)

    directories = path.split("\\")
    print(directories[3])

    output = os.path.join(path, f"Summarized {directories[3]} .csv")

    csvFiles = []
    for ext in extensions:
        csvFiles.extend(glob.glob(os.path.join(path, ext)))
    csvFiles = sorted(csvFiles)

    # Group files by name (excluding last 3 characters)
    grouped_files = defaultdict(list)
    
    for csv_file in csvFiles:
        name = Path(csv_file).stem
        # Group by name without last 3 characters
        group_key = name[:-3] if len(name) > 3 else name
        grouped_files[group_key].append(csv_file)

    header = True
    if os.path.exists(output) and os.path.getsize(output) > 0:
        header = False
    
    with open(output, 'a', newline='') as outfile:
        if header:
            outfile.write('Video Name,Linear Permeability,3rd Deg Permeability,Linear Std,3rd Deg Poly Std, File Count\n')

        # Process each group
        for group_name, file_list in sorted(grouped_files.items()):
            print(f'Processing group: {group_name} with {len(file_list)} files')
            
            linear_values = []
            poly_values = []
            i = 0
            for csv_file in file_list:
                name = Path(csv_file).stem
                print(f'  - {csv_file}')
                df = pd.read_csv(csv_file)
                linear = df.at[0, 'Permeability (avg DIB Rad)']
                poly = df.at[0, 'Permeability (slope)']
                linear = abs(linear)
                poly = abs(poly)
                linear_values.append(linear)
                poly_values.append(poly)
                outfile.write(f'{name}, {linear},{poly}\n')
                i+=1
                
            if len(file_list) > 1:
            # Calculate mean and standard deviation
                linear_mean = np.mean(linear_values)
                linear_std = np.std(linear_values, ddof= 1)
                poly_mean = np.mean(poly_values)
                poly_std = np.std(poly_values, ddof= 1)
                file_count = len(file_list)
                
                outfile.write(f'{group_name},{linear_mean},{poly_mean},{linear_std},{poly_std},{file_count}\n')
                if i == len(file_list):
                    outfile.write(f'\n')
            else:
                outfile.write(f'\n')
'''        

def main(csv_path):
    if os.path.isdir(csv_path):
        failed = 0
        failedFiles = []
        csvFiles = []
        for ext in extensions:
            csvFiles.extend(glob.glob(os.path.join(csv_path, ext)))
        for file in csvFiles:
            try: 
                    
                swiffer(file)
                newX = swapFileExt(file)
                print(newX)
                analyze(newX)
            except Exception as e:
                failed+=1
                failedFiles.append(file)
                print(f"✗ ERROR processing {file}: {str(e)}")
                print(f"Skipping to next file...\n")
                continue

            
        #analyzeFiles(csv_path)

    else:
        swiffer(csv_path)
        newX = swapFileExt(csv_path)
        print(newX)
        analyze(newX)
    print(f'Failed Files:{failed}')
    print(f'Failed Processing These Files: {failedFiles}')

if __name__=="__main__":
    csv_path = sys.argv[1]
    #opath = sys.argv[2]

    main(csv_path)
    
