# These are the only packages you are allowed import:
import os
import glob as glb
import math
import matplotlib.pyplot as plt

# "pass" indicates an empty block of code,
# remove it when filling in the functions.

def my_name() -> str:
    # replace this string with your real name.
    return "Thi Hanh Nguyen Ly"

#readFile function parse text file into a list
def readFile(file):
    handle = open(file)
    lines=[]
    for line in handle:
        line = line.strip()
        lines.append(line)
    handle.close()
    return lines


def parse_spectrum(filename: str) -> list[float]:
    lines=readFile(filename)
    read=[]
    for line in lines:
        #Spectrum data starts at lines without # mark
        #Storing lines without # in another list
        if not line[0]=='#':
            read.append(line.split())
    spectrum_result=[]
    for element in read:
        #Skipping the first column by skipping element with index 0
        for number in range(1, len(element)):
            spectrum_result.append(float(element[number]))
    return spectrum_result


def wavenumber(spectrum: list[float]) -> list[float]:
    x_value=[]
    #Finding range of wavenumber depending on the spectrum
    distance=math.ceil((4000-400)/len(spectrum))
    end_limit=4000+distance
    for num in range(400, end_limit, distance):
        x_value.append(float(num))
    return x_value #List from 400 to 4000, each number corresponding to each spectrum value
        
    
def absorbance(x: float) -> float:
    absorb=0.0
    #Check if x is in range 0 to 100, if not then set absorb=0.0 to avoid errors
    if x < 100 and x >0:
        absorb=(-1)*math.log10(x/100)
    else:
        absorb=0.0
        #print(f"Transmittance {x} is out of range")
    return absorb


def absorbance_list(spectrum: list[float]) -> list[float]:
    absorb_list=[absorbance(x) for x in spectrum]
    return absorb_list


def plot_transmittance(spectrum: list[float]) -> None:
    x=wavenumber(spectrum)
    y=spectrum
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, color='navy')
    plt.title('Plot of IR spectrum')
    plt.xlabel('Wavenumber (cm$^{-1}$)')
    plt.ylabel('Transmittance')
    
    #Inverting x_values in the plot (running from 4000 to 400)
    plt.gca().invert_xaxis()
    
    #Finding fingerprint region
    x_filtered = [val for val in x if 500 <= val <= 1500] #Filter x values between 500 and 1500
    y_filtered = 100
    plt.fill_between(x_filtered, y_filtered, color='lightblue', alpha=0.3)
    #Adding text indicates it is fingerprint region
    plt.text(1500, 20, 'Fingerprint region', fontsize=12, color='navy')
    
    plt.savefig("PureSpectrumFigure.png")
    plt.show()

def plot_sample(pure: list[float], sample: list[float], 
                sample_name=None) -> None:
    x1=wavenumber(pure)
    x2=wavenumber(sample)
    y1=pure
    y1_abs=absorbance_list(y1)
    y2=sample
    y2_abs=absorbance_list(y2)
    plt.figure(figsize=(8, 6))
    plt.plot(x1, y1, color='navy', label='pure antibiotic')
    plt.plot(x2, y2, color='r', label='sample')
    plt.xlabel('Wavenumber (cm$^{-1}$)')
    plt.ylabel('Transmittance')
    if sample_name:
        plt.title(f"Plot of IR spectrum of pure antibiotic and sample {sample_name}")
    else:
        plt.title(f"Plot of IR spectrum of pure antibiotic and sample")
    plt.text(2700, 7, f"Absorbance correlation: {correlation(y1_abs, y2_abs)}", fontsize=10, color='r')
    plt.legend()
    
    #Inverting x_values in the plot (running from 4000 to 400)
    plt.gca().invert_xaxis()
    plt.savefig(f"CorrelationSpectrum{sample_name}.png")
    plt.show()
    
def center_list(x: list[float]) -> list[float]:
    mean=sum(x)/len(x)
    center=[j-mean for j in x]
    return center


def correlation(x: list[float], y: list[float]) -> float:
    x_cent=center_list(x)
    y_cent=center_list(y)
    
    xy=list(map(lambda x_i,y_i: x_i*y_i, x_cent, y_cent)) #multiply each x_i for y_i
    x_sqr=list(map(lambda x_i: x_i**2, x_cent))
    y_sqr=list(map(lambda y_i: y_i**2, y_cent))
    
    numerator=sum(xy)
    denominator=math.sqrt(sum(x_sqr)*sum(y_sqr))
    return round(numerator/denominator, 3)


if __name__ == "__main__":
    print(f"My name is {my_name()}")
    spectrum=parse_spectrum("data/antibiotic_pure.dat")
    spectrum_abs=absorbance_list(spectrum)
    
    #Retrieve all filenames in data directory:
    files=[filename for filename in os.listdir("data/")]
    del files[0] #removing the pure data filename which is already parsed
    
    #Adjusting filenames in another list:
    strains=[]
    for filename in files:
        string1=filename.split('.')[0]        #this gives 'strainA_sample01'
        string3=string1.replace('strain', '') #this gives 'A_sample01'
        string4=string3.replace('sample', '') #this gives 'A_01'
        strains.append(string4)
        
    #Correlation and transmittance data will then be stored into dictionary corresponding to each sample
    sample_spectrum=[]
    corr=[]
    for i in range(0,6):
        for j in range(i*16, i*16 +16):
            sample=parse_spectrum(f"data/{files[j]}")
            sample_abs=absorbance_list(sample)
            #Correlation values are appended respectively from strainA_sample01 to strainF_sample16
            corr.append(correlation(spectrum_abs,sample_abs))
            #Transmittance values are also appended respectively from strainA_sample01 to strainF_sample16
            sample_spectrum.append(sample)
            
    data_corr={strains[i]: corr[i] for i in range(len(strains))}
    data={strains[i]: sample_spectrum[i] for i in range(len(strains))}
    #print(data_corr)
    
    #Write to text file as a table of correlation
    with open('corr_table.txt', 'w') as output:
        output.write(f"Samples\tCorrelation\n")
        for key, value in data_corr.items():
            output.write(f"{key}\t{value}\n")
    
    #Which samples have the highest correlation?
    max_corr=max(data_corr.values())
    max_strain = [sample for sample, value in data_corr.items() if value == max_corr]
    print(f"The strain {max_strain[0]} sample has the highest absorbance correlation which is {max_corr},"
            " thus, optimal strain for scaling up production.")

    #Plotting
    plot_transmittance(spectrum)
    optimal=data[max_strain[0]]
    plot_sample(spectrum, optimal, sample_name=max_strain[0])
    sample1=data['A_01']
    plot_sample(spectrum, sample1, sample_name='A_01')