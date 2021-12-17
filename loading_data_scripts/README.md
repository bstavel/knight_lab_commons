# Loading iEEG data across different recording system formats


### BCI2000

WashU, Loma Linda, Albany, Mayo Clinic all use this format.

Data collected at different time points can sometimes require different methods. For newer data, the most robust method seems to be using the mex files in the `BCI2000mex` folder within this folder. To use this method:

1. Download the `BCI2000mex` folder to wherever you are loading the data
2. If on a mac, allow not verified developer code by running this in the terminal:        
 `sudo xattr -r -d com.apple.quarantine /path/to/folder/BCI2000mex/*`
3. Open a clean instance of matlab. Make sure you do NOT have fieldtrip in your path
4. In matlab, run the following:
```
addpath(genpath('/path/to/folder/BCI2000mex'))
path_data = '/path/to/data/0000.dat'
[signal,states,parameters] = load_bcidat(path_data);
```

### Irvine

#### Neuralynx

#### Nihan Kohden data 

#### Aligning the two systems
Often, the two recording systems that must be aligned in order to use the data. The first is the Nihan Kohden data, which is the system that records the data used in any clinical decisions and often has more electrodes. Additionally, there are Neuralynx data which are set up by the Irvine team, often have single units, more noise, less channels, and the photodiode channel. In order to use the Nihan Kohden data, you have to do a cross correlation anaylsis to line up the nerual signals in time so that the photodiode channel is usable. [This script](https://github.com/bstavel/knight_lab_commons/blob/main/loading_data_scripts/Irvine_Align_NKT_2_NLX.m) identifies the proper way to lag the Neuralynx (NLX) data so that the photodiode is usable for either dataset.

Be careful with the script and any expectations you may have between the two datasets. Sometimes the +/- of the two datasets are switched, sometimes electrodes are mislabeled, there are breaks in the NLX data, etc. Don't run this straight through and take whatever number it returns. Check each step to make sure it makes sense and please add to it as you see improvements!

The script lives here: `main/loading_data_scripts/Irvine_Align_NKT_2_NLX.m`


